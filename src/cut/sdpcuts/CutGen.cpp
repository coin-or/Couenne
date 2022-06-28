/*
 *
 * Name:    CutGen.cpp
 * Authors: Andrea Qualizza
 *          Pietro Belotti
 * Purpose: Generation of all cuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneExprVar.hpp"
#include "CouenneExprClone.hpp"
#include "operators/CouenneExprMul.hpp"
#include "CouenneProblem.hpp"
#include "CouenneMatrix.hpp"
#include "CouenneSdpCuts.hpp"

#include "dsyevx_wrapper.hpp"

//#define DEBUG

const bool WISE_SPARSIFY = true;

#define SPARSIFY_MAX_CARD  10000
#define SPARSIFY_NEW_NZ_THRESHOLD 0.70

#define EV_TOL 1e-13

using namespace Couenne;


// sdpcut separator -- all minors
void CouenneSdpCuts::generateCuts (const OsiSolverInterface &si, OsiCuts &cs,
				   const CglTreeInfo info)
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
  const
#endif
  {

  // only at root node once

  if ((info . level + info . pass > 4)) return;

  problem_ -> domain () -> push (&si, &cs);

  for (std::vector <CouenneExprMatrix *>::const_iterator
	 minor  = minors_. begin ();
       minor   != minors_. end   (); ++minor)

    genCutSingle (*minor, si, cs, info);

  problem_ -> domain () -> pop ();
}


// sdpcut separator -- one minor at a time
void CouenneSdpCuts::genCutSingle (CouenneExprMatrix * const & minor,
				   const OsiSolverInterface &si,
				   OsiCuts &cs,
				   const CglTreeInfo info) const {
#ifdef DEBUG
  printf ("Generating cut on minor -> ");
  minor -> print ();
#endif

  std::vector <expression *> &varInd = minor -> varIndices ();

  int
    n = (int) (minor -> size ()),
    m,
    nVecs = (numEigVec_ < 0) ? n : numEigVec_,

    **indA  = new int * [n], // contains indices of matrix A below
    *indMap = new int [problem_ -> nVars ()],

    compressed_index = 0;

  double
    *A = new double [n * n];

  // Fill in inverse map
  for (std::vector <expression *>::const_iterator
	 i  = varInd . begin ();
       i   != varInd . end   (); ++i)
    indMap [(*i) -> Index ()] = compressed_index++;

  // ---------------------------------------------------------------------

  for (int i=0; i<n; ++i)
    CoinFillN (indA [i] = new int [n], n, -2);

  // build index and value matrices

  // Get minors_

  CoinFillN (A, n * n, 0.);

#ifdef DEBUG
  printf ("Components\n\n");
#endif

  int indI = 0;

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::const_iterator
	 i  = minor -> getRows () . begin ();
       i   != minor -> getRows () . end   (); ++i, ++indI) {

#ifdef DEBUG
    printf ("row %d [var %d]: ", indI, i -> first); fflush (stdout);
#endif

    int
      majInd = (i -> first == problem_ -> nVars ()) ? n-1 : indMap [i -> first],
      indJ   = 0;

    for (std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::const_iterator
	   j  = i -> second -> getElements () . begin ();
	 j   != i -> second -> getElements () . end   (); ++j, ++indJ) {

#ifdef DEBUG
      printf ("[r%d,v%d,%d,%g,", indJ,
	      (*j) -> getIndex (),
	      (*j) -> getElem () -> Index (),
	      (*((*j) -> getElem ()))  ());
      (*j) -> getElem () -> print (); printf ("] ");
      fflush (stdout);
#endif

      int minInd = ((*j) -> getIndex () == problem_ -> nVars ()) ? (n-1) : (indMap [(*j) -> getIndex ()]);

      expression *Elem = (*j) -> getElem ();

      A [majInd * n + minInd] = (*Elem) (); // evaluate variable at this position of matrix
      //A [indJ * n + indI] = (*Elem) ();

      indA [majInd] [minInd] = Elem -> Index ();
      //indA [indJ] [indI] =
    }

#ifdef DEBUG
    printf ("\n");
#endif
  }

  delete [] indMap;

  // fill in non-existing auxiliaries in X (if any) with their hypothetical product

  for   (int i=0; i<n-1; ++i) // get to n-1 since varIndices not defined afterward
    for (int j=i; j<n-1; ++j)
      if (indA [i] [j] == -2) // never happens on border row/column
	A [i * n + j] =
	A [j * n + i] =
	  (*(varInd [i])) () *
	  (*(varInd [j])) ();

#ifdef DEBUG
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      printf ("[%4d,%7.2g] ", indA [i][j], A [i * n + j]);
    printf ("\n");
  }
#endif

  double
    *Acopy = useSparsity_ ? CoinCopyOfArray (A, n * n) : NULL,
    *w = NULL,
    *z = NULL;

  //  printf ("calling dsyevx NONSPARSE\n");

  //------------------------------------------------------------------------------------------------------
  dsyevx_interface (n, A, m, w, z, EV_TOL,
		    -COIN_DBL_MAX, onlyNegEV_ ? 0. : COIN_DBL_MAX,
		    1, numEigVec_ < 0 ? n : numEigVec_);
  //------------------------------------------------------------------------------------------------------

  if (m < nVecs)
    nVecs = m;

  double
    **work_ev = new double * [m];

  for (int i=0; i < nVecs; i++) {

    if (w [i] >= 0) {
      nVecs = i; // updated nVecs
      break;
    }

    work_ev [i] = z + (i*n);

#ifdef SCALE_EIGENV
    double scaling_factor = sqrt (n);

    for (int j=0; j<n; j++)
      work_ev [i] [j] *= scaling_factor;
#endif
  }

  // Generate sdp (possibly dense) cuts //////////////////////////////////////////////////////////

  for (int i=0; i<nVecs; i++)
    genSDPcut (si, cs, minor, work_ev [i], work_ev [i], indA);

  int
    wise_evdec_num    = 0,
    card_sparse_v_mat = 0,
    min_nz;

  double **sparse_v_mat = NULL;

  if (useSparsity_) {

    sparse_v_mat = new double*[SPARSIFY_MAX_CARD];
    for (int i=0; i<SPARSIFY_MAX_CARD; i++)
      sparse_v_mat[i] = new double [n];

    min_nz = ceil (n * SPARSIFY_NEW_NZ_THRESHOLD);
    card_sparse_v_mat = 0;

    sparsify2 (n, A, sparse_v_mat, &card_sparse_v_mat, min_nz, &wise_evdec_num);

    for (int k=0; k < card_sparse_v_mat; k++)
      genSDPcut (si, cs, minor, sparse_v_mat [k], sparse_v_mat [k], indA);

    ////////////////////////////////////////////

    for (int i=0; i<nVecs; ++i) {

      card_sparse_v_mat = 0;
      double *v = work_ev[i];

      sparsify (WISE_SPARSIFY, i, w [i], v, n, Acopy, sparse_v_mat, &card_sparse_v_mat, &wise_evdec_num);

      for (int k=0; k < card_sparse_v_mat; k++) {

	genSDPcut (si, cs, minor, sparse_v_mat [k], sparse_v_mat [k], indA);

	if (useSparsity_)
	  additionalSDPcuts (si, cs, minor, n, Acopy, sparse_v_mat[k], indA);
      }
    }
  }

  for (int i=0; i<n; ++i)
    delete [] indA [i];
  delete [] indA;

  if (useSparsity_) {

    for (int i=0; i < SPARSIFY_MAX_CARD; i++)
      delete [] sparse_v_mat[i];

    delete [] sparse_v_mat;
    delete [] Acopy;
  }

  delete [] z;
  delete [] w;
  delete [] A;
  delete [] work_ev;
}


/************************************************************************/
void CouenneSdpCuts::genSDPcut (const OsiSolverInterface &si,
				OsiCuts &cs,
				CouenneExprMatrix *XX,
				double *v1, double *v2,
				int **indA) const {
  int
    nterms   = 0,
    n        = (int) (XX -> size ()),
    nvars    = problem_ -> nVars (),
    N        = n * n,
    *ind     = new int [N],
    *inverse = new int [nvars];

  double
    *coeff = new double [N],
    *xtraC = new double [nvars],
    rhs    = 0.;

  std::vector <expression *> &varIndices = XX -> varIndices ();

  CoinFillN (xtraC,   nvars, 0.);
  CoinFillN (inverse, nvars, -1);

  // for (int i=0; i<n; ++i) {
  //   if (fabs (v1 [i]) < COUENNE_EPS) v1 [i] = 0;
  //   if (fabs (v2 [i]) < COUENNE_EPS) v2 [i] = 0;
  // }

#ifdef DEBUG
  printf ("Solution: (");
  for (int i=0; i<n; i++) {
    if (i) printf (",");
    printf ("%g", v1 [i]);
  }
  printf (")\n");
#endif

  // ASSUMPTION: matrix is symmetric

  bool numerics_flag = false;

  // coefficients for X_ij
  for (int i=0; (i<n) && !numerics_flag; i++)

    for (int j=i; j<n; j++) {

      double coeff0 = v1 [i] * v2 [j] + v1 [j] * v2 [i];

      if (fabs (coeff0) < 1e-21) continue; // why 1e-21? Cbc/Clp related

      int index = indA [i] [j];

#ifdef DEBUG
      printf ("{%d,%g} ", index, coeff0);
#endif

      // Three cases:
      //
      // 1) index >= 0: corresponds to variable of the problem,
      //    proceed as normal
      //
      // 2) index == -1: constant, subtract it from right-hand side
      //
      // 3) index < -1: this variable (x_ij = x_i * x_j) is not in the
      //    problem, and must be replaced somehow.

#ifdef DEBUG
      if (index == -1) printf ("found constant: %g\n", (i==j) ? coeff0 / 2 : coeff0);
#endif
      if (index == -1)

	rhs -= ((i==j) ? coeff0 / 2 : coeff0);

      else if (index < -1) { // only happens in non-bordered matrix

	expression
	  *Xi = XX -> varIndices () [i],
	  *Xj = XX -> varIndices () [j];

	double
	  valXi = (*Xi) (),
	  valXj = (*Xj) (),
	  li, lj,
	  ui, uj,
	  L, U;

	Xi -> getBounds (li, ui);
	Xj -> getBounds (lj, uj);

#ifdef DEBUG
	printf ("expression: x%d [%g,%g] * x%d [%g,%g]\n",
		Xi -> Index (), li, ui,
		Xj -> Index (), lj, uj);
#endif

	if (i==j) {

	  // This variable (x_ii = x_i ^ 2) is not in the
	  // problem. Replace it with its upper envelope
	  //
	  // X_ii <= li^2 + (ui^2 - li^2) / (ui-li) * (xi-li) = (ui+li) * xi - li*ui

	  if ((fabs (li) > COUENNE_INFINITY) ||
	      (fabs (ui) > COUENNE_INFINITY)) {

	    // term would be X_ii <= inf, which adds inf to the
	    // inequality

	    numerics_flag = true;
	    break;
	  }

	  xtraC [varIndices [i] -> Index ()] += coeff0 / 2 * (li + ui);
#ifdef DEBUG
	  printf ("adding %g=%g*(%g+%g) (sq) to xtrac[%d=varInd[%d]]\n",
		  coeff0 / 2 * (li + ui), coeff0, li, ui, varIndices [i] -> Index (), i);
#endif
	  rhs += coeff0 / 2 * (li * ui);

	} else {

	  // This variable (x_ij = x_i * x_j) is not in the problem,
	  // and must be replaced somehow. Apply Fourier-Motzkin
	  // elimination using bounds and McCormick inequalities on
	  // the fictitious variable y_ij = x_i x_j:
	  //
	  //    y_ij >= L = min {l_i*l_j, u_i*u_j, l_i*u_j, u_i*l_j}
	  //    y_ij >= l_j x_i + l_i x_j - l_i l_j
	  //    y_ij >= u_j x_i + u_i x_j - u_i u_j
	  //
	  //    y_ij <= U = max {l_i*l_j, u_i*u_j, l_i*u_j, u_i*l_j}
	  //    y_ij <= u_j x_i + l_i x_j - u_i l_j
	  //    y_ij <= l_j x_i + u_i x_j - l_i u_j

	  exprMul Xij (new exprClone (Xi),
		       new exprClone (Xj));

	  Xij . getBounds (L, U);

	  double
	    rhsMll = lj * valXi + li * valXj - lj * li,
	    rhsMuu = uj * valXi + ui * valXj - uj * ui;

	  if (coeff0 < 0) {

	    if (L >= CoinMax (rhsMll, rhsMuu))

	      rhs -= coeff0 * L; // L is the tightest (greatest)

	    else if (rhsMuu > rhsMll) {  // rhsMuu is the tightest

	      if ((fabs (ui) > COUENNE_INFINITY) ||
		  (fabs (uj) > COUENNE_INFINITY)) {
		numerics_flag = true;
		break;
	      }

	      xtraC [varIndices [i] -> Index ()] += coeff0 * uj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * ui;
	      rhs += coeff0 * uj * ui;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMuu) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n",
		      coeff0 * uj, coeff0 * ui, coeff0, uj, ui, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif
	    } else { // rhsMll is the tightest

	      if ((fabs (li) > COUENNE_INFINITY) ||
		  (fabs (lj) > COUENNE_INFINITY)) {

		numerics_flag = true;
		break;
	      }

	      xtraC [varIndices [i] -> Index ()] += coeff0 * lj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * li;
	      rhs += coeff0 * lj * li;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMll) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n",
		      coeff0 * lj, coeff0 * li, coeff0, lj, li, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif
	    }

	  } else { // coeff0 > 0

	    double
	      rhsMlu = lj * valXi + ui * valXj - lj * ui,
	      rhsMul = uj * valXi + li * valXj - uj * li;

	    if (U <= CoinMin (rhsMlu, rhsMul))

	      rhs -= coeff0 * U; // U is the tightest (smallest)

	    else if (rhsMul < rhsMlu) { // rhsMul is the tightest

	      if ((fabs (li) > COUENNE_INFINITY) ||
		  (fabs (uj) > COUENNE_INFINITY)) {

		numerics_flag = true;
		break;
	      }

	      xtraC [varIndices [i] -> Index ()] += coeff0 * uj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * li;
	      rhs += coeff0 * uj * li;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMul) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n",
		      coeff0 * uj, coeff0 * li, coeff0, uj, li, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif

	    } else { // rhsMlu is the tightest

	      if ((fabs (ui) > COUENNE_INFINITY) ||
		  (fabs (lj) > COUENNE_INFINITY)) {

		numerics_flag = true;
		break;
	      }

	      xtraC [varIndices [i] -> Index ()] += coeff0 * lj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * ui;
	      rhs += coeff0 * lj * ui;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMlu) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n",
		      coeff0 * lj, coeff0 * ui, coeff0, lj, ui, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif
	    }
	  }
	}

	///////////////////////////////////////////////////////////

      } else {

#ifdef DEBUG
	printf ("normal term: %g x_%d [terms:%d]\n", (i==j) ? (0.5 * coeff0) : (coeff0), index, nterms);
#endif

	if      (inverse [index] >= 0)
	  coeff [inverse [index]] += (i==j) ? (0.5 * coeff0) : (coeff0);
	else {

	  coeff   [nterms]   = (i==j) ? (0.5 * coeff0) : (coeff0);
	  // if (inverse [index] >= 0)
	  //   printf ("duplicate index at nterms=%d: inverse[%d] = %d, value %g\n",
	  // 	  nterms, index, inverse [index], coeff [nterms]);
	  inverse [index]    = nterms;
	  ind     [nterms++] = index;
	}
      }
#ifdef DEBUG
      printf ("%d terms so far\n", nterms);
#endif
    }

  if (!numerics_flag)
    for (std::vector <expression *>::iterator
	 i  = varIndices . begin ();
       i   != varIndices . end   (); ++i) {

    int varInd = (*i) -> Index ();

    if ((inverse [varInd] >= 0) &&
	(fabs (xtraC [varInd]) > 1e-15)) {
#ifdef DEBUG
      printf ("now adding %g to coeff [inv [%d] = %d]\n", xtraC [varInd], varInd, inverse [varInd]);
#endif
      coeff [inverse [varInd]] += xtraC [varInd];
    } else if (fabs (xtraC [varInd]) > COUENNE_EPS) { // independent variable not appearing so far

      coeff   [nterms]   = xtraC [varInd];
      inverse [varInd]   = nterms; // probably useless
      ind     [nterms++] = varInd;
    }
  }

  delete [] inverse;
  delete [] xtraC;

  if (!numerics_flag) {

    OsiRowCut *cut = new OsiRowCut;
    cut -> setRow (nterms, ind, coeff);
    cut -> setLb (rhs);

    if (nterms > 0) {

#ifdef DEBUG
      printf ("SDP: separating ");
      cut -> print ();
#endif

      CoinAbsFltEq treatAsSame (COUENNE_EPS);
      cs.insertIfNotDuplicate (*cut, treatAsSame);

      double violation = 0.;

      if (problem_ -> bestSol () && ((violation = cut -> violated (problem_ -> bestSol ())) > 0.)) {

	printf ("Cut violates optimal solution by %g\n", violation);
	cut -> print ();
      }
    }

    delete cut;
  }

  delete [] ind;
  delete [] coeff;
}
