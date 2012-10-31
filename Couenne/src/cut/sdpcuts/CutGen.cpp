/* $Id$
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

#define DEBUG

//#define ONLY_NEG_EIGENV
//#define ONLY_MOST_NEG
//#define SPARSIFY

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

#define SPARSIFY_MAX_CARD  10000
#define WISE_SPARSIFY_GAP  0.0001

#define SPARSIFY_OLD_DELTA 0.50
#define SPARSIFY_NEW_DELTA 0.50

#define SPARSIFY_OLD_NZ_THRESHOLD 0.50
#define SPARSIFY_NEW_NZ_THRESHOLD 0.70

static int decomposition_counter;

#define EV_TOL 1e-13

using namespace Couenne;

/************************************************************************/
// sdpcut separator
void CouenneSdpCuts::generateCuts (const OsiSolverInterface &si, OsiCuts &cs, 
				   const CglTreeInfo info) const {

  problem_ -> domain () -> push (&si, &cs);

  for (std::vector <CouenneSparseMatrix *>::const_iterator 
	 minor  = minors_. begin ();
       minor   != minors_. end   (); ++minor)

    genCutSingle (*minor, si, cs, info);

  problem_ -> domain () -> pop ();
}

// sdpcut separator
void CouenneSdpCuts::genCutSingle (CouenneSparseMatrix * const & minor, 
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
    np = n,
    m,
    origsdpcuts_card = 0,
    duplicate_cuts   = 0,

    **indA  = new int * [np], // contains indices of matrix A below
    *indMap = new int [problem_ -> nVars ()];

  double 
    *A = new double [np * np];

  int compressed_index = 0;

  for (std::vector <expression *>::const_iterator 
	 i  = varInd . begin ();
       i   != varInd . end   (); ++i)
    indMap [(*i) -> Index ()] = compressed_index++;

  // ---------------------------------------------------------------------

  for (int i=0; i<np; ++i)
    CoinFillN (indA [i] = new int [np], np, -2);

  // build index and value matrices

  // Get minors_ 

  CoinFillN (A, np * np, 0.);

#ifdef DEBUG
  printf ("Components\n\n");
#endif

  int indI = 0;

  for (std::set <std::pair <int, CouenneSparseVector *> >::iterator 
	 i  = minor -> getRows () . begin ();
       i   != minor -> getRows () . end   (); ++i, ++indI) {

#ifdef DEBUG
    printf ("row %d [var %d]: ", indI, i -> first); fflush (stdout);
#endif

    int
      majInd = (i -> first == problem_ -> nVars ()) ? np-1 : indMap [i -> first],
      indJ   = 0;

    for (std::set <CouenneScalar *>::iterator 
	   j  = (*i) . second -> getElements () . begin ();
	 j   != (*i) . second -> getElements () . end   (); ++j, ++indJ) {

#ifdef DEBUG
      printf ("[r%d,v%d,%d,%g,", indJ, 
	      (*j) -> getIndex (), 
	      (*j) -> getElem () -> Index (), 
	      (*((*j) -> getElem ()))  ());
      (*j) -> getElem () -> print (); printf ("] ");
      fflush (stdout);
#endif

      int minInd = ((*j) -> getIndex () == problem_ -> nVars ()) ? (np-1) : (indMap [(*j) -> getIndex ()]);

      expression *Elem = (*j) -> getElem ();

      A [majInd * np + minInd] = (*Elem) ();
      //A [indJ * np + indI] = (*Elem) ();

      indA [majInd] [minInd] = Elem -> Index ();
      //indA [indJ] [indI] = 
    }

#ifdef DEBUG
    printf ("\n");
#endif
  }

  delete [] indMap;

  for   (int i=0; i<np-1; ++i) // get to np-1 since varIndices not defined afterward
    for (int j=0; j<np-1; ++j)
      if (indA [i] [j] == -2) // never happens on border row/column
	A [i * np + j] = 
	  //A [j * np + i] = 
	  (*(varInd [i])) () * 
	  (*(varInd [j])) ();

#ifdef DEBUG
  for (int i=0; i<np; ++i) {
    for (int j=0; j<np; ++j)
      printf ("[%4d,%7.2g] ", indA [i][j], A [i * np + j]);
    printf ("\n");
  }
#endif

#ifdef SPARSIFY_MINOR_SDP_CUTS
  double *Acopy = CoinCopyOfArray (A, np * np);
#endif

  double *w = NULL, *z = NULL;

  //------------------------------------------------------------------
#if defined ONLY_NEG_EIGENV && defined ONLY_MOST_NEG
  dsyevx_interface (np, A, m, w, z, EV_TOL, -COIN_DBL_MAX,           0., 1, 1);
#elif defined ONLY_NEG_EIGENV
  dsyevx_interface (np, A, m, w, z, EV_TOL, -COIN_DBL_MAX,           0., 1, np);
#elif defined ONLY_MOST_NEG
  dsyevx_interface (np, A, m, w, z, EV_TOL, -COIN_DBL_MAX, COIN_DBL_MAX, 1, 1);
#else
  dsyevx_interface (np, A, m, w, z, EV_TOL, -COIN_DBL_MAX, COIN_DBL_MAX, 1, np);
#endif
  //------------------------------------------------------------------

  double
    **work_ev = new double * [m];

  for (int i=0; i<m; i++) {

#ifdef ONLY_NEG_EIGENV
    if (w [i] > 0)
      break;
#endif

#ifdef ONLY_MOST_NEG
    if(i > 0)
      break;
#endif

    work_ev [i] = z + (i*np);

#ifdef SCALE_EIGENV
    double scaling_factor = sqrt (np);

    for (int j=0;j<np;j++)
      work_ev [i] [j] *= scaling_factor;
#endif
  }

  // Generate sdp cuts proper //////////////////////////////////////////////////////////

  for (int i=0;i<m;i++) {

#ifdef ONLY_NEG_EIGENV
    if (w [i] > 0)
      break;
#endif

#ifdef ONLY_MOST_NEG
    if(i > 0)
      break;
#endif

    genSDPcut (si, cs, minor, work_ev [i], work_ev [i], indA, removeduplicates_, &duplicate_cuts);
    origsdpcuts_card++;
  }

#if (defined SPARSIFY) || (defined SPARSIFY2) || (defined WISE_SPARSIFY)

  int wise_evdec_num = 0;

  int card_sparse_v_mat = 0;
  double **sparse_v_mat = new double*[SPARSIFY_MAX_CARD];
  for (int i=0; i<SPARSIFY_MAX_CARD; i++)
    sparse_v_mat[i] = new double [np];
#endif

#ifdef SPARSIFY2
  int min_nz;

  min_nz = ceil (np * SPARSIFY_NEW_NZ_THRESHOLD);
  card_sparse_v_mat = 0;

  sparsify2 (n,X,sparse_v_mat,&card_sparse_v_mat,min_nz,&wise_evdec_num);

  for (int k=0; k<card_sparse_v_mat; k++)
    genSDPcut (si, cs, minor, sparse_v_mat [k], sparse_v_mat [k], indA, removeduplicates_, &duplicate_cuts);
#endif

  for (int i=0;i<m;i++) {

#ifdef ONLY_NEG_EIGENV
    if (w [i] > 0) break;
#endif

#ifdef ONLY_MOST_NEG
    if(i > 0) break;
#endif

#if  (defined SPARSIFY) || (defined WISE_SPARSIFY)

    card_sparse_v_mat = 0;
    double *v = work_ev[i];

#if (defined WISE_SPARSIFY)
    sparsify (true,  i, w[i], v, n, X, sparse_v_mat, &card_sparse_v_mat, &wise_evdec_num);
#else
    sparsify (false, i, w[i], v, n, X, sparse_v_mat, &card_sparse_v_mat, &wise_evdec_num);
#endif

    for (int k=0; k<card_sparse_v_mat; k++) {

      genSDPcut (si, cs, minor, sparse_v_mat [k], sparse_v_mat [k], indA, removeduplicates_, &duplicate_cuts);

#ifdef SPARSIFY_MINOR_SDP_CUTS
      additionalSDPcuts (si,cs, np, Acopy, sparse_v_mat[k], indA, &duplicate_cuts);
#endif
    }
#endif // (defined SPARSIFY) || (defined WISE_SPARSIFY)
  }

  for (int i=0; i<np; ++i)
    delete [] indA [i];
  delete [] indA;

#if (defined SPARSIFY) || (defined SPARSIFY2) || (defined WISE_SPARSIFY)
  for(int i=0;i<SPARSIFY_MAX_CARD;i++)
    delete [] sparse_v_mat[i];
  delete [] sparse_v_mat;
#endif

  delete [] z;
  delete [] w;
  delete [] A;
  delete [] work_ev;

#ifdef SPARSIFY_MINOR_SDP_CUTS
  delete [] Acopy;
#endif
}


/************************************************************************/
void CouenneSdpCuts::genSDPcut (const OsiSolverInterface &si,
				OsiCuts &cs, 
				CouenneSparseMatrix *XX,
				double *v1, double *v2, 
				int **indA,
				bool checkduplicates, 
				int *duplicate_cuts) const {
  int
    nterms   = 0,
    n        = (int) (XX -> size ()),
    np       = n,
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

#ifdef DEBUG
  printf ("Solution: (");
  for (int i=0; i<np; i++) {
    if (i) printf (",");
    printf ("%g", v1[i]);
  }
  printf (")\n");
#endif

  // ASSUMPTION: matrix is symmetric

  // coefficients for X_ij
  for (int i=0; i<np; i++)

    for (int j=i; j<np; j++) {

      double coeff0 = v1 [i] * v2 [j] + v1 [j] * v2 [i];

      if (fabs (coeff0) < 1e-12) continue; // why 1e-21? Cbc/Clp related

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
      if (index==-1) printf ("found constant: %g\n", (i==j) ? coeff0 / 2 : coeff0);
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

	      xtraC [varIndices [i] -> Index ()] += coeff0 * uj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * ui;
	      rhs += coeff0 * uj * ui;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMuu) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n", 
		      coeff0 * uj, coeff0 * ui, coeff0, uj, ui, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif
	    } else { // rhsMll is the tightest

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

	      xtraC [varIndices [i] -> Index ()] += coeff0 * uj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * li;
	      rhs += coeff0 * uj * li;

#ifdef DEBUG
	      printf ("adding (%g,%g) = %g*(%g,%g) (rhsMul) to xtrac[%d=varInd[%d]] and xtrac[%d=varInd[%d]]\n", 
		      coeff0 * uj, coeff0 * li, coeff0, uj, li, varIndices [i] -> Index (), i, varIndices [j] -> Index (), j);
#endif

	    } else { // rhsMlu is the tightest

	      xtraC [varIndices [i] -> Index ()] += coeff0 * lj;
	      xtraC [varIndices [j] -> Index ()] += coeff0 * ui;
	      rhs += coeff0 * uj * li;

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

  for (std::vector <expression *>::iterator 
	 i  = varIndices . begin (); 
       i   != varIndices . end   (); ++i) {

    int varInd = (*i) -> Index ();

    if (inverse [varInd] >= 0) {
#ifdef DEBUG
      printf ("now adding %g to coeff [inv [%d] = %d]\n", xtraC [varInd], varInd, inverse [varInd]);
#endif
      coeff [inverse [varInd]] += xtraC [varInd];
    } else if (fabs (xtraC [varInd]) > 1e-8) { // independent variable not appearing so far

      coeff   [nterms]   = xtraC [varInd];
      inverse [varInd]   = nterms; // probably useless
      ind     [nterms++] = varInd;
    }
  }

  delete [] inverse;
  delete [] xtraC;

  OsiRowCut *cut = new OsiRowCut;
  cut -> setRow (nterms, ind, coeff);
  cut -> setLb (rhs);

  if (nterms > 0) {

#ifdef DEBUG
    printf ("SDP: separating ");
    cut -> print ();
#endif

    CoinAbsFltEq treatAsSame (1.0e-8);
    int initial = cs.sizeRowCuts();
    cs.insertIfNotDuplicate (*cut, treatAsSame);
    int final = cs.sizeRowCuts();

    if (initial == final) {

      ++ *duplicate_cuts;

      // if flag was false, we still add the duplicate cut
      if (!checkduplicates)
	cs.insert (cut);
    }
  }

  delete cut;
  delete [] ind;
  delete [] coeff;
}


/***********************************************************************/
void CouenneSdpCuts::myremoveBestOneRowCol (double *matrix, 
					    int n, int running_n, 
					    int min_nz,
					    bool *del_idx, 
					    double **sparse_v_mat, 
					    int *card_v_mat, 
					    int *evdec_num) const {

  if (running_n==1) 
    return;

  int 
    best_idx = -1,
    rnsq     = (running_n - 1) * (running_n - 1),
    card_ev_best;

  double

    best_val    = 1,

    *matrixCopy = CoinCopyOfArray (matrix, running_n * running_n),

    *T          = new double [rnsq],
    *Tcopy      = new double [rnsq],
    *Tbest      = new double [rnsq],

    *wbest      = new double [running_n - 1],
    *zbest      = new double [rnsq],

    *w          = new double [running_n - 1],
    *z          = new double [rnsq];

  //while (running_n > 1)

  for (int k=0; k<running_n; k++) {

    /// construct two matrices T and Tcopy that cut k-th row and k-th column

    for (int i=0, ii=0; i<running_n; i++) {

      if (i==k) continue;

      for(int j=0, jj=0; j<running_n;j++) {

	if (j==k) continue;

	int    idx1 = (running_n-1)*ii+jj;
	double val2 = matrixCopy [running_n*i+j];

	T     [idx1] = 
        Tcopy [idx1] = val2;

	jj++;
      }

      ii++;
    }

    int card_ev;

    (*evdec_num)++;

    dsyevx_interface (running_n - 1, T, card_ev, w, z, EV_TOL, -COIN_DBL_MAX, 0., 1, (running_n-1 == min_nz) ? (running_n - 1) : 1);

    double val = w [0]; // minimum eigenvalue

    if (val < 0 && val < best_val) {

      best_val=val;
      best_idx=k;

      std::memcpy (Tbest, Tcopy, rnsq            * sizeof (double));
      std::memcpy (zbest, z,     rnsq            * sizeof (double));
      std::memcpy (wbest, w,     (running_n - 1) * sizeof (double));

      card_ev_best = card_ev;
    }
  }

  if (best_idx >= 0) {
		
    if (del_idx == NULL) {
      del_idx = new bool[n];
      CoinFillN (del_idx, n, false);

      // for (int i=0;i<n;i++)
      // 	del_idx[i]=false;
    }

    int cnt_idx_orig = 0;
    int cnt_idx_minor = 0;

    while (cnt_idx_minor < running_n) {
      if (del_idx[cnt_idx_orig] == false) {
	if (cnt_idx_minor == best_idx) {
	  del_idx[cnt_idx_orig] = true;
	  break;
	}
	else 
	  cnt_idx_minor++;
      }
      cnt_idx_orig++;
    }

    if (running_n - 1 == min_nz) {
      for(int i=0;i<card_ev_best;i++) {
	if (wbest[i] < 0) {
	  double *curr_ev = zbest + (i*(running_n-1));
					
	  for(int j=0;j<n;j++)
	    sparse_v_mat[i][j]=0.0;
					
	  int idx_orig = 0;
	  int idx_minor = 0;

	  while (idx_orig < n) {
	    if (!(del_idx[idx_orig])) {
	      sparse_v_mat[i][idx_orig] = curr_ev[idx_minor];
	      idx_minor++;
	    }
	    idx_orig++;
	  }

	  (*card_v_mat)++;
	}
	else //no more negative eigenvalues
	  break;
      }
    }
    else {
      myremoveBestOneRowCol (Tbest, n, running_n - 1, min_nz, del_idx, sparse_v_mat, card_v_mat, evdec_num);
    }
  }

  // copy Tbest in matrixCopy

  // } // end while

  delete [] del_idx;

  delete [] z;
  delete [] w;

  delete [] T;
  delete [] Tcopy;
  delete [] matrixCopy;

  delete [] Tbest;
  delete [] zbest;
  delete [] wbest;

}// myremoveBestOneRowCol()


/************************************************************************/
void CouenneSdpCuts::sparsify2 (const int n,
				const double *sol, double **sparse_v_mat,
				int *card_v_mat, int min_nz, int *evdec_num) const {
	
  int np = n;
  double *matrix = new double[np*np];

  matrix[0] = 1;
  for (int i=0;i<n;i++)
    matrix[np*(i+1)] = sol[i];
  for (int i=0;i<n;i++) {
    for (int j=i;j<n;j++)
      matrix[(j+1)*np+(i+1)] = sol[indexQ (i, j, n)];
  }

  myremoveBestOneRowCol(matrix, np, np,min_nz,NULL,sparse_v_mat,card_v_mat,evdec_num);

  delete [] matrix;
}// sparsify2()


/************************************************************************/
void CouenneSdpCuts::additionalSDPcuts(const OsiSolverInterface &si,OsiCuts &cs, CouenneSparseMatrix *minor, int np, const double *A, const double *vector, int **indA, int *duplicate_cuts) const{

  int *indices;
  indices = new int[np];
  int cnt = 0;
  for(int i=0;i<np;i++) {
    if (vector[i] != 0.0)
      indices[i] = cnt++;
    else
      indices[i] = -1;
  }
	
  double *subA = new double[cnt*cnt];

  for (register int i=0; i<np; i++) {
    if (indices[i] >= 0) {
      for (register int j=0; j<np; j++) {
	if (indices[j] >= 0)
	  subA[cnt*indices[j]+indices[i]] = A[np*j+i];
      }
    }
  }

  double *w = NULL, *z = NULL;
  int m;

  //////////////////////////////////////////////////////
  dsyevx_interface (cnt, subA, m, w, z, EV_TOL, -COIN_DBL_MAX, 0., 1, cnt);
  //////////////////////////////////////////////////////

  double *v = new double[np];
  double *newv = new double[np];


  for (int k=0; k<m; k++) {
	
#ifdef ONLY_NEG_EIGENV
    if(w [k] > 0) {
      break;
    }
#endif

    double *zbase = z + k * cnt;
    for (int j=0; j<cnt; j++) {
      v [j] = *zbase++;
    }

    for(int j=0;j<np;j++) {
      if (indices[j] >= 0)
	newv[j] = v[indices[j]];
      else
	newv[j] = 0;
    }

    genSDPcut (si, cs, minor, newv, newv, indA, removeduplicates_,duplicate_cuts);
  }

  delete [] v;
  delete [] newv;

  delete [] w;
  delete [] z;

  delete [] subA;
  delete [] indices;
} // additionalSDPcuts


/************************************************************************/
void CouenneSdpCuts::update_sparsify_structures (const int np, const double *sol, double *v,
					 double* margin, double** mat, double *lhs, 
					 const int *zeroed, int evidx, bool decompose, 
					 int *evdec_num) const {

  // copy sol[] in mat[][]
  mat[0][0] = 1;
  for(int i=1; i<np; i++) {
    mat[0][i] = sol[i-1];
    mat[i][0] = sol[i-1];
  }
  for(int i=1; i<np; i++) {
    for(int j=i; j<np; j++) {
      int ind = indexQ(i-1, j-1, np-1);
      mat[i][j] = sol[ind];
      mat[j][i] = sol[ind];
    }
  }

  int minor_n = np;
  if (zeroed != NULL) {
    for(int i=0;i<np;i++)
      if (zeroed[i] == 0)
	minor_n--;
  }

  if ((decompose)  && (minor_n > 2)) {

    /*
      if (minor_n < np) {
      add_v_cut(np, loc_selected, loc_card_selected, locv, 
      init_card_selected, &has_init_vect,
      selected, &card_selected, &new_selected, 
      trace_bin, trace_bin_size,
      sparse_v_mat, card_v_mat);
      }
    */
    decomposition_counter++;
    (*evdec_num)++;
    double *minor_A = new double[np*np];
    double *minor_w = new double[np];
    double *minor_z = new double[np*np];

    //prepare active submatrix (induced by zeroed vector)

    int ii = 0;
    int jj = 0;
    for (int i=0;i<np;i++) {
      if (zeroed[i] == 0)
	continue;
      jj = 0;
      for (int j=0;j<np;j++) {
	if (zeroed[j] == 0)
	  continue;
	minor_A[(minor_n*ii) + jj] = mat[i][j];
	jj++;
      }
      ii++;
    }
    //eigendecomposition
    int m;
    //		dsyevx_wrapper_first_p (minor_n, minor_A, m, minor_w, minor_z,evidx+1,tracer_);
    dsyevx_interface (minor_n, minor_A, m, minor_w, minor_z, EV_TOL, -COIN_DBL_MAX, 0., 1, 1);

    //update v (reindex the evidx-th eigenvector entries)
    ii = 0;
    for (int i=0;i<np;i++) {
      v[i] = 0;
      if (zeroed[i] == 0)
	continue;
      v[i] = minor_z[ii];
      ii++;
    }
    delete [] minor_A;
    delete [] minor_w;
    delete [] minor_z;
  }

  for(int i=0; i<np; i++) {
    for(int j=0; j<np; j++) {
      mat[i][j] *= v[i] * v[j];
      if ((zeroed != NULL) && (zeroed[j] == 0)) {
	mat[i][j] = 0;
	mat[j][i] = 0;
      }
    }
  }

  (*lhs) = 0;
  for(int i=0; i<np; i++) {
    margin[i] = 0;
    for(int j=0; j<np; j++) {
      margin[i] += mat[i][j];	
    }
    (*lhs) += margin[i];
  }
}

/************************************************************************/
void CouenneSdpCuts::zero_comp(const int ind_i, const double delta,
			       const int np, const int *selected,
			       int *loc_selected, 
			       int *ploc_card_selected, int *ploc_card_new_selected, 
			       double *ploc_lhs, 
			       double *locmargin, double **locmat, 
			       const double *sol, double *locv, 
			       const int evidx, bool wise, int *evdec_num, double *recomp_gap, double *threshold) const {

  //double  curr_lhs = (*ploc_lhs);
  static int zerocount;
  bool local_wise = false;
  if ((wise) && ((*ploc_lhs)-delta > (*threshold))) {
    (*threshold) = (*ploc_lhs)-delta + (*recomp_gap);
    local_wise = true;
  }

  zerocount++;

  loc_selected[ind_i] = 0;
  (*ploc_card_selected)--;
  
  if(selected[ind_i] != 1) {
    (*ploc_card_new_selected)--;
  }
  (*ploc_lhs) -= delta;

  update_sparsify_structures (np,sol,locv,locmargin,locmat,ploc_lhs, loc_selected, evidx, local_wise, evdec_num);

} /* zero_comp */

/************************************************************************/
void CouenneSdpCuts::zero_valid_delta(const int np, const int *order,
				      const int * selected,
				      const int min_card_new_selected,
				      const double min_delta, const int start_point, 
				      const int curr_i, 
				      int *loc_selected, 
				      int *ploc_card_selected, 
				      int *ploc_card_new_selected, 
				      double *ploc_lhs, 
				      double *locmargin, double **locmat, 
				      int *pnchanged, 
				      const double *sol, double *locv, 
				      const int evidx, bool wise,double *recomp_gap, double *threshold,
				      int *pcard_selected,
				      int *pnew_selected,
				      double **sparse_v_mat,
				      int *pcard_v_mat,
				      const int init_card_selected, int *has_init_vect,
				      int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged = 0);
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];
    int skip = 0;

    if((selected[ind_i] == 0) && 
       (min_card_new_selected >= *ploc_card_new_selected)) {
      skip = 1;
    }

    if((skip) || (curr_ind == start_point) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(*ploc_lhs - delta < min_delta) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx, wise, evdec_num , recomp_gap,threshold);
      (*pnchanged)++;

    }
  }
} /* zero_valid_delta */

/************************************************************************/
void CouenneSdpCuts::zero_selected(const int np, const int *order,
				   const int *selected,
				   const int min_card_new_selected,
				   const double min_delta, const int start_point,
				   const int curr_i, 
				   int *loc_selected, int *ploc_card_selected, 
				   int *ploc_card_new_selected, 
				   double *ploc_lhs, 
				   double *locmargin, double **locmat, 
				   int *pnchanged, 
				   const double *sol, double *locv, 
				   const int evidx, bool wise,double *recomp_gap, double *threshold,
				   int *pcard_selected,
				   int *pnew_selected,
				   double **sparse_v_mat,
				   int *pcard_v_mat,
				   const int init_card_selected, int *has_init_vect,
				   int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged = 0);
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];

    if((selected[ind_i] == 0) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(*ploc_lhs - delta < min_delta) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx,wise,evdec_num,recomp_gap,threshold);
      (*pnchanged)++;
    } 
  }
} /* zero_selected */


/************************************************************************/
void CouenneSdpCuts::zero_pos_delta(const int np, const int *order,
				    const int *selected,
				    const int min_card_new_selected,
				    const int start_point, const int curr_i, 
				    int *loc_selected, int *ploc_card_selected, 
				    int *ploc_card_new_selected, 
				    double *ploc_lhs, 
				    double *locmargin, double **locmat, 
				    int *pnchanged, 
				    const double *sol, double *locv, 
				    const int evidx, bool wise, double *recomp_gap, double *threshold,
				    int *pcard_selected,
				    int *pnew_selected,
				    double **sparse_v_mat,
				    int *pcard_v_mat,
				    const int init_card_selected, int *has_init_vect,
				    int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged) = 0;
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];
    
    int skip = 0;

    if((selected[ind_i] == 0) && 
       (min_card_new_selected >= *ploc_card_new_selected)) {
      skip = 1;
    }

    if((skip) || (curr_ind == start_point) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(delta > 0) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx,wise,evdec_num,recomp_gap,threshold);
      (*pnchanged)++;

    } 
  }
} /* zero_pos_delta */


/************************************************************************/
void CouenneSdpCuts::add_v_cut(const int np,
			       const int *loc_selected, 
			       const int loc_card_selected,
			       const double *locv,
			       const int init_card_selected, int *has_init_vect,
			       int *selected, int *pcard_selected,
			       int *pnew_selected,
			       double **sparse_v_mat,
			       int *pcard_v_mat) const {

  (*pnew_selected) = 0;

  for(int i=0; i<np; i++) {
    if(loc_selected[i]) {
      sparse_v_mat[*pcard_v_mat][i] = locv[i];
      if(selected[i] == 0) {
	selected[i] = 1;
	(*pcard_selected)++;
	(*pnew_selected)++;
      }
    }
    else {
      sparse_v_mat[*pcard_v_mat][i] = 0;
    }
  }

#ifdef NORMALIZE_SPARSE_CUTS
  //normalization (setting vector norm to 1)
  double curr_norm = 0.0;
  for (int i=0;i<np;i++) {
    curr_norm += fabs(sparse_v_mat[*pcard_v_mat][i]);
  }
  for (int i=0;i<np;i++) {
    if (sparse_v_mat[*pcard_v_mat][i] != 0.0)
      sparse_v_mat[*pcard_v_mat][i] = sparse_v_mat[*pcard_v_mat][i]/curr_norm;
  }
#endif

  if(loc_card_selected + init_card_selected == np) {
    if(*has_init_vect == 1) {

      return;
    }
    else {
      (*has_init_vect) = 1;
    }
  }
	
  (*pcard_v_mat)++;
  
} /* add_v_cut */


/************************************************************************/
void CouenneSdpCuts::sparsify (bool use_new_sparsify,
			       const int evidx, const double eigen_val, 
			       const double *v, const int n,
			       const double *sol, double **sparse_v_mat,
			       int *card_v_mat, int *evdec_num) const {

  int i, j, np = n, nchanged = 0;
  double sq_np = sqrt((double)np);

  double min_delta;
  double is_zero = 1/(10 * sq_np);
  int min_number_new_per_cut = 1;
	
  int *selected = new int[np], card_selected = 0;
  int *loc_selected = new int[np], loc_card_selected = 0;
  int loc_card_new_selected = 0;
	
  double lhs = 0, loc_lhs = 0;
  double *margin = new double[np];
  double *locv = new double[np];
  double *locv_orig = new double[np];
  double *locmargin = new double[np];
  double **mat = new double*[np];
  double **locmat = new double*[np];
	
  //int seed = 225535;
  int *order = new int[np];
  double *rand_val = new double[np];

  *card_v_mat = 0;
	
  for (i=0; i<np; i++) {
    selected[i] = 0;
    mat[i] = new double[np];
    locmat[i] = new double[np];
    order[i] = i;
    rand_val[i] = drand48 (); //cpp_genalea(&seed);
		
    // zero small components in v
    if(fabs(v[i]) < is_zero) {

      locv_orig[i] = 0;
      selected[i] = -1; // -1: ind will be set to 0 in loc_selected
      card_selected++;
    } else {
      locv_orig[i] = v[i];
    }
  }

  /// Knuth shuffle for creating a random order
  for (int i=0; i<np; ++i) {

    int 
      newpos = i + (int) floor (((double) (np - i) - 1.e-3) * drand48 ()),
      tmp    = order [newpos];

    order [newpos] = order [i];
    order [i] = tmp;
  }

  update_sparsify_structures (np,sol,locv_orig,margin,mat,&lhs, NULL, evidx, false, evdec_num);

  int init_card_selected = card_selected; // to recognize if cut from original
  // vector is generated
  int has_init_vect = 0;
	
  min_delta = lhs * (use_new_sparsify ? SPARSIFY_NEW_DELTA : SPARSIFY_OLD_DELTA); // do not weaken the cut too much
  int start_point = -1; // order[start_point]: index that should not be removed
	
  while(card_selected < np) {

    for(i=0; i<np; i++) {
      if(selected[order[i]] == 0) {
	start_point = i;
	break;
      }
    }
    
    loc_card_selected = np;
    loc_card_new_selected = np;
    loc_lhs = lhs;
    double recomp_gap = fabs(lhs*WISE_SPARSIFY_GAP);
    double threshold = lhs + recomp_gap;

    // restore locv (might have been changed by WISE_SPARSIFY during the generation of the last sparse cut)
    for(i=0;i<np;i++)
      locv[i] = locv_orig[i];

    for(i=0; i<np; i++) {
      if(selected[i] == -1) {
	loc_selected[i] = 0;
	loc_card_selected--;
	loc_card_new_selected--;
      } else {
	loc_selected[i] = 1;
				
	if(selected[i] == 1) {
	  loc_card_new_selected--;
	}
      }
      locmargin[i] = margin[i];
      for(j=0; j<np; j++) {
	locmat[i][j] = mat[i][j];
      }
    }

    if(loc_lhs < min_delta) {

      int changed = 1;
			
      while(changed) {

	int curr_i = start_point;
				
	changed = 0;
				
	int sel_nchanged = -1; 

	while(sel_nchanged != 0) {
	  int new_selected = 0;
	  zero_selected(np, order, selected, min_number_new_per_cut,
			min_delta, start_point,
			curr_i, loc_selected, 
			&loc_card_selected, &loc_card_new_selected, 
			&loc_lhs, locmargin, locmat, 
			&sel_nchanged,sol,locv,evidx, use_new_sparsify, &recomp_gap,&threshold,
			&card_selected, &new_selected, 
			sparse_v_mat, card_v_mat,
			init_card_selected, &has_init_vect,evdec_num);
					
	  if(sel_nchanged) {
	    nchanged += sel_nchanged;
	    //	    changed = 1;
	  }
	} // while(sel_nchanged != 0) 

	int pos_nchanged = -1;
	
	while(pos_nchanged != 0) {
	  int new_selected = 0;
	  zero_pos_delta(np, order, selected, min_number_new_per_cut,
			 start_point, start_point, loc_selected, 
			 &loc_card_selected, &loc_card_new_selected, 
			 &loc_lhs, locmargin, locmat, 
			 &pos_nchanged,sol,locv,evidx,use_new_sparsify, &recomp_gap,&threshold,
			 &card_selected, &new_selected, 
			 sparse_v_mat, card_v_mat,
			 init_card_selected, &has_init_vect,evdec_num);
				
	  if(pos_nchanged) {
	    nchanged += pos_nchanged;
	    changed = 1;
	  }
	} /* while(pos_nchanged != 0) */

	if(changed) {
	  continue;
	}
				

	curr_i = start_point;
				
	int val_nchanged = -1;

	if(val_nchanged) {
	  int new_selected = 0;
	  zero_valid_delta(np, order, selected, min_number_new_per_cut,
			   min_delta, start_point,
			   curr_i, loc_selected, 
			   &loc_card_selected, &loc_card_new_selected, 
			   &loc_lhs, locmargin, locmat, 
			   &val_nchanged,sol,locv,evidx, use_new_sparsify, &recomp_gap,&threshold,
			   &card_selected, &new_selected, 
			   sparse_v_mat, card_v_mat,
			   init_card_selected, &has_init_vect,evdec_num);
					
	  if(val_nchanged) {
	    nchanged += val_nchanged;
	    changed = 1;
	  }
	}


      } /* while(changed) */

      if((loc_card_selected < np * (use_new_sparsify ? SPARSIFY_NEW_NZ_THRESHOLD : SPARSIFY_OLD_NZ_THRESHOLD)) || (*card_v_mat == 0)) {
				
	int new_selected = 0;
				
	add_v_cut(np, loc_selected, loc_card_selected, locv, 
		  init_card_selected, &has_init_vect,
		  selected, &card_selected, &new_selected, 
		  sparse_v_mat, card_v_mat);
      } else {
	selected[order[start_point]] = 1;
	card_selected++;
      }
    } else {
      // loc_lhs >= min_delta  use vector as is

      card_selected = np;

#ifndef ONLY_NEG_EIGENV
      int new_selected = 0;
			
      add_v_cut(np, loc_selected, loc_card_selected, locv, 
		init_card_selected, &has_init_vect,
		selected, &card_selected, &new_selected, 
		sparse_v_mat, card_v_mat);
#endif

    }
  } /* while(card_selected < np) */

  delete[] order;
  delete[] rand_val;
	
  for (i=0; i<np; i++) {
    delete [] mat[i];
    delete [] locmat[i];
  }
  delete [] mat;
  delete [] locmat;

  delete[] locv;
  delete[] locv_orig;
  delete[] margin;
  delete[] locmargin;
	
  delete[] selected;
  delete[] loc_selected;
}
