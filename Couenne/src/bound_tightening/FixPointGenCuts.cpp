/* $Id$
 *
 * Name:    FixPointGenCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: Fix point bound tightener -- separator
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneFixPoint.hpp"
#include "CouenneExprVar.hpp"

using namespace Couenne;

/// Cut Generator for FBBT fixpoint

void createRow (int,
		int,
		int,
		OsiSolverInterface *,
		const int    *,
		const double *,
		const double,
		const int,
		bool);

void CouenneFixPoint::generateCuts (const OsiSolverInterface &si,
				    OsiCuts &cs,
				    const CglTreeInfo treeInfo) const {

  /// An LP relaxation of a MINLP problem is available in the first
  /// parameter passed. Let us suppose that this LP relaxation is of the form
  ///
  /// LP = {x in R^n: Ax <= b}
  /// 
  /// For suitable nxm matrix A, rhs vector b, and variable vector
  /// x. Our purpose is that of creating a much larger LP that will
  /// help us find the interval [l,u] corresponding to the fixpoint of
  /// an FBBT algorithm. To this purpose, consider a single constraint
  /// of the above system:
  ///
  /// sum {i=1..n} a_ji x_i <= b_j
  ///
  /// According to two schools of thought (Leo's and mine), this
  /// single constraint can give rise to a number of FBBT
  /// constraints. The two schools of thoughts differ in the meaning
  /// of b: in mine, it is constant. In Leo's, it must be defined as a
  /// variable.

  OsiSolverInterface *fplp = NULL;

  if (true) { // placeholder for later selection of LP solver among
	      // those available

    fplp = new OsiClpSolverInterface;
  }

  /// We need to perform the following steps:
  ///
  /// define variables xL and xU
  /// define variables gL and gU for constraints (downward variables)
  ///
  /// add objective function sum_j (u_j - l_j)
  ///
  /// for each constraint a^j x <= b_j in Ax <= b:
  ///   for each variable x_i contained:
  ///     depending on sign of a_ji, add constraint on x_i^L or x_i^U
  ///     (*) add constraints on g_j as well
  ///
  /// solve LP
  ///
  /// if new bounds are better than si's old bounds
  ///   add OsiColCuts

  /// Get the original problem's coefficient matrix and rhs vector, A and b

  const CoinPackedMatrix *A = si. getMatrixByRow ();

  const int
       n = si.  getNumCols  (),
       m = si.  getNumRows  (),
    *ind = A -> getIndices  ();

  const double
    *lb  = si.  getColLower (),
    *ub  = si.  getColUpper (),
    *rlb = si.  getRowLower (),
    *rub = si.  getRowUpper (),
    *coe = A -> getElements ();

  //for (int i=0; i<n; i++) 
  //printf ("----------- x_%d in [%g,%g]\n", i, lb [i], ub [i]);

  // add lvars and uvars to the new problem
  for (int i=0; i<n; i++) {
    bool isActive = problem_ -> Var (i) -> Multiplicity () > 0;
    fplp -> addCol (0, NULL, NULL, isActive ? lb [i] : 0., isActive ? ub [i] : 0., -1.); // xL_i
  }

  for (int i=0; i<n; i++) {
    bool isActive = problem_ -> Var (i) -> Multiplicity () > 0;
    fplp -> addCol (0, NULL, NULL, isActive ? lb [i] : 0., isActive ? ub [i] : 0., +1.); // xU_i
  }

  if (extendedModel_) {

    for (int j=0; j<m; j++) fplp -> addCol (0, NULL, NULL, rlb [j],       COIN_DBL_MAX, 0.); // bL_j
    for (int j=0; j<m; j++) fplp -> addCol (0, NULL, NULL, -COIN_DBL_MAX, rub [j],      0.); // bU_j
  }

  // Scan each row of the matrix 
  
  for (int j=0; j<m; j++) { // for each row

    int nEl = A -> getVectorSize (j); // # elements in each row

    //printf ("row %4d: %4d elements\n", j, nEl);

    // for (int i=0; i<nEl; i++) {
    //   printf ("%+g x%d ", coe [i], ind [i]);
    //   fflush (stdout);
    // }

    // printf ("\n");

    if (!nEl)
      continue;

    // create cuts for the xL and xU elements //////////////////////

    if (extendedModel_ || rlb [j] > -COUENNE_INFINITY) 
      for (int i=0; i<nEl; i++) 
	createRow (-1, ind [i], n, fplp, ind, coe, rlb [j], nEl, extendedModel_); // upward constraints -- on xL_i

    if (extendedModel_ || rub [j] <  COUENNE_INFINITY) 
      for (int i=0; i<nEl; i++) 
	createRow (+1, ind [i], n, fplp, ind, coe, rub [j], nEl, extendedModel_); // upward constraints -- on xU_i

    // create (at most 2) cuts for the bL and bU elements //////////////////////

    if (extendedModel_) {
      createRow (-1, 2*n     + j, n, fplp, ind, coe, rlb [j], nEl, extendedModel_); // downward constraints -- on bL_i
      createRow (+1, 2*n + m + j, n, fplp, ind, coe, rub [j], nEl, extendedModel_); // downward constraints -- on bU_i
    }

    ind += nEl;
    coe += nEl;
  }

  // finally, add consistency cuts, bL <= bU

  if (extendedModel_)

    for (int j=0; j<m; j++) { // for each row

      int    ind [2] = {2*n + j, 2*n + m + j};
      double coe [2] = {1.,      -1.};

      CoinPackedVector row (2, ind, coe);
      fplp -> addRow (row, -COIN_DBL_MAX, 0.);
    }

  /// Now we have an fbbt-fixpoint LP problem. Solve it to get
  /// (possibly) better bounds

  fplp -> setObjSense (-1.); // we want to maximize 

  // printf ("writing lp\n");
  // fplp -> writeLp ("fplp");

  fplp -> initialSolve ();

  const double 
    *newLB = fplp -> getColSolution (),
    *newUB = newLB + n,
    *oldLB = si. getColLower (),
    *oldUB = si. getColUpper ();

  // check old and new bounds

  for (int i=0; i<n; i++) {

    // printf ("x%d: [%g,%g] --> [%g,%g]\n", i, 
    // 	    oldLB [i], oldUB [i], 
    // 	    newLB [i], newUB [i]);

    if (newLB [i] > oldLB [i] + COUENNE_EPS) {

      OsiColCut newBound;
      newBound.setLbs (1, &i, newLB + i);
      cs.insert (newBound);
    }

    if (newUB [i] < oldUB [i] - COUENNE_EPS) {

      OsiColCut newBound;
      newBound.setUbs (1, &i, newUB + i);
      cs.insert (newBound);
    }
  }
}


// single cut creation. Parameters:
//
// 1) sign:     tells us whether this is a <= or a >= (part of a) constraint.
// 2) indexVar: index of variable we want to do upward or downward bt
// 3) nVars:    number of variables in the original problems (original +
//              auxiliaries). Used to understand if we are adding an
//              up or a down constraint
// 4) p:        solver interface to which we are adding constraints
// 5) indices:  vector containing indices of the linearization constraint (the i's)
// 6) coe:                        coeff                                        a_ji's
// 7) rhs:      right-hand side of constraint
// 8) nEl:      number of elements of this linearization cut
// 9) extMod:   extendedModel_

void createRow (int sign,
		int indexVar,
		int nVars,
		OsiSolverInterface *p,
		const int *indices,
		const double *coe,
		const double rhs,
		const int nEl,
		bool extMod) {

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Define as I+ the set of indices in {1..n} such that a_ji > 0 and i != indexVar
  /// and I- viceversa. Consider now
  /// 
  /// sum {i=1..n} a_ji x_i = b_j,      (=)
  ///
  /// equivalent to the two constraints
  ///
  /// sum {i=1..n} a_ji x_i <= b_j      (<)
  /// sum {i=1..n} a_ji x_i >= b_j.     (>)
  ///
  /// Hence we should consider both (<) and (>) when we have an
  /// equality constraint. The resulting set of upward FBBT is as
  /// follows:
  ///
  /// sum {i in I+} a_ji xU_i + sum {i in I-} a_ji xL_i >= bU_j  (only if (<) enforced, i.e., sign ==  1)
  /// sum {i in I+} a_ji xL_i + sum {i in I-} a_ji xU_i <= bL_j  (only if (>) enforced, i.e., sign == -1)
  ///
  /// together with the constraints defining the initial bounding
  /// interval of the auxiliary variable (already included):
  ///
  /// bU_j <= bU0_j (<)
  /// bL_j >= bL0_j (>)
  ///
  /// The set of downward FBBT constraints is instead:
  ///
  /// xL_i >= (bL_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k) / a_ji   (if a_ji > 0 and (>))
  /// xU_i <= (bU_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k) / a_ji   (if a_ji > 0 and (<))
  ///
  /// xU_i <= (bL_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k) / a_ji   (if a_ji < 0 and (>))
  /// xL_i >= (bU_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k) / a_ji   (if a_ji < 0 and (<))
  ///
  /// probably better rewritten, to avoid (further) numerical issues, as
  ///
  ///   a_ji xL_i >=   bL_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k   (if a_ji > 0 and (>))
  ///   a_ji xU_i <=   bU_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k   (if a_ji > 0 and (<))
  ///
  /// - a_ji xU_i <= - bL_j + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k   (if a_ji < 0 and (>))
  /// - a_ji xL_i >= - bU_j + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k   (if a_ji < 0 and (<))
  ///
  /// or even better, to keep the old coefficients (but not the indices!), like this:
  ///
  /// a_ji xL_i + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k - bL_j >= 0  (if a_ji > 0 and (>))
  /// a_ji xU_i + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k - bU_j <= 0  (if a_ji > 0 and (<))
  ///
  /// a_ji xU_i + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k - bL_j >= 0  (if a_ji < 0 and (>))
  /// a_ji xL_i + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k - bU_j <= 0  (if a_ji < 0 and (<))
  ///
  /// and finally we need initial lower/upper bounds:
  ///
  /// xL_i >= xL0_i
  /// xU_i <= xU0_i
  ///
  /// and some consistency constraints
  ///
  /// bL_i <= bU_i
  ///
  /// (these and the two bound constraints above are already added in
  /// the main function above).
  ///
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// According to my school of thought, instead, there is no
  /// upward/downward directions to simulate. Hence, considering again
  ///
  /// sum {i=1..n} a_ji x_i <= b_j      (<) -- sign will be -1
  /// sum {i=1..n} a_ji x_i >= b_j      (>) -- sign will be +1
  ///
  /// we'll have similar constraints, except bL and bU no longer exist
  /// but are replaced by the original rhs.
  ///
  /// xL_i >= (b_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k) / a_ji   (if a_ji > 0 and (>))
  /// xU_i <= (b_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k) / a_ji   (if a_ji > 0 and (<))
  ///										     
  /// xU_i <= (b_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k) / a_ji   (if a_ji < 0 and (>))
  /// xL_i >= (b_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k) / a_ji   (if a_ji < 0 and (<))
  ///
  /// once again rewritten as 
  ///
  ///   a_ji xL_i >=   b_j - sum {k in I+} a_jk xU_k - sum {k in I-} a_jk xL_k   (if a_ji > 0 and (>))
  ///   a_ji xU_i <=   b_j - sum {k in I+} a_jk xL_k - sum {k in I-} a_jk xU_k   (if a_ji > 0 and (<))
  ///										     
  /// - a_ji xU_i <= - b_j + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k   (if a_ji < 0 and (>))
  /// - a_ji xL_i >= - b_j + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k   (if a_ji < 0 and (<))
  ///
  /// or even better:
  ///
  /// a_ji xL_i + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k >= b_j       (if a_ji > 0 and (>))
  /// a_ji xU_i + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k <= b_j       (if a_ji > 0 and (<))
  ///										     
  /// a_ji xU_i + sum {k in I+} a_jk xU_k + sum {k in I-} a_jk xL_k >= b_j       (if a_ji < 0 and (>))
  /// a_ji xL_i + sum {k in I+} a_jk xL_k + sum {k in I-} a_jk xU_k <= b_j       (if a_ji < 0 and (<))
  ///
  /// No other cuts are needed.
  ///
  /////////////////////////////////////////////////////////////////////////////////////////////////////


  // printf ("creating constraint from: ");
  // for (int i=0; i<nEl; i++)
  //   printf ("%+g x%d ", coe [i], indices [i]);
  // printf ("%c= %g\n", sign > 0 ? '<' : '>', rhs);

  int nTerms = nEl + (extMod ? 1 : 0); // always add one element when using extended model

  int    *iInd = new int    [nTerms];
  double *elem = new double [nTerms];

  // coefficients are really easy

  CoinCopyN (coe, nEl, elem);

  if (extMod) {
    elem [nEl] = -1.; // extended model, coefficient for bL or bU
    iInd [nEl] = indexVar;
  }

  // indices are not so easy...

  for (int k=0; k<nEl; k++) {

    int curInd = indices [k];

    iInd [k] = curInd; // Begin with xL_i, same index as x_i in the
                       // original model. Should add n depending on a
                       // few things... 

    if (curInd == indexVar) { // x_k is x_i itself
      if (((sign > 0) && (coe [indexVar] > 0.)) || 
	  ((sign < 0) && (coe [indexVar] < 0.)))

      iInd [k] += nVars;

    } else if (((coe [curInd] > 0.) && (sign < 0)) ||
	       ((coe [curInd] < 0.) && (sign > 0)))

      iInd [k] += nVars;
  }

  CoinPackedVector vec (nTerms, iInd, elem);

  double 
    lb = sign > 0 ? -COIN_DBL_MAX : extMod ? 0. : rhs,
    ub = sign < 0 ? +COIN_DBL_MAX : extMod ? 0. : rhs;

  p -> addRow (vec, lb, ub);

  // OsiRowCut *cut = new OsiRowCut (lb, ub,
  // 				  nTerms, nTerms,
  // 				  iInd, elem);

  // cut -> print ();

  // delete cut;
}
