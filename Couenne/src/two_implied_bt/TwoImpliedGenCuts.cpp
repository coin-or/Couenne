/* $Id$
 *
 * Name:    TwoImpliedGenCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate cuts using two inequalities from the LP relaxation
 * 
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#include <stdlib.h>

#include "BonCbc.hpp"
#include "BonBabInfos.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglCutGenerator.hpp"

#include "CouenneProblemElem.hpp"
#include "CouenneTwoImplied.hpp"
#include "CouenneExprVar.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

// do single pair of inequalities

int combine (OsiCuts &cs, 
	     int n1, int n2, 
	     const int *ind1, // indices
	     const int *ind2, 
	     double *sa1, // coeff (sparse array)
	     double *sa2,
	     const double *a1,  // coeff
	     const double *a2, 
	     double *clb, // variable bounds
	     double *cub,
	     double l1, // constraint bounds
	     double l2,  
	     double u1, 
	     double u2, 
	     int sign); // invert second constraint? -1: yes, +1: no

/// the main CglCutGenerator
void CouenneTwoImplied::generateCuts (const OsiSolverInterface &si, 
				      OsiCuts &cs, 
				      const CglTreeInfo info) const {

  static int nBadColMatWarnings = 0;

  std::set <std::pair <int, int> > pairs;

  const CoinPackedMatrix *colA = si. getMatrixByCol ();

  int 
    n = colA -> getMajorDim (), // # cols
    m = colA -> getMinorDim (); // # rows

  double 
    *sa1 = new double [n], // contains dense representation of a1 i.e. lots of zeros
    *sa2 = new double [n]; //                                  a2

  CoinFillN (sa1, n, 0.);
  CoinFillN (sa2, n, 0.);

  const double
    *A   = colA -> getElements (),
    *rlb = si.getRowLower (),
    *rub = si.getRowUpper ();

  const int *ind = colA -> getIndices      (),
            *sta = colA -> getVectorStarts ();

  // For every column i, compare pairs of rows j and k with nonzero
  // coefficients.
  //
  // If the coefficients have the same sign and the inequalities are
  // both >= or both <=, skip this pair -- no tightening can come from
  // this pair (this doesn't mean that for another variable the
  // opposite may happen).

  for (int i=0; i<n; i++, sta++) {

    int nEl = *(sta+1) - *sta;
    //printf ("column %d: %d elements %d -> %d\n", i, nEl, *sta, *(sta+1));

    for   (int jj = nEl,  j = *sta; jj--; j++)
      for (int kk = jj,   k = j+1;  kk--; k++) {

	register int 
	  indj = ind [j],
	  indk = ind [k];

	// should never happen, but if it does, just bail out
 	if ((indj > m) || (indj < 0) ||
 	    (indk > m) || (indk < 0)) {

	  if (nBadColMatWarnings++ <= 1)
	    printf ("Couenne: warning, matrix by row has nonsense indices. Skipping\n");

	  return; 
	}

	double prod = A [j] * A [k];

	if (prod > 0.) { // same sign -- skip unless finite lb1/ub2 OR
			 // finite ub1/lb2. This is to avoid a situation
			 // in which all coefficients in this pair
			 // have the same sign

	  if (!(
		((rlb [indj] > -COUENNE_INFINITY) && (rub [indk] <  COUENNE_INFINITY)) || 
		((rub [indj] <  COUENNE_INFINITY) && (rlb [indk] > -COUENNE_INFINITY))
		)
	      )
	    continue;

	} else

	if ((prod < 0.) && // opposite sign -- multiply second
			   // inequality by -1 and repeat
	    !(
	      ((rlb [indj] > -COUENNE_INFINITY) && (rlb [indk] > -COUENNE_INFINITY)) || 
	      ((rub [indj] <  COUENNE_INFINITY) && (rub [indk] <  COUENNE_INFINITY))
	     )
	    )
	  continue;

	pairs.insert (std::pair <int, int> (indj, indk));
      }
  }

  // pairs (h,k) are done. Now for each pair set new bounds, if possible ///////////////////////////

  const CoinPackedMatrix *rowA = si. getMatrixByRow ();

  const double
    *rA  = rowA -> getElements ();

  double
    *clb = CoinCopyOfArray (si.getColLower (), n),
    *cub = CoinCopyOfArray (si.getColUpper (), n);

  const int *rInd = rowA -> getIndices      (),
            *rSta = rowA -> getVectorStarts ();

  int 
    ntightened = 0,
    ntrials    = 0,
    nCurTightened;

  // info about LP problem: upper bound, dual bound

  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (si.getAuxiliaryInfo ());

  // data structure for FBBT

  t_chg_bounds *chg_bds = new t_chg_bounds [n];

  for (int i=0; i < n; i++) 
    if (problem_ -> Var (i) -> Multiplicity () <= 0) {
      chg_bds [i].setLower (t_chg_bounds::UNCHANGED);
      chg_bds [i].setUpper (t_chg_bounds::UNCHANGED);
    }

  do {

    nCurTightened = 0;

    // scan all pairs. All are potential pairs of inequalities that
    // can give a better (combined) implied bound

    for (std::set <std::pair <int, int> >:: iterator p = pairs.begin (); p != pairs.end (); ++p) {

      // indices of the two inequalities

      int 
	h = p -> first,
	k = p -> second;

      const double
	l1 = rlb [h], u1 = rub [h],
	l2 = rlb [k], u2 = rub [k],
	*a1 = rA + rSta [h],
	*a2 = rA + rSta [k];

      const int
	n1 = rSta [h+1] - rSta [h],
	n2 = rSta [k+1] - rSta [k],
	*ind1 = rInd + rSta [h],
	*ind2 = rInd + rSta [k];

      // fill in sa1 and sa2

      for (int i=n1; i--;) sa1 [ind1 [i]] = a1 [i];
      for (int i=n2; i--;) sa2 [ind2 [i]] = a2 [i];

      // A few cases here on the inequalities' bounds, which may
      // determine a change of sign in the second when combining the
      // two: for example, if we have -inf < ax < b and c < dx < +inf,
      // it makes sense to combine -inf < ax < b and -inf < - dx < -c
      // as the upper bound will create a finite convex combination.

      if ((u1 <   COUENNE_INFINITY && u2 <   COUENNE_INFINITY) ||
	  (l1 > - COUENNE_INFINITY && l2 > - COUENNE_INFINITY))

	nCurTightened += combine (cs, n1, n2, ind1, ind2, sa1, sa2, a1, a2, clb, cub, l1, l2, u1, u2, 1);

      // fill in sa2 with opposite values
      for (int i=n2; i--;) sa2 [ind2 [i]] = - a2 [i];

      if ((u1 <   COUENNE_INFINITY && l2 > - COUENNE_INFINITY) ||
	  (l1 > - COUENNE_INFINITY && u2 <   COUENNE_INFINITY))

	// do NOT invert l2 and u2, this is done in combine
	nCurTightened += combine (cs, n1, n2, ind1, ind2, sa1, sa2, a1, a2, clb, cub, l1, l2, u1, u2, -1);

      // clean sa1 and sa2

      for (int i=n1; i--;) sa1 [ind1 [i]] = 0.;
      for (int i=n2; i--;) sa2 [ind2 [i]] = 0.;
    }

    int objInd = problem_ -> Obj (0) -> Body () -> Index ();

    if (nCurTightened &&
	(objInd >= 0) && 
	babInfo && 
	babInfo -> babPtr ()) {

      printf ("FBBT\n");

      CouNumber
	UB      = babInfo -> babPtr () -> model (). getObjValue(),
	LB      = babInfo -> babPtr () -> model (). getBestPossibleObjValue (),
	primal0 = problem_ -> Ub (objInd), 
	dual0   = problem_ -> Lb (objInd);

      // Do one round of BT

      if ((UB < COUENNE_INFINITY) && 
	  (UB < primal0 - COUENNE_EPS)) { // update primal bound (MIP)

	problem_ -> Ub (objInd) = UB;
	chg_bds [objInd].setUpper (t_chg_bounds::CHANGED);
      }

      if ((LB > - COUENNE_INFINITY) && 
	  (LB > dual0 + COUENNE_EPS)) { // update dual bound
	problem_ -> Lb (objInd) = LB;
	chg_bds [objInd].setLower (t_chg_bounds::CHANGED);
      }
    
      if (!(problem_ -> btCore (chg_bds))) {

	//problem infeasible, add IIS of size 2

	OsiColCut *infeascut = new OsiColCut;

	if (infeascut) {
	  int i=0;
	  double upper = -1., lower = +1.;
	  infeascut -> setLbs (1, &i, &lower);
	  infeascut -> setUbs (1, &i, &upper);
	  cs.insert (infeascut);
	  delete infeascut;
	}
      }
    }

    ntightened += nCurTightened;

  } while (++ntrials < nMaxTrials_ && nCurTightened);

  // check if bounds improved, in that case create OsiColCuts

  if (ntightened) {

    const double 
      *oldLB = si. getColLower (),
      *oldUB = si. getColUpper ();

    // check old and new bounds

    int 
      *indLB = new int [n],
      *indUB = new int [n],
      ntightenedL = 0,
      ntightenedU = 0;

    double 
      *valLB = new double [n],
      *valUB = new double [n];

    for (int i=0; i<n; i++) {

      if (clb [i] > oldLB [i]) {
	indLB [ntightenedL]   = i;
	valLB [ntightenedL++] = clb [i];
      }

      if (cub [i] < oldUB [i]) {
	indUB [ntightenedU]   = i;
	valUB [ntightenedU++] = cub [i];
      }
    }

    if (ntightenedL || ntightenedU) {

      OsiColCut newBound;

      newBound.setLbs (ntightenedL, indLB, valLB);
      newBound.setUbs (ntightenedU, indUB, valUB);

      cs.insert (newBound);
    }

    delete [] indLB;
    delete [] indUB;
    delete [] valLB;
    delete [] valUB;
  }

  delete [] clb;
  delete [] cub;
  delete [] sa1;
  delete [] sa2; 
}
