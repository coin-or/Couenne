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

#include "CglCutGenerator.hpp"
#include "CouenneTwoImplied.hpp"
#include "CoinPackedMatrix.hpp"
#include "CouennePrecisions.hpp"

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
    *sa1 = new double [n], // contains dense representation of a1
    *sa2 = new double [n]; //                                  a2

  CoinFillN (sa1, n, 0.);
  CoinFillN (sa2, n, 0.);

  const double
    *A   = colA -> getElements (),
    *rlb = si.getRowLower (),
    *rub = si.getRowUpper ();

  const int *ind = colA -> getIndices      (),
            *sta = colA -> getVectorStarts ();

  for (int i=0; i<n; i++, sta++) {

    int nEl = *(sta+1) - *sta;
    //    printf ("column %d: %d elements %d -> %d\n", i, nEl, *sta, *(sta+1));

    for   (int jj = nEl,  j = *sta; jj--; j++)
      for (int kk = jj,   k = j+1;  kk--; k++) {

	// should never happen, but if it does, just bail out
	if (ind [j] > m || ind [j] < 0 ||
	    ind [k] > m || ind [k] < 0) {

	  if (nBadColMatWarnings++ <= 1)
	    printf ("Couenne: warning, matrix by row has nonsense indices. Skipping\n");

	  return; 
	}

	double prod = A [j] * A [k];

	if ((prod > 0.) &&
	    (((rlb [ind [j]] < -COUENNE_INFINITY) && (rlb [ind [k]] < -COUENNE_INFINITY)) ||
	     ((rub [ind [j]] >  COUENNE_INFINITY) && (rub [ind [k]] >  COUENNE_INFINITY))))
	  continue;

	if ((prod < 0.) &&
	    (((rlb [ind [j]] < -COUENNE_INFINITY) && (rub [ind [k]] >  COUENNE_INFINITY)) ||
	     ((rlb [ind [k]] < -COUENNE_INFINITY) && (rub [ind [j]] >  COUENNE_INFINITY))))
	  continue;

	pairs.insert (std::pair <int, int> (ind [j], ind [k]));
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

  do {

    nCurTightened = 0;

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

	nCurTightened += combine (cs, n1, n2, ind1, ind2, sa1, sa2, a1, a2, clb, cub, l1, l2, u1, u2, -1);

      // clean sa1 and sa2

      for (int i=n1; i--;) sa1 [ind1 [i]] = 0.;
      for (int i=n2; i--;) sa2 [ind2 [i]] = 0.;
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
