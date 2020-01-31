/* $Id$
 *
 * Name:    CouenneLPtightenBounds.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tighten LP bounds on all variables (including continuous)
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

namespace Couenne {

// Tighten bounds - lightweight. Returns -1 if infeasible, otherwise
// number of variables tightened.
template <class T> 
int CouenneSolverInterface<T>::tightenBounds (int lightweight) {

  if (!(cutgen_ -> enableLpImpliedBounds ()))
    return 0;

  int 
    ncols = T::getNumCols (),
    nTightened;

  double 
    *oldLower = new double [ncols],
    *oldUpper = new double [ncols];

  CoinCopyN (T::getColLower (), ncols, oldLower);
  CoinCopyN (T::getColUpper (), ncols, oldUpper);

//   printf ("-------- BOUNDS BEFORE ------------\n  ");
//   int j=0;
//   for (int i=0; i < ncols; i++) {
//     printf("x_%03d [%+15.8g %+15.8g] ", i, oldLower [i], oldUpper [i]);
//     if (!(++j % 6)) printf ("\n  ");
//   }
//   if (j % 6) printf ("\n");

  nTightened = tightenBoundsCLP (lightweight);

  if (nTightened < 0)
    return nTightened;

//   printf ("-------- BOUNDS DURING ------------\n  ");
//   j=0;
//   for (int i=0; i < ncols; i++) {
//     printf("x_%03d [%+15.8g %+15.8g] ", i, getColLower () [i], getColUpper () [i]);
//     if (!(++j % 6)) printf ("\n  ");
//   }
//   if (j % 6) printf ("\n");

  if (nTightened > 0) {

    // something was tightened. Run an extra btCore "por si las
    // moscas" (just in case)

    const double 
      *newLower = T::getColLower (),
      *newUpper = T::getColUpper ();

    t_chg_bounds *chgd = new t_chg_bounds [ncols];

    for (int i=0; i<ncols; i++) {
      if (newLower [i] > oldLower [i] + COUENNE_EPS) chgd [i].setLower (t_chg_bounds::CHANGED);
      if (newUpper [i] < oldUpper [i] - COUENNE_EPS) chgd [i].setUpper (t_chg_bounds::CHANGED);
    }

    cutgen_ -> Problem () -> domain () -> push (ncols, NULL, newLower, newUpper);

    if (!(cutgen_ -> Problem () -> btCore (chgd))) // infeasible
      nTightened = -1;

    else {

      const double 
	*newerLower = cutgen_ -> Problem () -> Lb (),
	*newerUpper = cutgen_ -> Problem () -> Ub ();

      for (int i=0; i<ncols; i++) {

	if (newerLower [i] > newLower [i] + COUENNE_EPS) {
	  T::setColLower (i, newerLower [i]);
	  if (newLower [i] < oldLower [i] + COUENNE_EPS) nTightened++; // extra tightening
	}

      	if (newerUpper [i] < newUpper [i] - COUENNE_EPS) {
	  T::setColUpper (i, newerUpper [i]);
	  if (newUpper [i] > oldUpper [i] - COUENNE_EPS) nTightened++; // extra tightening
	}
      }
    }

//     const double 
//       *newerLower = cutgen_ -> Problem () -> Lb (),
//       *newerUpper = cutgen_ -> Problem () -> Ub ();

//     printf ("-------- BOUNDS AFTER ------------\n  ");
//     for (int i=0; i < ncols; i++) {
//       printf("x_%03d [%+15.8g %+15.8g] ", i, newerLower [i], newerUpper [i]);
//       if (!(++j % 6)) printf ("\n  ");
//     }
//     if (j % 6) printf ("\n");

    cutgen_ -> Problem () -> domain () -> pop ();

    delete [] chgd;
  }

  delete [] oldLower;
  delete [] oldUpper;

  return nTightened;
}

}
