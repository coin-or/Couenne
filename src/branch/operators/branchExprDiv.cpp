/*
 *
 * Name:    branchExprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: select branch for divisions
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprDiv.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneObject.hpp"

using namespace Couenne;

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprDiv::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts,
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {

  if (brDist) {free (brDist); brDist = NULL;} // clear it, computeMulBrDist will fill it

  int xi = arglist_ [0]        -> Index (),
      yi = arglist_ [1]        -> Index (),
      wi = obj -> Reference () -> Index ();

  assert ((xi >= 0) && (yi >= 0) && (wi >= 0));

  brpts  = (double *) realloc (brpts,  sizeof (double));

  // choosing branching variable and point is difficult, use
  // proportion in bound intervals

  CouNumber yl = info -> lower_    [yi],
            yu = info -> upper_    [yi],
            y0 = info -> solution_ [yi];

  // if [yl,yu] contains 0, create two nodes

  if ((yl < -COUENNE_EPS) &&
      (yu >  COUENNE_EPS)) {

    var = arglist_ [1];

    *brpts = 0.;

    way = (y0 > *brpts) ? TWO_RIGHT : TWO_LEFT;

    brDist = computeMulBrDist (info, wi, yi, xi, yi, brpts);

    return CoinMin (brDist [0], brDist [1]);
    /*return ((fabs (y0) < COUENNE_EPS) ? 1. :
	    fabs (info -> solution_ [xi] / y0 -
	    info -> solution_ [wi]));*/
  }

  // From now on, [yl,yu] may be unlimited in one sense only, and
  // interval does not contain 0.
  //
  // As convexification is still carried out by applying McCormick
  // rules to x=w*y (where original operator is w=x/y), try to get
  // closer to a situation where both y and w are bounded, if
  // necessary by branching on w.
  //
  // Branch first on y if unbounded, and then on w. As a result of
  // bound tightening, if both y and w are bounded, x will be, too.

  if ((yl < -COUENNE_INFINITY) || // and yu < 0
      (yu >  COUENNE_INFINITY)) { // and yl > 0

    var = arglist_ [1];

    // if y0 close to bounds, branch away from it
    if      (fabs (y0-yl) < COUENNE_NEAR_BOUND) *brpts = y0 + 1. + yl*10.;
    else if (fabs (y0-yu) < COUENNE_NEAR_BOUND) *brpts = y0 - 1. + yu*10.;
    else                                        *brpts = y0;

    way = (y0 > 0.) ? TWO_LEFT : TWO_RIGHT;

    brDist = computeMulBrDist (info, wi, yi, xi, yi, brpts);

    return CoinMin (brDist [0], brDist [1]);

    //return ((fabs (y0) < COUENNE_EPS) ? 1. :
    //fabs (info -> solution_ [xi] / y0 - info -> solution_ [wi]));
  }

  // y is bounded, and y0 nonzero; if w is unbounded, bound w first as
  // x will be too.

  CouNumber wl = info -> lower_    [wi],
            wu = info -> upper_    [wi],
            w0 = info -> solution_ [wi],
            x0 = info -> solution_ [xi];

  // w unbounded in >= one direction ///////////////////////////////////////

  if ((wl < -COUENNE_INFINITY) ||
      (wu >  COUENNE_INFINITY)) {

    var = obj -> Reference ();

    if ((wl < -COUENNE_INFINITY) &&
	(wu >  COUENNE_INFINITY)) { // unbounded in two directions

      CouNumber
	wreal = x0 / y0,
	wmin  = w0,
	wmax  = wreal; // assume (x0,y0,w0) is below w=x/y

      if (wreal < w0) { // but swap if (x0,y0,w0) is above w=x/y
	wmin = wreal;
	wmax = w0;
      }

      *brpts = wreal;
      way = (w0 < wreal) ? TWO_LEFT : TWO_RIGHT;

      brDist = computeMulBrDist (info, wi, yi, xi, wi, brpts);
      return CoinMin (brDist [0], brDist [1]);

    } else {

      // unbounded in one direction only, use two way branching

      // if y0 close to bounds, branch away from it
      if      (fabs (w0-wl) < COUENNE_NEAR_BOUND) *brpts = w0 + 1 + wl*10;
      else if (fabs (w0-wu) < COUENNE_NEAR_BOUND) *brpts = w0 - 1 + wu*10;
      else                                        *brpts = w0;

      way = (wl < - COUENNE_INFINITY) ? TWO_RIGHT : TWO_LEFT;

      brDist = computeMulBrDist (info, wi, yi, xi, wi, brpts);
      return CoinMin (brDist [0], brDist [1]);
    }
    //return ((fabs (y0) < COUENNE_EPS) ? 1. : fabs (x0/y0 - w0));
  }

  // w and y are bounded (and so is x). Choose between x, y, z
  // depending on intervals first and then to vicinity to midpoint
  CouNumber
    xl = info -> lower_ [xi],
    xu = info -> upper_ [xi],
    dx = xu-xl,
    dy = yu-yl,
    dw = wu-wl;

  // Check largest interval and choose branch variable accordingly.
  // Branching point depends on where the current point is, but for
  // now just focus on the width of the intervals

  way = TWO_RAND;

  if (dx > dy)
    if (dx > dw) {var = arglist_[0];      *brpts = (xl+xu)/2.; /*return fabs (x0-y0*w0);*/}
    else         {var = obj->Reference(); *brpts = (wl+wu)/2.; /*return fabs (w0-x0/y0);*/}
  else
    if (dy > dw) {var = arglist_[1];      *brpts = (yl+yu)/2.; /*return fabs (y0-x0/w0);*/}
    else         {var = obj->Reference(); *brpts = (wl+wu)/2.; /*return fabs (w0-x0/y0);*/}

  brDist = computeMulBrDist (info, wi, yi, xi, var -> Index (), brpts);
  return CoinMin (brDist [0], brDist [1]);
}
