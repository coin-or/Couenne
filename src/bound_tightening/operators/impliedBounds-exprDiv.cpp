/*
 *
 * Name:    impliedBounds-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bounds for division operators
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprDiv.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

/// implied bound processing for expression w = x/y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprDiv::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  bool resx, resy = resx = false;

  CouNumber
    wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
    wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  // y is a constant
  if (arglist_ [1] -> Type () == CONST) {

    int ind = arglist_ [0] -> Index ();

    if (ind < 0) {
      printf ("exprDiv::impliedBound: Error, w = c/d constants\n");
      exit (-1);
    }

    CouNumber c = arglist_ [1] -> Value ();

    if (fabs (c) < COUENNE_EPS) {
      printf ("exprDiv::impliedBound: Error, division by zero\n");
      exit (-1);
    }

    // a copy of exprMul::impliedBound for the case where y is a constant

    bool xInt = arglist_ [0] -> isInteger ();

    if (c > COUENNE_EPS) {

      if (updateBound (-1, l+ind, xInt ? ceil (wl*c - COUENNE_EPS) : (wl*c))) {
	resx = true;
	chg [ind].setLower (t_chg_bounds::CHANGED);
      }

      if (updateBound (+1, u+ind, xInt ? floor (wu*c + COUENNE_EPS) : (wu*c))) {
	resx = true;
	chg [ind].setUpper (t_chg_bounds::CHANGED);
      }
    }
    else if (c < - COUENNE_EPS) {

      if (updateBound (-1, l+ind, xInt ? ceil  (wu*c - COUENNE_EPS) : (wu*c))) {
	resx = true;
	chg [ind].setLower (t_chg_bounds::CHANGED);
      }

      if (updateBound (+1, u+ind, xInt ? floor (wl*c + COUENNE_EPS) : (wl*c))) {
	resx = true;
	chg [ind].setUpper (t_chg_bounds::CHANGED);
      }
    }
  } else {

    // deal with all other cases

    // Each bound on w is represented on the xy plane with two cones,
    // such that joining the extreme rays of both one obtains two
    // lines, or in other words, the second cone is obtained through
    // the transformation (x',y') = (-x,-y) applied to the first cone.
    //
    // Bounds can be tightened according to four different cases,
    // depending on which corner of the bounding box belongs to which
    // of the cones and on the sign of the bounds.
    //
    // Define wl <= w <= wu,
    //        xl <= x <= xu,
    //        yl <= y <= yu.
    //
    // Then the "tightenable" bounds are, depending on the corner:
    //
    //             _______________________________________________________
    //            |           w >= wl         |          w <= wu          |
    //            |___________________________|___________________________|
    //            |     l<0     |    l>0      |     u<0     |     u>0     |
    //            |_____________|_____________|_____________|_____________|
    //   Cone     |upper |lower |upper |lower |upper |lower |upper |lower |
    //            |      |      |      |      |      |      |      |      |
    // 1 xl,yl    | INF  |  -   |  xl  |  yl  |xl,yl?|  yl  |  yl  |  -   |
    // 2 xl,yu    |  yu  |  xl  |  -   | INF  |  -   |  yu  |  yu  |yu,xl?|
    // 3 xu,yl    |  xu  |  yl  | INF  |  -   |  yl  |  -   |yl,xu?|  yl  |
    // 4 xu,yu    |  -   | INF  |  yu  |  xu  |  yu  |xu,yu?|  -   |  yu  |
    //            |______|______|______|______|______|______|______|______|
    //
    // Where "INF" stands for "infeasible subproblem", "-" for
    // "nothing to improve", and the rest is improved (those with "?"
    // may improve).

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber x0 = 0,
      *xl = l + xi, *yl = l + yi,
      *xu = u + xi, *yu = u + yi;

    /*printf ("from              : w[%d] [%e %e], x%d [%e %e] / y%d [%e %e]",
      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);*/

    // avoid changing bounds if x is constant
    if (xi == -1)
      xl = xu = &x0;

    //////////// deal with lower bound of w=x/y /////////////////////////////////////////////

    // simple case wl = wu = 0

    if ((fabs (wl) < COUENNE_EPS) &&
	(fabs (wu) < COUENNE_EPS)) {

      resx = updateBound (-1, xl, 0.) || resx;
      resx = updateBound (+1, xl, 0.) || resx;
      return resx || resy;
    }

    bool resxL, resxU, resyL,
      resyU = resxL = resxU = resyL = false;

    // general case

    if        (wl < - COUENNE_EPS && wl > - COUENNE_INFINITY / 10) { // w >= wl, wl negative

      // point C: (xl,yl)

      resyL = ((*yl<0) && (*yl > *xl/wl + COUENNE_EPS) && updateBound (-1, yl,      0)) || resyL;
      resxL = ((*yl>0) && (*yu < *xl/wl)               && updateBound (-1, xl, *yu*wl)) || resxL; // new
      resyL = ((*yl>0) && (*xu < *yl*wl)               && updateBound (-1, yl, *xu/wl)) || resyL; // new

      // if ((*yl>0) && (*yl < *xl/wl)) { // point C violates x/y >= wl, down
      // 	resxL = updateBound (-1, xl, *yu*wl) || resxL; //
      // 	resyL = updateBound (-1, yl, *xu/wl) || resyL; //
      // }
      // point B: (xu,yu)
      // if ((*yu<0) && (*yu > *xu/wl)) { // point B violates x/y >= wl, down
      // 	resxU = updateBound (+1, xu, *yl*wl) || resxU;
      // 	resyU = updateBound (+1, yu, *xl/wl) || resyU;
      // }

      resyU = ((*yu>0) && (*yu < *xu/wl - COUENNE_EPS) && updateBound (+1, yu,      0)) || resyU;
      resxU = ((*yu<0) && (*yl > *xu/wl)               && updateBound (+1, xu, *yl*wl)) || resxU; // new
      resyU = ((*yu<0) && (*xl > *yu*wl)               && updateBound (+1, yu, *xl/wl)) || resyU; // new

    } else if (wl >   COUENNE_EPS) { // w >= wl, wl positive

      //

      resyL = ((*yl<0) && (*yl < *xl/wl) && updateBound (-1, yl, CoinMin (*xl/wl, 0.))) || resyL;
      resxL = ((*yl>0) && (*yl > *xl/wl) && updateBound (-1, xl, *yl*wl))               || resxL;
      //
      resyU = ((*yu>0) && (*yu > *xu/wl) && updateBound (+1, yu, CoinMax (*xu/wl, 0.))) || resyU;
      resxU = ((*yu<0) && (*yu < *xu/wl) && updateBound (+1, xu, *yu*wl))               || resxU;
    }

    //////////// deal with upper bound of w=x/y /////////////////////////////////////////////


    if        (wu >   COUENNE_EPS && wu < COUENNE_INFINITY / 10) { // w <= wu, wu positive

      //
      resyL = ((*yl<0) && (*yl > *xu/wu + COUENNE_EPS) && updateBound (-1, yl,      0)) || resyL;
      resxU = ((*yl>0) && (*xu > *yu*wu)               && updateBound (+1, xu, *yu*wu)) || resxU;
      resyL = ((*yl>0) && (*yl < *xl/wu)               && updateBound (-1, yl, *xl/wu)) || resyL;

      // if ((*yl>0) && (*yl < *xu/wu)) {
      // 	resxU = updateBound (+1, xu, *yu*wu) || resxU;
      // 	resyL = updateBound (-1, yl, *xl/wu) || resyL;
      // }

      // //

      // if ((*yu<0) && (*yu > *xl/wu)) {
      // 	resxL = updateBound (-1, xl, *yl*wu) || resxL;
      // 	resyU = updateBound (+1, yu, *xu/wu) || resyU;
      // }

      resyU = ((*yu>0) && (*yu < *xl/wu - COUENNE_EPS) && updateBound (+1, yu,      0)) || resyU;
      resxL = ((*yu<0) && (*xl < *yl*wu)               && updateBound (-1, xl, *yl*wu)) || resxL;
      resyU = ((*yu<0) && (*yu > *xu/wu)               && updateBound (+1, yu, *xu/wu)) || resyU;

    } else if (wu < - COUENNE_EPS) { // w <= wu, wu negative

      //

      resyL = ((*yl<0) && (*yl < *xu/wu) && updateBound (-1, yl, CoinMin (*xu/wu,0.)))  || resyL;//
      resxL = ((*yu<0) && (*yu < *xl/wu) && updateBound (-1, xl, *yu*wu))               || resxL;

      //

      resyU = ((*yu>0) && (*yu > *xl/wu) && updateBound (+1, yu, CoinMax (*xl/wu,0.)))  || resyU;
      resxU = ((*yl>0) && (*yl > *xu/wu) && updateBound (+1, xu, *yl*wu))               || resxU;
    }

    if (resxL) chg [xi].setLower(t_chg_bounds::CHANGED);
    if (resxU) chg [xi].setUpper(t_chg_bounds::CHANGED);
    if (resyL) chg [yi].setLower(t_chg_bounds::CHANGED);
    if (resyU) chg [yi].setUpper(t_chg_bounds::CHANGED);

    /*if (resx || resy)
      printf ("                 \ntightened division: w[%d] [%e %e], x%d [%e %e] / y%d [%e %e]\n",
	      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);
	      else printf ("                                                 \r");*/

    resx = resxL || resxU;
    resy = resyL || resyU;
  }

  bool
    xInt = arglist_ [0] -> isInteger (),
    yInt = arglist_ [1] -> isInteger ();

  if (resx && xInt) {
    int xi = arglist_ [0] -> Index ();
    assert (xi >= 0);
    u [xi] = floor (u [xi] + COUENNE_EPS);
    l [xi] = ceil  (l [xi] - COUENNE_EPS);
  }

  if (resy && yInt) {
    int yi = arglist_ [1] -> Index ();
    assert (yi >= 0);
    u [yi] = floor (u [yi] + COUENNE_EPS);
    l [yi] = ceil  (l [yi] - COUENNE_EPS);
  }

  return (resx || resy);
}
