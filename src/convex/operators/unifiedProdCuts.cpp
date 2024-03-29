/*
 *
 * Name:    unifiedProdCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: unified convexification of products and divisions
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneFunTriplets.hpp"

using namespace Couenne;

/// Add cut around curve x*y=k

void contourCut (const CouenneCutGenerator *cg,
		 OsiCuts &cs,
		 CouNumber xp, CouNumber yp, // current point
		 CouNumber wb,               // bound on w
		 int sign,                   // is wb lower or upper?
		 CouNumber x0, CouNumber y0, // (allegedly) outside point
		 CouNumber x1, CouNumber y1, //              inside
		 int xi, int yi, int wi) {   // indices of the variables

  // Upper right corner of the bounding box of (x,y) is feasible,
  // the opposite corner is not, hence there is a cut violated by
  // (x0,y0).

  // If (x0,y0) is not in the same orthant as the contour in
  // question, move it in so that we can apply a Newton step to
  // find closest point on contour.

  int xsign = (x1 >= 0) ? 1 : -1, // define orthant where the "inside
      ysign = (y1 >= 0) ? 1 : -1; // point" lies

  if      (((xsign > 0) ? xp : -xp) <= COUENNE_EPS)
    if    (((ysign > 0) ? yp : -yp) <= COUENNE_EPS) {

      // opposite orthant, put in the right one where constraint is violated
      xp = yp = sqrt (fabs (wb))/2;
      if (xsign<0) xp = -xp;
      if (ysign<0) yp = -yp;
    }                                                // otherwise, must cross one axis only:
    else                                            {xp = sqrt (fabs(wb/yp)); if (xsign<0) xp=-xp;}//y
  else if (((ysign > 0) ? yp : -yp) <= COUENNE_EPS) {yp = sqrt (fabs(wb/xp)); if (ysign<0) yp=-yp;}//x

  // pt here describes a function of the form wb*x^(-1)
  kpowertriplet pt (-1, wb);

  CouNumber
    // tangent point closest to current point
    xt = powNewton (xp, yp, &pt),
    *lb = cg -> Problem () -> Lb (),
    *ub = cg -> Problem () -> Ub (),
    xL = lb [xi], xU = ub [xi],
    yL = lb [yi], yU = ub [yi];

  if (xt == 0.) // no tangents are possible
    return;

  // check if (xt,wb/xt) is outside of bounds. If so, project it back
  // into the bounding box
  if ((xt < xL) && (xL != 0.)) xt = xL;
  if ((xt > xU) && (xU != 0.)) xt = xU;

  if ((wb / xt < yL) && (yL != 0.)) xt = wb / yL;
  if ((wb / xt > yU) && (yU != 0.)) xt = wb / yU;

  // coefficient of w in the lifted cut
  CouNumber
    alpha = ((fabs (x1) < COUENNE_INFINITY) &&
	     (fabs (y1) < COUENNE_INFINITY) &&
	     (fabs (x1*y1 - wb) > 0.)) ?
    ((2*wb/xt - y1 - wb*x1 / (xt*xt)) / (x1*y1 - wb)) : 0;

  //printf ("+++++ %d %d %d. [%c] xp (%g,%g) wb %g out(%g,%g) in(%g,%g) --> [%g,%g] alpha %g\n",
  //xi, yi, wi, (sign<0) ? '-' : '+', xp, yp, wb, x0, y0, x1, y1, xt, wb/xt, alpha);

  if (alpha != 0)
    cg     -> createCut (cs, alpha*wb + 2*wb/xt, sign, wi, alpha, yi, 1., xi, wb/(xt*xt));
  else  cg -> createCut (cs,            2*wb/xt, sign,            yi, 1., xi, wb/(xt*xt));
}



// Unified procedure to create convexification cuts for an expression of the form w = x*y
void Couenne::unifiedProdCuts (const CouenneCutGenerator *cg, OsiCuts &cs,
			       int xi, CouNumber x0, CouNumber xl, CouNumber xu,
			       int yi, CouNumber y0, CouNumber yl, CouNumber yu,
			       int wi, CouNumber w0, CouNumber wl, CouNumber wu,
			       t_chg_bounds *chg, enum expression::auxSign sign) {

  bool cLX,  cRX,  cLY,  cRY,  cLW,  cRW =
       cLX = cRX = cLY = cRY = cLW = true;

  if (!(cg -> isFirst ()) && chg) {
    cLX= chg[xi].lower() != t_chg_bounds::UNCHANGED; cRX= chg[xi].upper() != t_chg_bounds::UNCHANGED;
    cLY= chg[yi].lower() != t_chg_bounds::UNCHANGED; cRY= chg[yi].upper() != t_chg_bounds::UNCHANGED;
    cLW= chg[wi].lower() != t_chg_bounds::UNCHANGED; cRW= chg[wi].upper() != t_chg_bounds::UNCHANGED;
  }

  // Add McCormick convexification cuts:
  //
  // 1) w >= yl x + xl y - yl xl
  // 2) w >= yu x + xu y - yu xu
  //
  // 3) w <= yl x + xu y - yl xu
  // 4) w <= yu x + xl y - yu xl
  //
  // These cuts are added if the corresponding bounds are finite

  if (sign != expression::AUX_LEQ) {
    if ((cLX || cLY) && is_boundbox_regular (yl, xl)) cg -> createCut (cs,yl*xl,-1,wi,-1.,xi,yl,yi,xl);
    if ((cRX || cRY) && is_boundbox_regular (yu, xu)) cg -> createCut (cs,yu*xu,-1,wi,-1.,xi,yu,yi,xu);
  }

  if (sign != expression::AUX_GEQ) {
    if ((cRX || cLY) && is_boundbox_regular (yl, xu)) cg -> createCut (cs,yl*xu,+1,wi,-1.,xi,yl,yi,xu);
    if ((cLX || cRY) && is_boundbox_regular (yu, xl)) cg -> createCut (cs,yu*xl,+1,wi,-1.,xi,yu,yi,xl);
  }

  // If w=xy and w >= l > 0 (resp. w <= u < 0) are "tight" bounds
  // (i.e. they are tighter than those obtained through propagation of
  // x and y's bounds), McCormick's convexification is not tight as
  // the surface has a curve contour at w=l (resp. w=u).
  //
  // McCormick rules induce a tangent to this contour at the bounds of
  // both variables, but it may be useful to add further cuts along
  // the contour to eliminate infeasible point (x0,y0,w0), which may
  // be in the convexification but out of the contour (on its "convex"
  // side, or "out of the belly").
  //
  // Suppose P (xt,l/xt) (resp. (xt,u/xt) is the point on the contour
  // closest to (x0,y0), found through a Newton method. The cut is
  // tangent to the contour in P and has the form
  //
  //        y - l/xt >= -l/(xt^2) (x-xt)   if xl*yl < l and xu*yu > l
  //        y - l/xt <= -l/(xt^2) (x-xt)   if xl*yl > l and xu*yu < l
  //
  // (resp. y - u/xt <= -u/(xt^2) (x-xt)   if xl*yu > u and xu*yl < u
  //        y - u/xt >= -u/(xt^2) (x-xt)   if xl*yu < u and xu*yl > u)
  //
  // These can be lifted to satisfy, at equality, the point
  // (xu,yu,wu=xu*yu) (resp. (xl,yl,wl=xl*yl)), where xl and xu are
  // lower and upper bound of x, etc.
  //
  //        alpha (w - l) + y - l/xt >= -l/(xt^2) (x-xt) ...
  //
  // where alpha is such that the relation holds at equality at the
  // point (xu,yu,xu*yu):
  //
  //    alpha = [-yu + l/xt - l/(xt^2)(xu-xt)] / (xu*yu - l)

  if (cg -> Problem () -> MultilinSep () == CouenneProblem::MulSepSimple ||
      fabs (wu - wl) < COUENNE_EPS) {

    if ((x0 > xl + COUENNE_EPS) && (y0 > yl + COUENNE_EPS) &&
	(x0 < xu + COUENNE_EPS) && (y0 < yu + COUENNE_EPS)) {

      if (cLW && (wl > 0) && (x0*y0 < wl) && (sign != expression::AUX_GEQ)) { // that is, if (x0,y0) is out of the contour

	CouNumber xyl = xl * yl;

	// first and third orthant
	if      ((xyl <  wl) && (xu*yu >= wl)) contourCut (cg,cs, x0,y0, wl, +1, xl,yl, xu,yu, xi,yi,wi);
	else if ((xyl >= wl) && (xu*yu <  wl)) contourCut (cg,cs, x0,y0, wl, -1, xu,yu, xl,yl, xi,yi,wi);
      }

      // Similarly for w <= u < 0

      if (cRW && (wu < 0) && (x0*y0 > wu) && (sign != expression::AUX_LEQ)) { // that is, if (x0,y0) is out of the contour

	CouNumber xuyl = xl * yu;

	// second and fourth orthant
	if      ((xuyl > wu) && (xl*yu <= wu)) contourCut (cg,cs, x0,y0, wu, +1, xu,yl, xl,yu, xi,yi,wi);
	else if ((xuyl <=wu) && (xl*yu >  wu)) contourCut (cg,cs, x0,y0, wu, -1, xl,yu, xu,yl, xi,yi,wi);
      }
    }
  } else
    if (cg -> Problem () -> MultilinSep () == CouenneProblem::MulSepTight)
      upperEnvHull (cg, cs,
		    xi, x0, xl, xu,
		    yi, y0, yl, yu,
		    wi, w0, wl, wu);
}
