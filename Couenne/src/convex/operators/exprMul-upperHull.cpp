/* $Id$ 
 *
 * Name:    exprMul-upperHull.cpp
 * Author:  Pietro Belotti
 * Purpose: generates upper envelope of a product
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneExprMul.hpp"
#include "CouennePrecisions.hpp"

namespace Couenne {

//#define DEBUG

int findIntersection (CouNumber  x0, CouNumber  y0,
		      CouNumber  x1, CouNumber  y1,
		      CouNumber *wl, CouNumber *wu,
		      CouNumber *xA, CouNumber *yA,
		      CouNumber *xB, CouNumber *yB);

int genMulCoeff (CouNumber x1, CouNumber y1, 
		  CouNumber x2, CouNumber y2, 
		  char whichUse,
		  CouNumber &cX, CouNumber &cY, CouNumber &cW);


// invert interval bounds and current point
inline void invertInterval (register double &l, register double &u, register double x) {

  register double tmp = l; 
  l = -u; 
  u = -tmp;

  x = -x;
}

void upperEnvHull (const CouenneCutGenerator *cg, OsiCuts &cs, 
		   int xi, CouNumber x0, CouNumber xl, CouNumber xu,
		   int yi, CouNumber y0, CouNumber yl, CouNumber yu,
		   int wi, CouNumber w0, CouNumber wl, CouNumber wu) {

  //
  // See forthcoming paper for explanation ;-)
  //

#ifdef DEBUG
  printf ("entering points: ===================================================\n");
  printf ("x [%d]: %9g\t[%9g\t%9g]\n", xi, x0, xl, xu);
  printf ("y [%d]: %9g\t[%9g\t%9g]\n", yi, y0, yl, yu);
  printf ("w [%d]: %9g\t[%9g\t%9g]\n", wi, w0, wl, wu);
#endif

  // Preprocess to reduce everything to a first-orthant problem

  if ((wl < 0 && wu > 0)) // nothing to tighten
    return;

  // project back into bounding box
  if (x0 < xl) x0 = xl;  if (x0 > xu) x0 = xu;
  if (y0 < yl) y0 = yl;  if (y0 > yu) y0 = yu;

  // preliminary bound tightening
  if (wl >= 0) {
    if        ((xl >= 0) || (yl >= 0) || (xl * yl < wl - COUENNE_EPS)) {
      if (xl < 0) xl = 0;
      if (yl < 0) yl = 0;
    } else if ((xu <= 0) || (yu <= 0) || (xu * yu < wl - COUENNE_EPS)) {
      if (xu > 0) xu = 0;
      if (yu > 0) yu = 0;
    }
  } else {
    if        ((xl >= 0) || (yu <= 0) || (xl * yu > wu + COUENNE_EPS)) {
      if (xl < 0) xl = 0;
      if (yu > 0) yu = 0;
    } else if ((xu <= 0) || (yl >= 0) || (xu * yl > wu + COUENNE_EPS)) {
      if (xu > 0) xu = 0;
      if (yl < 0) yl = 0;
    }
  }

  // Reduce

  bool 
    flipX = xl < 0,
    flipY = yl < 0,
    flipW = false;

  if (flipX && flipY) { // swap bounds on x and y

    invertInterval (xl,xu,x0);
    invertInterval (yl,yu,y0);

  } else if (flipX) { // swap bounds on x and w only

    invertInterval (xl,xu,x0);
    invertInterval (wl,wu,w0);

    flipW = true;

  } else if (flipY) { // swap bounds on y and w only

    invertInterval (yl,yu,y0);
    invertInterval (wl,wu,w0);

    flipW = true;
  }

#ifdef DEBUG
  printf ("reduced points:\n");
  printf ("x: %9g\t[%9g\t%9g]\n", x0, xl, xu);
  printf ("y: %9g\t[%9g\t%9g]\n", y0, yl, yu);
  printf ("w: %9g\t[%9g\t%9g]\n", w0, wl, wu);
#endif

  // Check whether lower and upper curve both contain bounding box of
  // x,y. If so, there is nothing to separate

  if (((xl*yl >= wl) &&  // b.box totally contained between two curves
       (xu*yu <= wu)) || //
      (x0*y0 >= w0))     // or current point below curve
    return;

  // Find intersections of halfline from origin

  CouNumber xLow, yLow, xUpp, yUpp;
  if (findIntersection (0., 0., x0, y0, &wl, &wu, &xLow, &yLow, &xUpp, &yUpp))
    return; // couldn't find proper point

#ifdef DEBUG
  printf ("intersections:\n");
  printf ("lower: %9g\t%9g\tprod %9g\n", xLow, yLow, xLow*yLow);
  printf ("upper: %9g\t%9g\tprod %9g\n", xUpp, yUpp, xUpp*yUpp);
#endif

  // Case 1: both are outside of bounding box, either NW or SE. McCormick's suffice.

  if ((xLow <= xl && yUpp >= yu) ||
      (yLow <= yl && xUpp >= xu))
    return;

  // There will be at least one cut. Define coefficients and rhs ---
  // will have to change them back if (flipX || flipY)

  CouNumber
    cX,  cY,  cW,  c0,  c0X,  c0Y,  c0W;

  if (xLow >= xl && xUpp <= xu &&
      yLow >= yl && yUpp <= yu) {

#ifdef DEBUG
    printf ("easy lifting:\n");
#endif

    // Case 2: both are inside. Easy lifting...
    if (genMulCoeff (xLow, yLow, xUpp, yUpp, 0, cX, cY, cW))
      return;

    c0X = cX * xLow;
    c0Y = cY * yLow;
    c0W = cW * wl;

  } else if (xLow >= xl && yLow >= yl && 
	     (xUpp > xu || yUpp > yu)) {

#ifdef DEBUG
    printf ("through lower, not upper:\n");
#endif

    // Case 3a and 3b: through lower curve, but not upper. 

    if (yUpp > yu) { // upper intersect is North; place it within box
      yUpp = yu;
      xUpp = wu / yu;
    } else {         //                    East
      xUpp = xu;
      yUpp = wu / xu;
    }

    // find intersection on low curve on half line through new point and (x0,y0)
    if ((findIntersection (xUpp, yUpp, x0, y0, &wl, NULL, &xLow, &yLow, NULL, NULL)) ||
	(xLow < xl || yLow < yl) ||                            // McCormick's suffice
	(genMulCoeff (xLow, yLow, xUpp, yUpp, 0, cX, cY, cW))) // Otherwise, lift inequality on lower point
      return;

    c0X = cX * xLow;
    c0Y = cY * yLow;
    c0W = cW * wl;

  } else if (xUpp <= xu && yUpp <= yu && 
	     (xLow < xl || yLow < yl)) {

#ifdef DEBUG
    printf ("through upper, not lower:\n");
#endif

    // Case 4a and 4b: viceversa (lift for validity)

    if (yLow < yl) { // upper intersect is South; place it within box
      yLow = yl;
      xLow = wl / yl;
    } else {         //                    West
      xLow = xl;
      yLow = wl / xl;
    }

    // find intersection on low curve on half line through new point and (x0,y0)
    if (findIntersection (xLow, yLow, x0, y0, NULL, &wu, NULL, NULL, &xUpp, &yUpp))
      return;

    if (xUpp > xu || yUpp > yu) // McCormick's suffice
      return;

    // Otherwise, lift inequality on UPPER point
    if (genMulCoeff (xLow, yLow, xUpp, yUpp, 1, cX, cY, cW))
      return;

    c0X = cX * xUpp;
    c0Y = cY * yUpp;
    c0W = cW * wu;

  } else if ((xLow < xl && xUpp > xu) ||
	     (yLow < yl && yUpp > yu)) {

#ifdef DEBUG
    printf ("between lower and upper:\n");
#endif

    // Case 5: both outside of bounding box, N and S or W and
    //         E. Separate both from lower and from upper

    if (yLow < yl) { // upper intersect is South; place it within box
      yLow = yl;      xLow = wl / yl;
      yUpp = yu;      xUpp = wu / yu;
    } else {         //                    West
      xLow = xl;      yLow = wl / xl;
      xUpp = xu;      yUpp = wu / xu;
    }

#ifdef DEBUG
    printf ("New intersections:\n");
    printf ("lower: %9g\t%9g\tprod %9g\n", xLow, yLow, xLow*yLow);
    printf ("upper: %9g\t%9g\tprod %9g\n", xUpp, yUpp, xUpp*yUpp);
#endif

    // Nothing to find. Just separate two inequalities at the same
    // point, just using different support
    //
    // if (genMulCoeff (xLow, yLow, xUpp, yUpp, 0, cX,  cY,  cW) ||
    //     genMulCoeff (xLow, yLow, xUpp, yUpp, 1, cXp, cYp, cWp))
    //   return;

    // A more clever way (courtesy of Andrew J. Miller): find the
    // intersect on the lower (upper) curve on the line through xLP
    // and the upper (lower) point

    CouNumber xLow2, yLow2, xUpp2, yUpp2;

    if ((findIntersection (xLow, yLow, x0, y0, NULL, &wu, NULL,   NULL,   &xUpp2, &yUpp2) || genMulCoeff (xLow, yLow, xUpp, yUpp, 0, cX, cY, cW)) &&
	(findIntersection (xUpp, yUpp, x0, y0, &wl, NULL, &xLow2, &yLow2, NULL,   NULL)   || genMulCoeff (xLow, yLow, xUpp, yUpp, 1, cX, cY, cW)))
      return;

#ifdef DEBUG
    printf ("coeffs: (%g,%g,%g)\n", // [(%g,%g,%g)]\n", 
	    cX,cY,cW);
#endif

    c0X = cX * xLow; //   c0Xp = cXp * xUpp;
    c0Y = cY * yLow; //   c0Yp = cYp * yUpp;
    c0W = cW * wl;   //   c0Wp = cWp * wu;  

//     twoIneqs = true;

  } else {

#ifdef DEBUG
    printf ("points are in a weird position:\n");
    printf ("lower: %9g\t%9g\tprod %9g\n", xLow, yLow, xLow*yLow);
    printf ("upper: %9g\t%9g\tprod %9g\n", xUpp, yUpp, xUpp*yUpp);
#endif

    return;
  }

  // Re-transform back into original variables

  if (flipX) {cX = -cX; c0X = -c0X;}
  if (flipY) {cY = -cY; c0Y = -c0Y;}
  if (flipW) {cW = -cW; c0W = -c0W;}

  c0 = c0X + c0Y + c0W;

#ifdef DEBUG
  printf ("there are cuts\n");
#endif

  //cg -> createCut (cs, alpha*wb + 2*wb/xt, sign, wi, alpha, yi, 1., xi, wb/(xt*xt));
  cg   -> createCut (cs, c0,  1, wi, cW,  yi, cY,  xi, cX);
}


// finds intersections of a parametric line (x,y) = (x0,y0) + t [(x1,y1) - (x0,y0)] 
// on curves xy = wl and xy = yu

int findIntersection (CouNumber  x0, CouNumber  y0,
		      CouNumber  x1, CouNumber  y1,
		      CouNumber *wl, CouNumber *wu,
		      CouNumber *xA, CouNumber *yA,
		      CouNumber *xB, CouNumber *yB) {

  // The parametric line is of the form
  //
  //  x = x0 + t (x1-x0)
  //  y = y0 + t (y1-y0)
  // 
  // and for that to satisfy xy = wl and xy = wu we must have
  //
  // x0 y0 + t [x0 (y1-y0) + y0 (x1-x0)] + t^2 (x1-x0) (y1-y0) = wl
  //                                                           = wu
  //
  // or a t^2 + b t + c - wl = 0 for proper values of a,b,c.
  //    a t^2 + b t + c - wu = 0
  //
  // Because of the way this procedure will be used, of the two
  // solutions found we must always use the minimum nonnegative one

  CouNumber
    a = (x1-x0) * (y1-y0),
    c = x0 * y0,
    b = x0*y1 + y0*x1 - 2*c; // x0 * (y1-y0) + y0 * (x1-x0)

  if (fabs (a) < COUENNE_EPS)
    return 1;

  // These are the solution to the above equation.

  CouNumber tL1, tL2, tU1, tU2, tL=0., tU=0.;

  if (wl) {
    tL1 = (- b - sqrt (b*b - 4*a*(c-*wl))) / (2*a);
    tL2 = (- b + sqrt (b*b - 4*a*(c-*wl))) / (2*a);
    //printf ("solutions L: %g %g (b2-4ac=%g)\n", tL1, tL2, b*b - 4*a*(c-*wl));
    tL = (tL1 < 0) ? tL2 : tL1;
  }

  if (wu) {
    tU1 = (- b - sqrt (b*b - 4*a*(c-*wu))) / (2*a);
    tU2 = (- b + sqrt (b*b - 4*a*(c-*wu))) / (2*a);
    //printf ("solutions U: %g %g (b2-4ac=%g)\n", tU1, tU2, b*b - 4*a*(c-*wu));
    tU = (tU1 < 0) ? tU2 : tU1;
  }

  if (xA) *xA = x0 + tL * (x1-x0);   if (yA) *yA = y0 + tL * (y1-y0);
  if (xB) *xB = x0 + tU * (x1-x0);   if (yB) *yB = y0 + tU * (y1-y0);

  return 0;
}


// generate coefficients for a plane through points (x1, y1, x1 y1)
// and (x2, y2, x2 y2) such that intersecting it with one of them (the
// first if whichUse==0, the second otherwise) gives a tangent to the
// curve xy = k.

int genMulCoeff (CouNumber x1, CouNumber y1, 
		  CouNumber x2, CouNumber y2, 
		  char whichUse,
		  CouNumber &cX, CouNumber &cY, CouNumber &cW) {

  // the x-y slope of this constraint must be tangent to a curve xy=k
  // at (xD,yD). Easy:

  CouNumber xD, yD, xO, yO;

  if (!whichUse) {
    xD = x1; xO = x2;
    yD = y1; yO = y2;
  } else {
    xD = x2; xO = x1;
    yD = y2; yO = y1;
  }

  cX = yD;
  cY = xD;
  //c0 = 2*xD*yD;

#ifdef DEBUG
    printf ("points: (%g,%g) (%g,%g), cW = (%g - %g) / %g = %g\n", 
	    xD,yD, xO,yO, 2*xD*yD, (cX*xO + cY*yO), (xO*yO - xD*yD),
	    (2*xD*yD - (cX*xO + cY*yO)) / (xO*yO - xD*yD));
#endif

  // Now the hard part: lift it so that it touches the other curve

  if (fabs (xO*yO - xD*xD) < COUENNE_EPS) 
    return 1; // no cut if the two points are on the same curve

  // should ALWAYS be negative
  cW = (2*xD*yD - (cX*xO + cY*yO)) / (xO*yO - xD*yD);

  return 0;
}

}
