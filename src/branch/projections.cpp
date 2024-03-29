/*
 *
 * Name:    projections.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tools for projecting points on lines/planes
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneTypes.hpp"
#include "CouennePrecisions.hpp"
//#include "CouenneProjections.hpp"

namespace Couenne {

/*  compute projection of point (x0, y0) on the segment defined by
 *  line ax + by + c <>= 0 (sign provided by parameter sign) and
 *  bounds [lb, ub] on x. Return distance from segment, 0 if satisfied
 */

CouNumber project (CouNumber a, CouNumber b, CouNumber c,
		   CouNumber x0, CouNumber y0,
		   CouNumber lb, CouNumber ub, int sign,
		   CouNumber *xp, CouNumber *yp) {

  /* compute projection of (x0,y0) onto line ax+by+c=0 */
  CouNumber
    t  = - (a*x0 + b*y0 + c);

  /* projection coordinates */
  CouNumber xpr, ypr;

  /* does point lie on line? */
  if (fabs (t) < COUENNE_EPS) return 0.;

  /* check if point satisfies inequality */
  if      (sign > 0) {if (t < 0.) return 0.;}
  else if (sign < 0) {if (t > 0.) return 0.;}

  /* t corresponding to intersection point */
  t /= sqrt (a*a + b*b);

  /* compute projection coordinates */
  xpr = x0 + a*t;
  ypr = y0 + b*t;

  /* don't need sign any longer, take its absolute value */
  if (t < 0.) t = -t;

  /* if projected point is outside [lb,ub], set xp to closest bound
     and yp accordingly, and compute distance to (x0,y0) */
  if ((xpr < lb) || (xpr > ub)) {

    if      (xpr < lb) xpr = lb;
    else if (xpr > ub) xpr = ub;

    ypr = (- c - a * xpr) / b - y0;
    xpr -= x0;

    t = sqrt (xpr * xpr + ypr * ypr);
  }

  /* update output parameters */
  if (xp) *xp = xpr;
  if (yp) *yp = ypr;

  /* return distance */
  return t;
}


/** Compute projection of point (x0, y0) on the segment defined by two
 *  points (x1,y1), (x2, y2) -- sign provided by parameter
 *  sign. Return distance from segment, 0 if on it.
 */

CouNumber projectSeg (CouNumber x0,  CouNumber y0,
		      CouNumber x1,  CouNumber y1,
		      CouNumber x2,  CouNumber y2,
		      int sign,
		      CouNumber *xp, CouNumber *yp) {
  CouNumber
    dx = x2-x1,
    dy = y2-y1,
    a =  -dy,
    b =  dx,
    c = x1*dy - y1*dx;

  return project (a, b, c, x0, y0, x1, x2, sign, xp, yp);
}

}

/*
int main (int argc, char **argv) {

  CouNumber
    a = atof (argv [1]),
    b = atof (argv [2]),
    c = atof (argv [3]),
    x0 = atof (argv [4]),
    y0 = atof (argv [5]),
    lb = atof (argv [6]),
    ub = atof (argv [7]);

  int sign = atoi (argv [8]);

  char sig = (sign < 0) ? '<' : (sign > 0) ? '>' : ' ';

  printf ("projecting (%.3f,%.3f) on %.3f x + %.3f y + %.3f %c= 0, xp in [%.3f,%.3f]\n",
	  x0, y0, a, b, c, sig, lb, ub);

  printf (" ==> distance is %.4f\n",
	  project (a, b, c, x0, y0, lb, ub, sign));

}
*/
