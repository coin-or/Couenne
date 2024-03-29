/*
 *
 * Name:    computeMulBrDist.cpp
 * Author:  Pietro Belotti
 * Purpose: compute distance to new convexifications generated by branching on product
 *
 * (C) Carnegie-Mellon University, 2008-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"

#include "CouenneExprMul.hpp"
#include "CouenneFunTriplets.hpp"
#include "CouenneProjections.hpp"

namespace Couenne {

// compute distance from future convexifications in set \f$\{(x,y,w):
// w = xy\}\f$ with x,y,w bounded
double *computeMulBrDist (const OsiBranchingInformation *info,
			  int xi, int yi, int wi,
			  int brind, double *brpt, int nPts) {

  // use rule of thumb to compute distance: fix two of the three
  // variables and compute distance between current LP point and curve
  // w=x/y (z is the branching variable, x or y)
  //
  // 1) fix x,y: distances are w - x/y and ||(w-x/y, z-zb)||_2
  // 2) fix x,w:               y - x/w and ||(y-x/w, z-zb)||_2
  // 3) fix y,w:               x - y*w and ||(x-y*w, z-zb)||_2

  CouNumber
    x0 = info -> solution_ [xi], //xl = info -> lower_ [xi], xu = info -> upper_ [xi],
    y0 = info -> solution_ [yi], //yl = info -> lower_ [yi], yu = info -> upper_ [yi],
    w0 = info -> solution_ [wi]; //wl = info -> lower_ [wi], wu = info -> upper_ [wi];

  double *dist = (double *) malloc (2 * sizeof (double));

  // two main cases:
  //
  // 1) wi is the branching index: the bounding box is divided in two
  // by the rule w <= wb and w >= wb. Finding the distances from the
  // current point (x0,y0,w0) to the two semi-surfaces depends on
  // which side w0 stands.
  //
  // 2) xi or yi are the branching index: reduce to the problem of
  // finding the distance from a point (x0,y0,w0) to the same surface
  // within the (two) element of a partition of the bounding box.


  // case 1 ////////////////////////////////////////////////////////////////////
  //
  // Depending on whether
  //
  // a) w0 <=/>= wb
  // b) w0 <=/>= x0*y0
  // c) wb <=/>= 0
  //
  // there are eight (!) cases, each with a similar computation for
  // the two distances. We unify it below.

  // [for now give a simplified version with rough distance computation]

  if (brind == wi) {

    // easy implementation: for own side, w0 - x0 y0; for other side,
    // horiz/vert distance to curve x0 w0 = wb

    bool side = (((x0*y0 > *brpt) && (*brpt > 0)) ||
		 ((x0*y0 < *brpt) && (*brpt < 0)));

    dist [side ? 1 : 0] = CoinMax
      (fabs (w0 - x0*y0), CoinMin
       ((fabs (y0) > COUENNE_EPS) ? fabs (x0 - *brpt / y0) : 0.,
	(fabs (x0) > COUENNE_EPS) ? fabs (y0 - *brpt / x0) : 0.));

    dist [side ? 0 : 1] = fabs (w0 - x0*y0);
  }

  // case 2 ////////////////////////////////////////////////////////////////////

  else {

    CouNumber diff = info -> solution_ [brind] - *brpt;
    bool side = (diff > 0.);

    dist [side ? 0 : 1] = CoinMax (fabs (w0 - x0 * y0), fabs (diff));
    dist [side ? 1 : 0] = fabs (w0 - x0 *y0);
  }

  //dist [0] = dist [1] = 1;
  return dist;
}

}
