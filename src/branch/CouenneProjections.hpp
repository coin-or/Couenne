/*
 *
 * Name:    projections.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tools for projecting points on lines/planes
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneProjections_hpp
#define CouenneProjections_hpp

#include <stdio.h>

#include "CouenneConfig.h"
#include "CouennePrecisions.hpp"

namespace Couenne {

/** Compute projection of point (x0, y0) on the segment defined by
 *  line ax + by + c <>= 0 (sign provided by parameter sign) and
 *  bounds [lb, ub] on x. Return distance from segment, 0 if satisfied
 */

COUENNELIB_EXPORT
CouNumber project (CouNumber a,   CouNumber b, CouNumber c,
		   CouNumber x0,  CouNumber y0,
		   CouNumber lb,  CouNumber ub,
		   int sign,
		   CouNumber *xp = NULL, CouNumber *yp = NULL);

/** Compute projection of point (x0, y0) on the segment defined by two
 *  points (x1,y1), (x2, y2) -- sign provided by parameter
 *  sign. Return distance from segment, 0 if on it.
 */

COUENNELIB_EXPORT
CouNumber projectSeg (CouNumber x0,  CouNumber y0,
		      CouNumber x1,  CouNumber y1,
		      CouNumber x2,  CouNumber y2,
		      int sign,
		      CouNumber *xp = NULL, CouNumber *yp = NULL);
}

#endif
