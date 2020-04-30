/* $Id$
 *
 * Name:    branchExprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch suggestion for exprAbs
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <math.h>

// apparently unknown to OS
#ifndef M_SQRT2l
#define M_SQRT2l (sqrt (2.))
#endif

#include "CoinHelperFunctions.hpp"
#include "CouenneExprAbs.hpp"
#include "CouenneObject.hpp"

using namespace Couenne;

//static const double sqrt_2 = sqrt (2.);

/// set up branching object by evaluating branching points for each
/// expression's arguments. For an exprAbs, simply branch at zero.
CouNumber exprAbs::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression * &var,
				 double * &brpts,
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {
  var = argument_;

  int ind = var -> Index ();

  assert ((ind >= 0) && (obj -> Reference () -> Index () >= 0));

  CouNumber x0 = info -> solution_ [ind],
            y0 = info -> solution_ [obj -> Reference () -> Index ()];

  brpts = (double *) realloc (brpts, sizeof (double));

  // the best branching point for |x| is 0, as the two subproblems
  // will have exact convexifications (lines)
  *brpts = 0.;

  way = TWO_RAND; // don't care which subtree to visit first

  // no need to compute two distances for pseudocost, as this object
  // will only branch once...

  brDist = (double *) realloc (brDist, 2 * sizeof (double));

  assert ((y0 >=  x0 - COUENNE_EPS) && 
	  (y0 >= -x0 - COUENNE_EPS));

  brDist [0] = (x0 + y0) / M_SQRT2l;
  brDist [1] = (y0 - x0) / M_SQRT2l;

  // exact distance between current point and the two subsequent
  // convexifications
  return CoinMin (brDist [0], brDist [1]);
}
