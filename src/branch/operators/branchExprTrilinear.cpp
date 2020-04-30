/* $Id$
 *
 * Name:    branchExprTrilinear.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for trilinear terms
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"

#include "CouenneExprTrilinear.hpp"
#include "CouenneFunTriplets.hpp"
#include "CouenneProjections.hpp"

using namespace Couenne;

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprTrilinear::selectBranch (const CouenneObject *obj,
				       const OsiBranchingInformation *info,
				       expression *&var,
				       double * &brpts, 
				       double * &brDist, // distance of current LP point 
				       			 // to new convexifications
				       int &way) {

  if (brDist) {free (brDist); brDist = NULL;} // clear it, computeMulBrDist will fill it

  int
    xi = arglist_ [0] -> Index (),
    yi = arglist_ [1] -> Index (),
    zi = arglist_ [2] -> Index ();

  assert ((xi >= 0) && (yi >= 0) && (zi >= 0));

  CouNumber 
    xl = info -> lower_     [xi], yl = info -> lower_     [yi], zl = info -> lower_     [zi],
    xu = info -> upper_     [xi], yu = info -> upper_     [yi], zu = info -> upper_     [zi];

  brpts  = (double *) realloc (brpts,      sizeof (double));
  brDist = (double *) realloc (brDist, 2 * sizeof (double));

  // case 0: term is fixed

  if ((fabs (xu - xl) < COUENNE_EPS) &&
      (fabs (xu - xl) < COUENNE_EPS) &&
      (fabs (xu - xl) < COUENNE_EPS)) {

    var = NULL;
    return 0.;
  }

  // a very simple branching scheme: 
  //
  //            if any bound interval is (-inf,+inf), break it by branching at zero
  // otherwise, if any bound interval is [a,+inf) or (-inf,b], 
  //              break it by branching at zero (if interval includes zero) or 
  //                                    at 2a-1 (resp. 2b+1)
  // otherwise, branch on the largest bound, in the middle 
  
  if ((xl < -COUENNE_INFINITY) && (xu > COUENNE_INFINITY)) {*brpts = 0.; brDist [0] = brDist [1] = 1.; var = arglist_ [0]; return 1.;}
  if ((yl < -COUENNE_INFINITY) && (yu > COUENNE_INFINITY)) {*brpts = 0.; brDist [0] = brDist [1] = 1.; var = arglist_ [1]; return 1.;}
  if ((zl < -COUENNE_INFINITY) && (zu > COUENNE_INFINITY)) {*brpts = 0.; brDist [0] = brDist [1] = 1.; var = arglist_ [2]; return 1.;}


#define SETBNDS(l,u,ind) {			\
\
  if (l < -COUENNE_INFINITY) {\
    if (u > 1.) {*brpts = 0.;               brDist [0] = brDist [1] = 1.; var = arglist_ [ind]; return 1.;}\
    else        {*brpts = 2*-fabs (u) - 1.; brDist [0] = brDist [1] = 1.; var = arglist_ [ind]; return 1.;}\
  }\
\
  if (u >  COUENNE_INFINITY) {\
    if (l < -1.) {*brpts = 0.;              brDist [0] = brDist [1] = 1.; var = arglist_ [ind]; return 1.;}\
    else         {*brpts = 2*fabs (u) + 1.; brDist [0] = brDist [1] = 1.; var = arglist_ [ind]; return 1.;}\
  }\
}

  SETBNDS (xl, xu, 0);
  SETBNDS (xl, xu, 1);
  SETBNDS (xl, xu, 2);

  // bounds all finite, choose largest bound interval width
  if      ((xu - xl > yu - yl) && (xu - xl > zu - zl)) {*brpts = .5 * (xl + xu); brDist [0] = brDist [1] = 1.; var = arglist_ [0]; return 1.;}
  else if ((yu - yl > zu - zl))                        {*brpts = .5 * (yl + yu); brDist [0] = brDist [1] = 1.; var = arglist_ [1]; return 1.;}
  else                                                 {*brpts = .5 * (zl + zu); brDist [0] = brDist [1] = 1.; var = arglist_ [2]; return 1.;}
}
