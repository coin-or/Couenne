/*
 *
 * Name:    branchExprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for powers
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneExprPow.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneProjections.hpp"
#include "CouenneFunTriplets.hpp"

using namespace Couenne;

/// generic approach for negative powers (used by expr{Pow,Inv}::selectBranch())
CouNumber negPowSelectBranch (const CouenneObject *obj,
			      const OsiBranchingInformation *info,
			      double * &brpts,
			      double * &brDist, // distance of current LP
				                // point to new convexifications
			      int &way,
			      CouNumber k,
			      CouNumber x0, CouNumber y0,
			      CouNumber l,  CouNumber u);


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprPow::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts,
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {

  // return branching point and branching policies for an expression
  // of the form x^k

  var = arglist_ [0];

  int
    ind = arglist_ [0]        -> Index (),
    wi  = obj -> Reference () -> Index ();

  assert ((ind >= 0) && (wi >= 0) && (arglist_ [1] -> Type () == CONST));

  double k = arglist_ [1] -> Value ();

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  /*printf ("selbra for "); print ();
  printf ("%g [%g,%g] -> %g [%g,%g]\n",
	  x0, l, u, y0,
	  info -> lower_    [wi],
	  info -> upper_    [wi]);*/

  if      (x0 < l) x0 = l;
  else if (x0 > u) x0 = u;

  // bounds coincide (happens within setupList)
  if (fabs (u-l) < COUENNE_EPS) {
    brpts = (CouNumber *) realloc (brpts, sizeof (CouNumber));
    *brpts = 0.5*(l+u);
    way = TWO_RAND;
    brDist = (double *) realloc (brDist, 2 * sizeof (double));
    return (brDist [0] = brDist [1] = fabs (y0 - (issignpower_ ? COUENNE_sign(x0) * pow(fabs(x0), k) : pow (x0, k))));
  }

  assert (!issignpower_ || k > 0);

  // case 1: k negative, resort to method similar to exprInv:: ///////////////////////////////
  if (k<0) return negPowSelectBranch (obj, info, brpts, brDist, way, k, x0, y0, l, u);

  brDist = (double *) realloc (brDist, 2 * sizeof (double));

  int intk = 0;

  bool isInt    =                        fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (k != 0.) && fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS;

  // case 2: k is positive and even /////////////////////////////////////////////////////////

  if (isInt && !((1 == intk % 2) || issignpower_)) {

    if ((l < - COUENNE_INFINITY) &&
	(u >   COUENNE_INFINITY)) { // infinite bounds

      if (y0 < pow (x0, k)) { // good side

	brpts = (double *) realloc (brpts, sizeof (double));
	powertriplet pt (k);
	*brpts = powNewton (x0, y0, &pt);

	way = (x0 > 0) ? TWO_RIGHT : TWO_LEFT; // priority to same side as x0

	x0 -= *brpts;
	y0 -= pow (*brpts, k);

	return (brDist [0] = brDist [1] = sqrt (x0*x0 + y0*y0)); // exact distance
      }

      // on the bad side /////////////////////

      // TODO: restore when we can do three-way branching

      /*double xp = pow (y0, 1./k);

      brpts = (double *) realloc (brpts, 2 * sizeof (double));
      brpts [0] = 0.5 * (x0 - xp);
      brpts [1] = 0.5 * (x0 + xp);

      way = THREE_CENTER;

      return CoinMin (x0 - brpts [0],
		      CoinMin (brpts [1] - x0,
			       projectSeg (x0, y0,
					   brpts [0], pow (brpts [0], k),
					   brpts [1], pow (brpts [1], k), -1)));*/

      // in the meantime, branch on current point (next branch won't
      // have unbounded x)

      brpts = (double *) realloc (brpts, sizeof (double));
      *brpts = x0;
      way = TWO_RAND;

      return (brDist [0] = brDist [1] = fabs (y0 - pow (x0,k)));

      // no bounds on x
      /*      double alpha = pow ((y0 + pow (x0, k))/2, 1./k),
             yroot = pow (y0, 1./k);
      brpts = (double *) realloc (brpts, 2 * sizeof (double));
      double lambdaL = (-x0 / yroot), lambdaR = 0;
      if (lambdaL < 0) {
	lambdaR = -lambdaL;
	lambdaL = 0;
      }
      CouNumber // approx distance
	b0 = brpts [0] = -alpha + lambdaL * (alpha - yroot),
	b1 = brpts [1] =  alpha + lambdaR * (yroot - alpha);
      way = THREE_CENTER;
      return CoinMin (projectSeg (x0, y0, b0, pow (b0, k), b1, pow (b1, k), -1),
		      CoinMin (x0 - b0, b1 - x0));      */
    }

    // at least one bound is finite //////////////////////////////////////////////

    brpts = (double *) realloc (brpts, sizeof (double));

    if (l < -COUENNE_INFINITY) {

      // if y0 is huge, try to make it close to 0
      *brpts = obj -> midInterval (-safe_pow (y0, 1. / k, issignpower_), l, u, info);
      way = TWO_RIGHT;

      //printf ("  ----> brptPow %g\n", *brpts);

      return CoinMin (brDist [0] = x0 - *brpts,
		      brDist [1] = projectSeg
		      (x0,     y0,
		       *brpts, safe_pow (*brpts, k, issignpower_),
		       u,      safe_pow (u, k, issignpower_), -1));
    }

    if (u >  COUENNE_INFINITY) {

      // if y0 is huge, try to make it close to 0
      //*brpts = CoinMin (safe_pow (y0, 1. / k), COU_MAX_COEFF / k);
      *brpts = obj -> midInterval (safe_pow (y0, 1. / k, issignpower_), l, u, info);
      way = TWO_LEFT;

      //printf ("  ----> brptPow %g\n", *brpts);

      return CoinMin (brDist [1] = *brpts - x0,
		      brDist [0] = projectSeg
		      (x0,     y0,
		       l,      safe_pow (l,      k, issignpower_),
		       *brpts, safe_pow (*brpts, k, issignpower_), -1));
    }

    // both bounds are finite ///////////////////////////////////////////////

    powertriplet ft (k);
    //*brpts = maxHeight (&ft, x0, y0, l, u);
    *brpts = obj -> getBrPoint (&ft, x0, l, u, info);

    way = (x0 < *brpts) ? TWO_LEFT : TWO_RIGHT;

    //    w -> print (); printf (" := "); dynamic_cast <exprAux *> (w) -> Image () -> print ();
    /*print ();
    printf (" (%g,%g) [%g,%g] --> %g, distances = %g,%g\n",
	    x0, y0, l, u, *brpts,
	    projectSeg (x0, y0, l,      safe_pow (l,k), *brpts, safe_pow (*brpts,k), 0),
	    projectSeg (x0, y0, *brpts, safe_pow (*brpts,k),      u, safe_pow (u,k), 0));*/

    // Min area  -- exact distance
    CouNumber retval =
      CoinMin (brDist [0] = projectSeg (x0, y0, l,      safe_pow (l,      k, issignpower_), *brpts, safe_pow (*brpts, k, issignpower_), 0),
	       brDist [1] = projectSeg (x0, y0, *brpts, safe_pow (*brpts, k, issignpower_),      u, safe_pow (u,      k, issignpower_), 0));

    //printf ("returning %.10e\n", retval);

    return retval;

    /*      brpts = (double *) realloc (brpts, sizeof (double));
     *brpts = midInterval (x0, l, u, info);
     way =
     (l < - COUENNE_INFINITY) ? TWO_RIGHT :
     (u >   COUENNE_INFINITY) ? TWO_LEFT  : TWO_RAND;

     return fabs (y0 - pow (x0, k)); // not an exact measure
    */
  }

  // from here on, we use two-way branch

  brpts = (double *) realloc (brpts, sizeof (double));
  *brpts = x0; // just in case none of the ifs below is satisfied...
  CouNumber pow0 = issignpower_ ? COUENNE_sign(x0) * pow(fabs(x0), k) : pow (x0, k);

  // case 3: k>1 and odd or signpower ////////////////////////////////////////////////////////////

  if ((isInt && (1 == intk % 2)) || issignpower_) {

    way = (x0 > 0.) ? TWO_RIGHT : TWO_LEFT;

    if (((l < - COUENNE_INFINITY) && (u > COUENNE_INFINITY)) || // [-inf,+inf[
	((l < - COUENNE_INFINITY) && (y0 < pow0))            ||
	((u >   COUENNE_INFINITY) && (y0 > pow0))) {

      if (((y0 > 0) && (y0 < pow0)) ||
	  ((y0 < 0) && (y0 > pow0))) {

	*brpts = 0;
	return (brDist [0] = brDist [1] = fabs (pow0 - y0));

      } else {

	*brpts = COUENNE_sign(y0) * pow(fabs(y0), 1./k);

	return (brDist [0] = brDist [1] = (y0 > 0) ? // approx distance
		projectSeg (x0, y0, x0, CoinMax (pow0, 0.), *brpts, y0, 0) :
		projectSeg (x0, y0, x0, CoinMin (pow0, 0.), *brpts, y0, 0));
      }
    }

    // otherwise, on the side of the current point the convexification
    // is bounded.

    powertriplet pt (k, issignpower_);
    *brpts = obj -> getBrPoint (&pt, x0, l, u, info);

    // in min-area and balanced strategy, point returned is
    // positive. Put the right sign

    if (y0 < pow0)
      *brpts = -*brpts;

    if (*brpts > x0) {

      brDist [0] = y0 - safe_pow (*brpts, k, issignpower_);
      brDist [1] = sqrt (brDist [0] * brDist [0] + (x0 - *brpts) * (x0 - *brpts));

    } else {

      brDist [1] = y0 - safe_pow (*brpts, k, issignpower_);
      brDist [0] = sqrt (brDist [1] * brDist [1] + (x0 - *brpts) * (x0 - *brpts));
    }

    // otherwise, the convexification is surely bounded.
    //
    // apply minarea

    /*if (u > l + COUENNE_EPS) {
      *brpts = safe_pow ((safe_pow (u, k) - safe_pow (l,k)) / (k* (u-l)), 1./(k-1.));
      if (u<0)
	*brpts = -*brpts;
	}*/
  }

  // case 4: k positive, in ]0,1[ and 1/k is integer and odd ////////////////////////

  if (isInvInt && (intk % 2) && !issignpower_) {

    way = (x0 > 0.) ? TWO_RIGHT : TWO_LEFT;

    if (((l < - COUENNE_INFINITY) && (u > COUENNE_INFINITY)) || // ]-inf,+inf[
	((l < - COUENNE_INFINITY) && (y0 < pow0))            ||
	((u >   COUENNE_INFINITY) && (y0 > pow0))) {

      if (((x0 > 0.) && (y0 > pow0)) ||
	  ((x0 < 0.) && (y0 < pow0))) { // in the same orthant as the
					// curve (i.e. first or third
					// orthant)

	*brpts = 0.;
	return (brDist [0] = brDist [1] = fabs (pow0 - y0));

      } else {

	*brpts = x0; // safe for any x0 as (l != x0 != u)

	//*brpts = obj -> midInterval (x0, l, u, info);

	return (brDist [0] = brDist [1] = (x0 > 0) ? // approx distance
	  projectSeg (x0,y0,x0, pow0, CoinMax (0., pow (y0, 1./k)), y0, 0) :
	  projectSeg (x0,y0,x0, pow0, CoinMin (0., pow (y0, 1./k)), y0, 0));
      }
    }

    // otherwise, on the side of the current point the convexification
    // is bounded.

    powertriplet pt (k);
    *brpts = obj -> getBrPoint (&pt, x0, l, u, info);

    // in min-area and balanced strategy, point returned is
    // positive. Put the right sign

    if      (y0 > pow0) *brpts = -*brpts;
    else if (y0 < pow0) *brpts = -*brpts;

    if (*brpts > x0) {

      brDist [0] = y0 - safe_pow (*brpts, k);
      brDist [1] = sqrt (brDist [0] * brDist [0] + (x0 - *brpts) * (x0 - *brpts));

    } else {

      brDist [1] = y0 - safe_pow (*brpts, k);
      brDist [0] = sqrt (brDist [1] * brDist [1] + (x0 - *brpts) * (x0 - *brpts));
    }

    return CoinMax (brDist [0], brDist [1]);

    /*if (u > l + COUENNE_EPS) {
      *brpts = safe_pow ((safe_pow (u, k) - safe_pow (l,k)) / (k* (u-l)), 1./(k-1.));
      if (u<0)
	*brpts = -*brpts;
    } else *brpts = 0.5*(l+u);*/
  }

  if (k>1) { // case 5: k>1, but not integer //////////////////////////////////////

    if (y0 < pow0) { // on the good side, i.e. out of the convex
		     // side. We don't care if u is infinity

      powertriplet pt (k, issignpower_);
      *brpts = powNewton (x0, y0, &pt);

      way = TWO_LEFT;

      x0 -= *brpts;
      y0 -= pow (*brpts, k);

      return (brDist [0] = brDist [1] = sqrt (x0*x0 + y0*y0));

    } else { // on the bad (concave) side, above curve

      // as a rule of thumb, take the x coordinate of the midpoint of
      // horizontal segment between current point and curve

      if (obj -> Strategy () == CouenneObject::MID_INTERVAL) {
	assert (!issignpower_);
	*brpts = 0.5 * (x0 + pow (y0, 1. / k));
      } else {
	powertriplet pt (k, issignpower_);
	*brpts = obj -> getBrPoint (&pt, x0, l, u, info);
      }

      way = TWO_LEFT;

      if (l < 0. && !issignpower_) l = 0.;

      CouNumber
	powbpt = issignpower_ ? COUENNE_sign(*brpts) * pow (fabs (*brpts), k) : pow (*brpts, k),
	projL  = projectSeg (x0, y0, l, issignpower_ ? COUENNE_sign(u) * pow(fabs(u),k) : pow (l, k), *brpts, powbpt, -1);

      return (brDist[0] = brDist[1] = (u > COUENNE_INFINITY) ?
	CoinMin (projL, *brpts - x0) :
	CoinMin (projL, projectSeg (x0,y0, *brpts, powbpt, u, pow(u,k), -1)));
    }

  } else { // case 6: 0<k<1 //////////////////////////////////////////////////////////////////

    if (y0 < pow0) { // on the bad side, below

      // same rule of thumb as above, take the x coordinate of the
      // midpoint of horizontal segment between current point and
      // curve

      if (obj -> Strategy () == CouenneObject::MID_INTERVAL)
	*brpts = 0.5 * (x0 + pow (y0, 1. / k));
      else {
	powertriplet pt (k);
	*brpts = obj -> getBrPoint (&pt, x0, l, u, info);
      }

      //*brpts = 0.5 * (x0 + pow (x0, 1. / k));
      way = TWO_LEFT;

      if (l < 0 && !issignpower_) l = 0.;

      CouNumber
	powbpt = pow (*brpts, k),
	projL  = projectSeg (x0, y0, l, pow (l, k), *brpts, powbpt, +1);

      return (brDist[0] = brDist[1] = (u > COUENNE_INFINITY) ?
	CoinMin (projL, powbpt - y0) :
	CoinMin (projL, projectSeg (x0, y0, *brpts, powbpt, u, pow (u,k), +1)));

    } else { // on the convex side. We don't care if u is infinite

      powertriplet pt (k);
      *brpts = powNewton (x0, y0, &pt);

      way = TWO_LEFT;

      x0 -= *brpts;
      y0 -= pow (*brpts, k);

      return (brDist [0] = brDist [1] = sqrt (x0*x0 + y0*y0));
    }
  }

  // failsafe: return null, so that CouenneObject picks the default
  // variable/branchingpoint for the expression

  assert (false);

  var = NULL;
  return 0.;
}
