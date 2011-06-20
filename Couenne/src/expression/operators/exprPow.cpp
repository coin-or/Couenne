/* $Id$
 *
 * Name:    exprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <math.h>
#include <assert.h>

#include "CouennePrecisions.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprLog.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneProblem.hpp"

#include "CouenneConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

/// simplify power f(x) ^ g(x)

expression *exprPow::simplify () {

  exprOp:: simplify ();

  if ((*arglist_) -> Type () == CONST) { // expr = c1 ^ g(x)

    CouNumber c0 = (*arglist_) -> Value ();

    if (arglist_ [1] -> Type () == CONST) { // expr = c1 ^ c2

      CouNumber c1 = arglist_ [1] -> Value ();

      delete arglist_ [0]; 
      delete arglist_ [1];

      arglist_ [0] = arglist_ [1] = NULL;

      return new exprConst (pow (c0, c1));
    }
    else 
      if (fabs (c0) <= COUENNE_EPS_SIMPL) 
	return new exprConst (0.);
  }
  else // only need to check if g(x) == 0

    if (arglist_ [1] -> Type () == CONST) {

      CouNumber expon = arglist_ [1] -> Value ();

      if (fabs (expon) <= COUENNE_EPS_SIMPL) // expr = x ^ 0 = 1
	return new exprConst (1.);

      else if (fabs (expon - 1.) <= COUENNE_EPS_SIMPL) { // expr = x ^ 1 = x

	delete arglist_ [1];
	expression *ret = arglist_ [0];
	arglist_ [0] = arglist_ [1] = NULL;
	return ret;
      }

      else if (fabs (expon + 1.) <= COUENNE_EPS_SIMPL) { // expr = x ^ -1 = 1/x

	delete arglist_ [1];
	expression *ret = new exprInv (arglist_ [0]);
	arglist_ [0] = arglist_ [1] = NULL;
	return ret;
      }

      //
      // x^k = x for x binary. Too bad we don't know bounds yet, so the code below will give segfault

      //       // is it an integer variable with bounds [-1,0] or [0,1]
      //       else if ((arglist_ [0] -> Type () == VAR) && (arglist_ [0] -> isDefinedInteger ())) {

      // 	CouNumber lb, ub;
      // 	arglist_ [0] -> getBounds (lb, ub);

      // 	if ((fabs (lb)      < COUENNE_EPS) &&
      // 	    (fabs (ub - 1.) < COUENNE_EPS)) {  // {0,1}

      // 	  delete arglist_ [1];
      // 	  expression *ret = arglist_ [0];
      // 	  arglist_ [0] = arglist_ [1] = NULL;
      // 	  return ret;

      // 	} else if ((fabs (lb + 1.) < COUENNE_EPS) &&
      // 		   (fabs (ub)      < COUENNE_EPS)) { // {-1,0}

      // 	  delete arglist_ [1];
      // 	  expression *ret = new exprOpp (arglist_ [0]);
      // 	  arglist_ [0] = arglist_ [1] = NULL;
      // 	  return ret;
      // 	}
      //       }

    }

  return NULL;
}


/// differentiate power of expressions

expression *exprPow::differentiate (int index) {

  if (!(arglist_ [0] -> dependsOn (index))  &&
      !(arglist_ [1] -> dependsOn (index)))
    return new exprConst (0.);

  if (arglist_ [0] -> Type () == CONST) { // k^f(x), k constant

    CouNumber base = arglist_ [0] -> Value ();

    if (base == 0.)
      return new exprConst (0.);

    return new exprMul (new exprConst (log (base)),
			new exprMul (new exprPow (new exprConst (base), 
						  arglist_ [1] -> clone ()),
				     arglist_ [1] -> differentiate (index)));

  } else if (arglist_ [1] -> Type () == CONST) { // f(x)^k, k constant

    CouNumber exponent = arglist_ [1] -> Value ();

    return new exprMul (new exprConst (exponent),
			new exprMul (new exprPow (arglist_ [0] -> clone (),
						  new exprConst (exponent - 1.)),
				     arglist_ [0] -> differentiate (index)));
  }

  // all other cases: f(x)^g(x)

  expression **alm  = new expression * [2];
  expression **alp  = new expression * [2];
  expression **als  = new expression * [2];
  expression **alm1 = new expression * [2];
  expression **alm2 = new expression * [2];
  expression **ald  = new expression * [2];

  alp [0] = new exprClone (arglist_ [0]);
  alp [1] = new exprClone (arglist_ [1]);

  alm [0] = new exprPow (alp, 2);

  alm1 [0] = arglist_ [1] -> differentiate (index);
  alm1 [1] = new exprLog (new exprClone (arglist_ [0]));

  als [0] = new exprMul (alm1, 2);

  ald [0] = new exprClone (arglist_ [1]);
  ald [1] = new exprClone (arglist_ [0]);

  alm2 [0] = new exprDiv (ald, 2);
  alm2 [1] = arglist_ [0] -> differentiate (index);

  als [1] = new exprMul (alm2, 2);

  alm [1] = new exprSum (als, 2);

  return new exprMul (alm, 2);
}


/// get a measure of "how linear" the expression is:
///
/// ZERO      = 0: a zero
/// CONSTANT  = 1: a constant
/// LINEAR    = 2: linear
/// QUADRATIC = 3: quadratic
/// NONLINER  = 4: nonlinear non-quadratic

int exprPow::Linearity () {

  if (arglist_ [0] -> Type () == CONST) {

    if (arglist_ [1] -> Type () == CONST) return CONSTANT;
    else                                  return NONLINEAR;
  }
  else {

    double exponent = arglist_ [1] -> Value ();

    if (fabs (exponent - COUENNE_round (exponent)) > COUENNE_EPS)
      return NONLINEAR;

    if (arglist_ [1] -> Type () == CONST) { 

      int expInt = (int) COUENNE_round (exponent);

      if (arglist_ [0] -> Linearity () == LINEAR) {

	switch (expInt) {

	case 0:  return CONSTANT;
	case 1:  return LINEAR;
	case 2:  return QUADRATIC;

	default: return NONLINEAR;
	}
      }
      else 
	if (arglist_ [0] -> Linearity () == QUADRATIC) 
	  switch (expInt) {

	  case 0:  return CONSTANT;
	  case 1:  return QUADRATIC;

	  default: return NONLINEAR;
	  }
	else return NONLINEAR;
    }
    else return NONLINEAR;
  }
}


/// is this expression integer?
bool exprPow::isInteger () {

  // base

  if (!(arglist_ [0] -> isInteger ())) { 

    // base not integer: check if constant and integer
    CouNumber lb, ub;
    arglist_ [0] -> getBounds (lb, ub);

    if ((fabs (lb - ub) > COUENNE_EPS) ||
	!::isInteger (lb))
      return false;
  }

  // exponent

  if (!(arglist_ [1] -> isInteger ())) { 

    // exponent not defined integer: check if constant and at integer
    // value (and positive, or base negative integer)

    CouNumber lb, ub;
    arglist_ [1] -> getBounds (lb, ub);

    if ((fabs (lb - ub) > COUENNE_EPS) ||
	!::isInteger (lb))
      return false;

    if (lb < 0) { // exponent negative, check again base

      arglist_ [0] -> getBounds (lb, ub);

      if ((fabs (lb - ub) > COUENNE_EPS) ||
	  (fabs (lb) < COUENNE_EPS) ||
	  !::isInteger (1. / lb))
	return false;
    }
  } else {

    // if base integer and exponent integer, must check that exponent
    // is nonnegative

    CouNumber lb, ub;
    arglist_ [1] -> getBounds (lb, ub);

    if (lb < .5)
      return false;
  }

  return true;
}


/// compute $y^{lv}$ and  $y^{uv}$ for Violation Transfer algorithm
void exprPow::closestFeasible (expression *varind,
			       expression *vardep, 
			       CouNumber &left,
			       CouNumber &right) const {
  CouNumber
    x  = (*varind) (),
    y  = (*vardep) (),
    k  = arglist_ [1] -> Value (),
    xk = safe_pow (x, k),
    yk = safe_pow (y, 1./k);

  int intk = 0;

  bool isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  // three cases: 
  // 1) k or  1/k odd,        => have either left or right
  // 2) k or  1/k even,       => may have both
  // 3) k and 1/k fractional  => have either left or right

  if (isInt || isInvInt)

    if (intk % 2) // case 1

      if (k > 0) 
	((y < xk) ? left : right) = yk; // easy, x^k is continuous

      else

	if      (y < 0.)          // third, fourth orthant
	  if (y < xk) right = yk; // in convex region y < 1/x within third orthant
	  else        left  = yk; // remaining non-convex area

	else                      // first, second orthant
	  if (y > xk) left  = yk; // in convex region y > 1/x within first orthant
	  else        right = yk; // remaining non-convex area

    else // case 2

      if (y <= 0.) // third, fourth orthant => no solution
	left = - (right = COIN_DBL_MAX);

      else

	if (k > 0) 

	  if (k < 1) // roots, have x >= 0

	    if (x > yk) left  = yk;
	    else        right = yk;

	  else

	    if (x > yk)       left  =  yk;
	    else if (x < -yk) right = -yk;
	    else              left  = - (right = yk);

	else // k negative
	  if (y < xk) // between asymptotes
	    left = - (right = yk);
	  else  // in one of the two convex areas
	    if (x > 0) left  =  yk;
	    else       right = -yk;

  else // case 3: assume x bounded from below by 0

    if (k > 0) ((y < xk) ? left : right) = yk;
    else       ((y > xk) ? left : right) = yk;
}


/// return l-2 norm of gradient at given point
CouNumber exprPow::gradientNorm (const double *x) {

  int ind0 = arglist_ [0] -> Index ();
  CouNumber exponent = arglist_ [1] -> Value ();
  return (ind0 < 0) ? 0. : fabs (exponent * safe_pow (x [ind0], exponent - 1));
}


/// can this expression be further linearized or are we on its
/// concave ("bad") side
bool exprPow::isCuttable (CouenneProblem *problem, int index) const {

  CouNumber exponent = arglist_ [1] -> Value ();

  bool
    isInt    = ::isInteger (exponent),
    isInvInt = (exponent != 0.) && ::isInteger (1. / exponent);

  int intExp = (isInt ? COUENNE_round (exponent) : (isInvInt ? COUENNE_round (1. / exponent) : 0));

  if (exponent > 0.) {

    if (isInt || isInvInt) {

      if (intExp % 2) return false; // exponent odd or 1/odd

      CouNumber 
	x = problem -> X (arglist_ [0] -> Index ()),
	y = problem -> X (index);

      if (isInt) return (y <= safe_pow (x, exponent)); // below convex curve ==> cuttable

      return (y >= safe_pow (x, exponent)); // above concave k-th root curve ==> cuttable
    } else {

      // non-integer exponent
      CouNumber 
	x = problem -> X (arglist_ [0] -> Index ()),
	y = problem -> X (index);

      return (((exponent <= 1.) && (y >= safe_pow (x, exponent))) ||
	      ((exponent >= 1.) && (y <= safe_pow (x, exponent))));
    }      
  } else {

    // non-integer exponent
    CouNumber 
      x  = problem -> X (arglist_ [0] -> Index ()),
      y  = problem -> X (index),
      lb = problem -> Lb (index),
      ub = problem -> Ub (index);

    if (isInt || isInvInt)

      if (!(intExp % 2)) return (((lb > 0) || (ub < 0)) && (y * safe_pow (fabs (x), -exponent) <= 1.));
      else               return (((lb > 0) || (ub < 0)) && (y * safe_pow (x,        -exponent) <= 1.));
    else                 return                            (y * safe_pow (x,        -exponent) <= 1.);
  }
}
