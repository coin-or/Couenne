/*
 *
 * Name:    exprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h> // ! must go

#include "CouenneExprInv.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExpression.hpp"

#include "CoinFinite.hpp"

using namespace Couenne;

// differentiation
expression *exprInv::differentiate (int index) {

  return new exprOpp (new exprDiv (argument_ -> differentiate (index),
				   new exprPow (new exprClone (argument_),
						new exprConst (2.))));
}


// printing
void exprInv::print (std::ostream &out,
		     bool descend) const {
  out << "(1/";
  argument_ -> print (out, descend);
  out << ")";
}
//  exprUnary::print (out, "1/", PRE);}


/// general function to tighten implied bounds of a function w = x^k,
/// k negative, integer or inverse integer, and odd
void invPowImplBounds (int wind, int index,
		       CouNumber *l, CouNumber *u, CouNumber k,
		       bool &resL, bool &resU,
		       enum expression::auxSign sign) {

  CouNumber wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
            wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  // 0 <= l <= w <= u

  if (wl >= 0.) {
    if (wu > COUENNE_EPS) {
      if (wu < COUENNE_INFINITY) resL = updateBound (-1, l + index, pow (wu, k));
      else                       resL = updateBound (-1, l + index, 0.);
    }
    if (wl > COUENNE_EPS)        resU = updateBound (+1, u + index, pow (wl, k));
  }

  // l <= w <= u <= 0

  if (wu <= -0.) {
    if (wl < - COUENNE_EPS) {
      if (wl > - COUENNE_INFINITY) resU = updateBound (+1, u + index, pow (wl, k)) || resU;
      else                         resU = updateBound (+1, u + index, 0.)          || resU;
    }
    if (wu < - COUENNE_EPS)        resL = updateBound (-1, l + index, pow (wu, k)) || resL;
  }
}


/// implied bound processing for expression w = 1/x, upon change in
/// lower- and/or upper bound of w, whose index is wind
bool exprInv::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  // Expression w = 1/x: we can only improve the bounds if
  //
  //    0 <= l <= w <= u         or
  //         l <= w <= u <= 0.
  //
  // Then 1/u <= x <= 1/l (given l, u finite and nonzero)

  int index = argument_ -> Index ();

  bool resL, resU = resL = false;

  invPowImplBounds (wind, index, l, u, -1., resL, resU, sign);

  bool argInt = argument_ -> isInteger ();

  if (resL) {
    chg [index].setLower(t_chg_bounds::CHANGED);
    if (argInt) l [index] = ceil  (l [index] - COUENNE_EPS);
  }

  if (resU) {
    chg [index].setUpper(t_chg_bounds::CHANGED);
    if (argInt) u [index] = floor (u [index] + COUENNE_EPS);
  }

  return (resL || resU);
}


/// return l-2 norm of gradient at given point
CouNumber exprInv::gradientNorm (const double *x) {
  int ind = argument_ -> Index ();
  CouNumber xx;
  if (ind < 0) return 0.;
  else {
    xx = x [argument_ -> Index ()];
    return 1. / (xx*xx);
  }
}


/// can this expression be further linearized or are we on its
/// concave ("bad") side
bool exprInv::isCuttable (CouenneProblem *problem, int index) const {

  int xind = argument_ -> Index ();

  double
    x = problem -> X (xind),
    y = problem -> X (index);

  return (((problem -> Lb (xind) >= 0) && (x > 0) && (y*x <= 1)) ||
	  ((problem -> Ub (xind) <= 0) && (x < 0) && (y*x <= 1)));
}
