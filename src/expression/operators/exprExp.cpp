/*
 *
 * Name:    exprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprExp.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneProblem.hpp"

#include "CoinFinite.hpp"

using namespace Couenne;

// differentiation
expression *exprExp::differentiate (int index) {

  return new exprMul (new exprExp (new exprClone (argument_)),
		      argument_ -> differentiate (index));
}


// Get expressions of lower and upper bound of an expression (if any)
void exprExp::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  lb = new exprExp (lba);
  ub = new exprExp (uba);
}


// Get value of lower and upper bound of an expression (if any)
void exprExp::getBounds (CouNumber &lb, CouNumber&ub) {

  CouNumber lba, uba;
  argument_ -> getBounds (lba, uba);

  lb = exp (lba);
  ub = exp (uba);
}


/// implied bound processing for expression w = exp(x), upon change in
/// lower- and/or upper bound of w, whose index is wind
bool exprExp::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  bool resU, resL = resU = false;
  int ind = argument_ -> Index ();

  CouNumber b;

  if ((b = sign == expression::AUX_GEQ ? 0.           : l [wind]) > 0.) // lower bound
    resL = updateBound (-1, l + ind, argument_->isInteger () ? ceil  (log (b)) : log (b));

  if ((b = sign == expression::AUX_LEQ ? COIN_DBL_MAX : u [wind]) < COIN_DBL_MAX / 10.) { // upper bound

    if ((b >= -0.) && (b < COUENNE_EPS)) // to prevent infeasibilities due to numerics
      b = COUENNE_EPS;

    resU = updateBound ( 1, u + ind, argument_ -> isInteger () ? floor (log (b)) : log (b));
  }

  if (b < - COUENNE_EPS) {
    // make it infeasible
    resU = updateBound ( 1, u + ind, -1.) || true;
    resL = updateBound (-1, l + ind,  1.) || true;
  }

  if (resL) chg [ind].setLower (t_chg_bounds::CHANGED);
  if (resU) chg [ind].setUpper (t_chg_bounds::CHANGED);

  return (resL || resU);
}


/// can this expression be further linearized or are we on its
/// concave ("bad") side
bool exprExp::isCuttable (CouenneProblem *problem, int index) const {

  double
    x = problem -> X (argument_ -> Index ()),
    y = problem -> X (index);

  return (y <= exp (x));
}
