/*
 *
 * Name:    exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for class logarithm
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <math.h>

#include "CouenneExprLog.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMax.hpp"
#include "CouenneExprMin.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneProblem.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

/// get bounds of log (x) based on bounds of x

void exprLog::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  // [low|upp]er bound of w=log(x) is log (max (0, [low|upp]er (x)))
  //  lb = new exprLog (new exprMax (new exprConst (1e-100), lba));
  //  ub = new exprLog (new exprMax (new exprConst (1e-100), uba));

  expression **all  = new expression * [4];

  all [0] = new exprClone (lba); all [1] = new exprLog (lba);
  all [2] = new exprConst (0);   all [3] = new exprConst (- COUENNE_INFINITY);
  lb = new exprMax (all, 4);

  expression **alu  = new expression * [4],
             **alum = new expression * [4];

  alum [0] = new exprConst (COUENNE_INFINITY);
  alum [1] = new exprConst (COUENNE_INFINITY);
  alum [2] = new exprClone (uba);
  alum [3] = new exprLog (uba);

  alu [0] = new exprClone (uba); alu [1] = new exprMin (alum, 4);
  alu [2] = new exprConst (0.);  alu [3] = new exprConst (- COUENNE_INFINITY);
  ub = new exprMax (alu, 4);
}


/// get bounds of log (x) based on bounds of x

void exprLog::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lba, uba;
  argument_ -> getBounds (lba, uba);

  lb = log (CoinMax (1e-50, lba));
  ub = log (CoinMax (1e-50, uba));
}


/// differentiation
expression *exprLog::differentiate (int index) {
  return new exprDiv (argument_ -> differentiate (index),
		      new exprClone (argument_));
}


/// implied bound processing for expression w = log(x), upon change in
/// lower- and/or upper bound of w, whose index is wind
bool exprLog::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  int ind = argument_ -> Index ();

  bool
    res   = false,
    isint = argument_ -> isInteger();

  CouNumber
    wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
    wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  if (updateBound (-1, l+ind, isint ? ceil (exp (wl) - COUENNE_EPS) : exp (wl))) {
    res = true;
    chg [ind].setLower (t_chg_bounds::CHANGED);
  }

  if (updateBound ( 1, u+ind, isint? floor (exp (wu) + COUENNE_EPS) : exp (wu))) {
    res = true;
    chg [ind].setUpper (t_chg_bounds::CHANGED);
  }

  return res;
}


/// return l-2 norm of gradient at given point
CouNumber exprLog::gradientNorm (const double *x) {
  return (argument_ -> Index () < 0) ? 0. :
    1. / (CoinMax (1. / COUENNE_INFINITY, x [argument_ -> Index ()]));
}


/// can this expression be further linearized or are we on its
/// concave ("bad") side
bool exprLog::isCuttable (CouenneProblem *problem, int index) const {

  double
    x = problem -> X (argument_ -> Index ()),
    y = problem -> X (index);

  return ((x == 0.) || (y > log (x)));
}
