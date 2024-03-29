/*
 *
 * Name:    exprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the absulute value of a function
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprAbs.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMin.hpp"
#include "CouenneExprMax.hpp"
#include "CouenneExprOpp.hpp"

#include "CouenneProblem.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

/// find lower and upper bound of a given expression

void exprAbs::getBounds (expression *&lb, expression *&ub) {

  // get bounds of the argument
  expression *lba, *uba;

  argument_ -> getBounds (lba, uba);

  // lower bound = max (0, lb, -ub)

  expression **all = new expression * [6];
  all [0] = new exprConst (0.);  all [1] = new exprConst (0.);
  all [2] = new exprOpp (uba);   all [3] = new exprOpp (new exprClone (uba));
  all [4] = lba;                 all [5] = new exprClone (lba);

  lb = new exprMax (all, 6);

  // upper bound = max (|lb|, |ub|)

  ub = new exprMax (new exprAbs (new exprClone (lba)),
		    new exprAbs (new exprClone (uba)));
}


/// find value of lower and upper bound of a given expression

void exprAbs::getBounds (CouNumber &lb, CouNumber &ub) {

  // get bounds of the argument
  CouNumber lba, uba;

  argument_ -> getBounds (lba, uba);

  if (lba > 0) {
    lb = lba;
    ub = uba;
  } else if (uba < 0) {
    lb = -uba;
    ub = -lba;
  } else {
    lb = 0.;
    ub = CoinMax (-lba, uba);
  }
}


/// differentiation

expression *exprAbs::differentiate (int index) {

  expression **arglist = new expression * [4];
  expression  *diffarg = argument_ -> differentiate (index);

  arglist [0] = new exprConst (0.);
  arglist [1] = new exprClone (diffarg);
  arglist [2] = new exprOpp (new exprClone (argument_));
  arglist [3] = new exprOpp (diffarg);

  return new exprMin (arglist, 4);
}


/// implied bound processing for expression w = |x|, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprAbs::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  int index = argument_ -> Index ();

  CouNumber *xl = l + index, wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
            *xu = u + index, wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  // for w >= b > 0, we can only improve xlb if it is at least  b
  //                                     xub             most  -b

  bool tighter = false;

  if (wl > 0) {
    if      (*xl > 0) {
      if (updateBound (-1, xl, argument_ -> isInteger () ? ceil   (wl - COUENNE_EPS) :  wl)) {
	tighter = true;
	chg [index].setLower(t_chg_bounds::CHANGED);
      }
    }
    else if (*xu < 0) {
      if (updateBound (+1, xu, argument_ -> isInteger () ? floor (-wl + COUENNE_EPS) : -wl)) {
	tighter = true;
	chg [index].setUpper(t_chg_bounds::CHANGED);
      }
    }
  }

  // w <= u (if u < 0 the problem is infeasible)

  if (wu < COUENNE_INFINITY) {
    if (updateBound (-1, xl, argument_ -> isInteger () ? ceil (-wu - COUENNE_EPS) : -wu)) {
      tighter = true;
      chg [index].setLower(t_chg_bounds::CHANGED);
    }
    if (updateBound (+1, xu, argument_ -> isInteger () ? floor (wu + COUENNE_EPS) :  wu)) {
      tighter = true;
      chg [index].setUpper(t_chg_bounds::CHANGED);
    }
  }

  return tighter;
}


/// closest feasible points in function in both directions
void exprAbs::closestFeasible (expression *varind, expression *vardep,
			       CouNumber& left, CouNumber& right) const
{
  CouNumber valdep = (*vardep)();
  CouNumber curr = (*varind)();

  if (valdep < 0.) { // no way to restore feasibility, set infinite interval
    left = -COUENNE_INFINITY;
    right = COUENNE_INFINITY;
  }
  else if (curr < -valdep) { // hence curr negative and behind left part of |x|
    left = curr;
    right = -valdep;
  }
  else if (curr > valdep) { // hence curr positive and after right half-line of |x|
    left = valdep;
    right = curr;
  }
  else { // between the two half-lines
    left = -valdep;
    right = valdep;
  }

}


/// can this expression be further linearized or are we on its
/// concave ("bad") side
bool exprAbs::isCuttable (CouenneProblem *problem, int index) const {

  double
    x = problem -> X (argument_ -> Index ()),
    y = problem -> X (index);

  return ((y <= x) || (y <= -x));
}


/// either CONVEX, CONCAVE, AFFINE, or NONCONVEX
//enum convexity exprAbs::convexity () const {

//  CouNumber lb, ub;
//  getBounds (lb, ub);
//  return (((lb >= 0.) || (ub <= 0.)) ? AFFINE : CONVEX);
//}
