/*
 *
 * Name:    exprCos.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for cosines
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <math.h>

#include "CouenneExprCos.hpp"
#include "CouenneExprSin.hpp"
#include "CouenneExprBCos.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprClone.hpp"

using namespace Couenne;

static const CouNumber
  pi  = M_PI,
  pi2 = M_PI * 2.,
  pih = M_PI / 2.;


// return an expression -sin (argument), the derivative of cos (argument)
expression *exprCos::differentiate (int index) {

  return new exprOpp (new exprMul (new exprSin (new exprClone (argument_)),
				   argument_ -> differentiate (index)));
}


// compute bounds of sin x given bounds of x
void exprCos::getBounds (expression *&lb, expression *&ub) {

  expression *xl, *xu;
  argument_ -> getBounds (xl, xu);

  lb = new exprLBCos (xl, xu);
  ub = new exprUBCos (new exprClone (xl), new exprClone (xu));
}

// compute value of bounds of cos x given bounds of x
void exprCos::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber l, u;
  argument_ -> getBounds (l, u);

  if ((u - l >= pi2) ||      // 1) interval spans whole cycle
      (floor (l/pi2 - 0.5) < // 2) there is a pi + 2k pi between l and u
       floor (u/pi2 - 0.5)))
    lb = -1.;
  else lb = CoinMin (cos (l), cos (u));

  if ((u - l >= pi2) || // 1) interval spans whole cycle
      (floor (l/pi2) <  // 2) there is a 3/2 pi + 2k pi between l and u
       floor (u/pi2)))
    ub = 1.;
  else ub = CoinMax (cos (l), cos (u));
}


/// closest feasible points in function in both directions
void exprCos::closestFeasible (expression *varind, expression *vardep,
			       CouNumber& left, CouNumber& right) const
{
  CouNumber curr = (*varind)();
  int period = (int)(curr/pi2);
  CouNumber curr_noperiod = curr - pi2*period;
  CouNumber inv = acos((*vardep)());

  if (curr_noperiod < inv) {
    left = pi2*period - inv;
    right = pi2*period + inv;
  }
  else if (curr_noperiod < pi2-inv) {
    left = pi2*period + inv;
    right = pi2*(period+1) - inv;
  }
  else {
    left = pi2*(period+1) - inv;
    right = pi2*(period+1) + inv;
  }
}
