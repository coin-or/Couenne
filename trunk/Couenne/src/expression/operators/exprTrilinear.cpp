/* $Id$
 *
 * Name:    exprTrilinear.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of trilinear terms
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdlib.h>
#include <assert.h>

#include "CoinHelperFunctions.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprTrilinear.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprClone.hpp"
#include "CouennePrecisions.hpp"

using namespace Couenne;

/// Constructors, destructor
exprTrilinear::exprTrilinear  (expression **al, int n): 
  exprMul (al, n) {}


/// constructor with only two factors
exprTrilinear::exprTrilinear (expression *arg0, expression *arg1, expression *arg2):

  exprMul (NULL, 0) {

  nargs_ = 3;
  arglist_ = new expression * [nargs_];

  arglist_ [0] = arg0;
  arglist_ [1] = arg1;
  arglist_ [2] = arg2;

  qsort (arglist_, nargs_, sizeof (expression*), compareExpr);
}


/// compute $y^{lv}$ and  $y^{uv}$ for Violation Transfer algorithm
void exprTrilinear::closestFeasible (expression *varind,
				     expression *vardep,
				     CouNumber &left,
				     CouNumber &right) const {

  printf ("using VT and trilinear terms: not implemented yet\n");
  exit (-1);

  // TODO: if use quadrilinear and VT, fill this

  // expression *varoth = arglist_ [0]; // suppose $w = cy$;

  // if (varoth -> Index () == varind -> Index ())
  //   varoth = arglist_ [1]; // actually no, it's $w = x*c$

  // assert (varoth -> Index () >= 0);

  // CouNumber
  //   x = (*varind) (),
  //   y = (*vardep) (),
  //   c = (*varoth) ();

  // if (c < 0.)
  //   if (y < c*x) {assert (y/c >= right); right = y/c;}
  //   else         {assert (y/c <= left);  left  = y/c;}
  // else if (c > 0.)
  //   if (y < c*x) {assert (y/c <= left);  left  = y/c;}
  //   else         {assert (y/c >= right); right = y/c;}
  // else left = - (right = COIN_DBL_MAX);

}


/// return l-2 norm of gradient at given point
CouNumber exprTrilinear::gradientNorm (const double *x) {

  int 
    ind0 = arglist_ [0] -> Index (),
    ind1 = arglist_ [1] -> Index (),
    ind2 = arglist_ [2] -> Index ();

  CouNumber 
    x0 = (ind0 < 0) ? arglist_ [0] -> Value () : x [ind0],
    x1 = (ind1 < 0) ? arglist_ [1] -> Value () : x [ind1],
    x2 = (ind1 < 0) ? arglist_ [2] -> Value () : x [ind2];

  if (ind0 < 0)
    if (ind1 < 0)
      if (ind2 < 0) return 0.;                            // c*d*e
      else          return fabs (x0*x1);                  // c*d*y
    else
      if (ind2 < 0) return fabs (x0*x2);                  // c*y*e
      else          return fabs (x0*sqrt(x1*x1 + x2*x2)); // c*y*z
  else 
    if (ind1 < 0)
      if (ind2 < 0) return fabs (x1*x2);                  // x*d*e
      else          return fabs (x1*sqrt(x0*x0 + x2*x2)); // x*d*z
    else
      if (ind2 < 0) return fabs (x2*sqrt(x0*x0 + x1*x1)); // x*y*e
      else          return sqrt (x0*x0 + x1*x1 + x2*x2);  // x*y*z
}
