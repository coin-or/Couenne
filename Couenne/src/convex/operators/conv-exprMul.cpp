/* $Id$
 *
 * Name:    conv-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: utility methods to convexify multiplications
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <queue>

#include "CouenneTypes.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprTrilinear.hpp"
#include "CouenneExprBMul.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (expression *&lb, expression *&ub) {

  int i;

  if ((arglist_ [i=0] -> Type () == CONST) ||
      (arglist_ [i=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [i] -> Value ();

    if (!i && (arglist_ [1] -> Type () == CONST)) { 

      // !i means i==0, or the first is constant. If you are here,
      // both are constant, which should not happen...

      CouNumber prod = c * arglist_ [1] -> Value ();

      lb = new exprConst (prod);
      ub = new exprConst (prod);

      return;
    }
    else {

      // expression is of the type c*x

      expression *lbi, *ubi;
      arglist_ [1-i] -> getBounds (lbi, ubi);

      if (c >= 0) {
	lb = new exprMul (new exprConst (c), lbi);
	ub = new exprMul (new exprConst (c), ubi);
      } else {
	lb = new exprMul (new exprConst (c), ubi);
	ub = new exprMul (new exprConst (c), lbi);
      }
    }
  }
  else {

    // expression is of the type x*y

    expression **almin = new expression * [4];
    expression **almax = new expression * [4];

    arglist_ [0] -> getBounds (almin [0], almin [1]);
    arglist_ [1] -> getBounds (almin [2], almin [3]);

    almax [0] = new exprClone (almin [0]);
    almax [1] = new exprClone (almin [1]);
    almax [2] = new exprClone (almin [2]);
    almax [3] = new exprClone (almin [3]);

    lb = new exprLBMul (almin, 4);
    ub = new exprUBMul (almax, 4);
  }
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lb1, ub1, lb2, ub2;

  arglist_ [0] -> getBounds (lb1, ub1);
  arglist_ [1] -> getBounds (lb2, ub2);

  if (ub1 < 0) { // use lb1, dominant
    if      (ub2 < 0) {lb = safeProd(ub1,ub2); ub = safeProd(lb1,lb2);}
    else if (lb2 > 0) {lb = safeProd(lb1,ub2); ub = safeProd(ub1,lb2);}
    else              {lb = safeProd(lb1,ub2); ub = safeProd(lb1,lb2);}
  } else if (lb1 > 0) { // use ub1, dominant
    if      (ub2 < 0) {lb = safeProd(ub1,lb2); ub = safeProd(lb1,ub2);}
    else if (lb2 > 0) {lb = safeProd(lb1,lb2); ub = safeProd(ub1,ub2);}
    else              {lb = safeProd(ub1,lb2); ub = safeProd(ub1,ub2);}
  } else { // there is a zero to consider
    if      (ub2 < 0) {lb = safeProd(ub1,lb2); ub = safeProd(lb1,lb2);}
    else if (lb2 > 0) {lb = safeProd(lb1,ub2); ub = safeProd(ub1,ub2);}
    else              {lb = CoinMin (safeProd(lb1,ub2), safeProd(lb2,ub1));
                       ub = CoinMax (safeProd(lb1,lb2), safeProd(ub1,ub2));}
  }
}
