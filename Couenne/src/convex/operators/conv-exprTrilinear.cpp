/* $Id$
 *
 * Name:    conv-exprTrilinear.cpp
 * Author:  Pietro Belotti
 * Purpose: convexify trilinear terms
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"
#include "CouenneExprTrilinear.hpp"
#include "CouenneExprMin.hpp"
#include "CouenneExprMax.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

/// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {

  register int t1 = v1 -> Type (), t2;
  return (((t1 == VAR) || (t1 == AUX)) &&
	  (((t2 = v2 -> Type ()) == VAR) || (t2 == AUX)) && 
	  (v1 -> Index () == v2 -> Index ()));
}

/// get lower/upper bounds of product f(x) g(x) in expression form
void exprTrilinear::getBounds (expression *&lb, expression *&ub) {

  expression
    **arglistMax = new expression* [16],
    **arglistMin = new expression* [16],
    **lbA        = new expression* [3],
    **ubA        = new expression* [3];

  for (int i=0; i<3; i++)
    arglist_ [i] -> getBounds (lbA [i], ubA [i]);

  for     (int i0 = 0; i0 < 2; i0++)
    for   (int i1 = 0; i1 < 2; i1++)
      for (int i2 = 0; i2 < 2; i2++) {

	int indexTerm = i0*8 + i1*4 + i2*2;

	arglistMax [indexTerm] = new exprTrilinear (new exprClone (i0 ? ubA [0] : lbA [0]),
						    new exprClone (i1 ? ubA [1] : lbA [1]),
						    new exprClone (i2 ? ubA [2] : lbA [2]));

	arglistMin [indexTerm] = new exprClone (arglistMax [indexTerm]);

	arglistMax [indexTerm + 1] = new exprStore (arglistMax [indexTerm]); // evaluated at the end, safe to just copy
	arglistMin [indexTerm + 1] = new exprStore (arglistMax [indexTerm]); // evaluated at the end, safe to just copy	
      }

  lb = new exprMin (arglistMin, 16);
  ub = new exprMax (arglistMax, 16);
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprTrilinear::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber 
    *lbA = new CouNumber [3],
    *ubA = new CouNumber [3];

  for (int i=0; i<3; i++)
    arglist_ [i] -> getBounds (lbA [i], ubA [i]);

  lb =  COUENNE_INFINITY;
  ub = -COUENNE_INFINITY;

  for     (int i0 = 0; i0 < 2; i0++)
    for   (int i1 = 0; i1 < 2; i1++)
      for (int i2 = 0; i2 < 2; i2++) {

	double curbound = 
	  (i0 ? ubA [0] : lbA [0]) * 
	  (i1 ? ubA [1] : lbA [1]) *
	  (i2 ? ubA [2] : lbA [2]);

	if (lb > curbound) lb = curbound;
	if (ub < curbound) ub = curbound;
      }
}
