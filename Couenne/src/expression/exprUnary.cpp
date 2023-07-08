/* $Id$
 *
 * Name:    exprUnary.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the unary expression class
 *
 * (C) Carnegie-Mellon University, 2006-2009.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneTypes.hpp"
#include "CouenneExprUnary.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"

using namespace Couenne;

// print unary expression
void exprUnary::print (std::ostream &out, 
		       bool descend) const {

  if (printPos () == PRE)  out << printOp ();
  out << "("; 
  argument_ -> print (out, descend); 
  out << ")";
  if (printPos () == POST) out << printOp ();
}


/// comparison when looking for duplicates
int exprUnary::compare (exprUnary  &e1) { 

  int c0 = code (),
      c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;
  else // have to compare arguments 
    return argument_ -> compare (*(e1.argument_));
}


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)
exprAux *exprUnary::standardize (CouenneProblem *p, bool addAux) {

  exprAux *subst;

  if ((subst = argument_ -> standardize (p))) {

    if ((subst -> Type () == AUX) ||
	(subst -> Type () == VAR)) 
      argument_ = new exprClone (subst);
    else argument_ = subst;
  }

  return (addAux ? (p -> addAuxiliary (this)) : new exprAux (this, p -> domain ()));
}


/// replace variable with other
void exprUnary::replace (exprVar *x, exprVar *w) {

  if (argument_ -> Type () == VAR) {
    if (argument_ -> Index () == x -> Index ()) {
      delete argument_;
      argument_ = new exprClone (w);
    }
  } else argument_ -> replace (x, w);
}


/// is this expression integer?
bool exprUnary::isInteger () {

  // only check if argument is, *at this point in the algorithm*,
  // constant -- due to branching rules, for instance. If so, check if
  // the corresponding evaluated expression is integer.

  CouNumber al, au;
  argument_ -> getBounds (al, au);

  if ((al > -COUENNE_INFINITY) &&     // Funny: if al=-(au=-1.7e308) returns true...
      (au <  COUENNE_INFINITY) &&
      fabs (al - au) < COUENNE_EPS) { // argument is constant

    CouNumber fval = (F ()) (al); 

    // check if f(lb=ub) is integer
    if (fabs (COUENNE_round (fval) - fval) < COUENNE_EPS)
      return true;
  }

  return false;
}


// simplify unary operators
expression *exprUnary:: simplify () {

  expression *subst;

  // Simplify argument g(x) of this expression f(g(x))
  if ((subst = argument_ -> simplify ())) {

    delete argument_;
    argument_ = subst;

    // g(x) is a constant k, therefore return f (k)
    if (subst -> Type () == CONST) {

      expression *ret = new exprConst (operator () ());
      argument_ = NULL;
      delete subst;

      return ret;
    } 
  }

  return NULL;
}
