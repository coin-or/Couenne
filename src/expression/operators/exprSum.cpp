/*
 *
 * Name:    exprSum.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdlib.h>

#include "CouenneExprSum.hpp"
#include "CouenneExprCopy.hpp"
#include "CouenneExprConst.hpp"

using namespace Couenne;

/// Constructor
exprSum::exprSum  (expression **al, int n):
  exprOp (al, n) { //< non-leaf expression, with argument list

  if (al==NULL) {
    arglist_ = new expression * [1];
    *arglist_ = new exprConst (0);
    nargs_ = 1;
  }

  // commutative operator, sort elements
  if (nargs_ > 1)
    qsort (arglist_, nargs_, sizeof (expression*), compareExpr);
}


/// Copy constructor
exprSum::exprSum (expression *arg0, expression *arg1):
  exprOp (arg0, arg1) {

  //  qsort (arglist_, nargs_, sizeof (expression*), compareExpr); // useless...

  if (arg0 -> compare (*arg1) > 0) { // swap elements
    expression *swap = arglist_ [0];
    arglist_ [0] = arglist_ [1];
    arglist_ [1] = swap;
  }
}


/// simplify sums

expression *exprSum::simplify () {

  exprOp:: simplify ();

  if (nargs_ == 1) {

    expression *ret = arglist_ [0];
    arglist_ [0] = NULL;
    return ret;
  }

  // from here on we assume the operands have been simplified

  CouNumber total     = 0;
  bool      found_one = false;

  for (int i=0; i<nargs_; i++) {

    // check for constant operands in multiplications

    if (arglist_ [i] -> Type () == CONST) {

      total += arglist_ [i] -> Value ();
      found_one = true;
      delete arglist_ [i];
      arglist_ [i] = NULL;
    }
  }

  /*
  if (found_one && shrink_arglist (total, 0))
    return new exprConst (arglist_ [0] -> Value ());
  else return NULL;
  */

  if (found_one && shrink_arglist (total, 0.)) {
    expression *ret = arglist_ [0];
    arglist_ [0] = NULL;
    return ret;
  }
  else return NULL;
}


/// differentiate sum of expressions

expression *exprSum:: differentiate (int index) {

  //  exprOp::differentiate (index); // FIX IT: why is it called?

  expression **arglist = new expression * [nargs_];

  int nonconst = 0;

  for (int i = 0; i < nargs_; i++)
    if (arglist_ [i] -> dependsOn (index))
      arglist [nonconst++] = arglist_ [i] -> differentiate (index);

  if (!nonconst) {
    delete [] arglist;
    return new exprConst (0.);
  }
  else return new exprSum (arglist, nonconst);
}


/// Get lower and upper bound of an expression (if any)
void exprSum::getBounds (expression *&lb, expression *&ub) {

  expression **all = new expression * [nargs_];
  expression **alu = new expression * [nargs_];

  for (int i=0; i<nargs_; i++)
    arglist_ [i] -> getBounds (all [i], alu [i]);

  lb = new exprSum (all, nargs_);
  ub = new exprSum (alu, nargs_);
}


/// Get lower and upper bound of an expression (if any)
void exprSum::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber tlb, tub;

  lb = ub = 0;

  for (int i=0; i<nargs_; i++) {
    arglist_ [i] -> getBounds (tlb, tub);
    lb += tlb;
    ub += tub;
  }
}


/// get a measure of "how linear" the expression is (see CouenneTypes.h)
int exprSum::Linearity () {

  int linmax = arglist_ [0] -> Linearity ();

  for (int i=1; i<nargs_; i++) {
    int lin = arglist_ [i] -> Linearity ();
    if (lin > linmax) linmax = lin;
  }
  return linmax;
}
