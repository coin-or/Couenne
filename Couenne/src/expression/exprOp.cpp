/* $Id$
 *
 * Name:    exprOp.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for multivariate function class
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExpression.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprOp.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprConst.hpp"
#include "CouennePrecisions.hpp"

#include <cstdio>

using namespace Couenne;

namespace Couenne {
class CouenneProblem;
class Domain;
}

// General N-ary function destructor
exprOp::~exprOp () {

  if (arglist_) {
    for (expression **alist = arglist_; nargs_--; alist++)
      if (*alist) delete (*alist);
    delete [] arglist_;
  }
}


// print expression

void exprOp::print (std::ostream &out, 
		    bool descend) const {
  
  if (printPos () == PRE)
    out << printOp ();

  if (nargs_ > 1)
    {out << "("; fflush (stdout);}
  for (int i=0; i<nargs_; i++) {
    if (arglist_ [i])
      arglist_ [i] -> print (out, descend); 
    fflush (stdout);
    if (i < nargs_ - 1) {
      if (printPos () == INSIDE) out << printOp ();
      else                       out << ",";
    }
    if (!((i + 1) % MAX_ARG_LINE))
      out << std::endl;
    fflush (stdout);
  }
  if (nargs_ > 1) {
    out << ")";
    fflush (stdout);
  }
}


/// compare general n-ary expressions

int exprOp::compare (exprOp &e1) {

  int c0 =     code (),
      c1 = e1. code ();

  if (c0 < c1) return -1;
  if (c0 > c1) return  1;

  // have to compare arguments one by one
  if (nargs_ < e1.nargs_) return -1;
  if (nargs_ > e1.nargs_) return  1;

  // not an exprGroup, compare arguments
  for (int i = nargs_; i--;) {

    int res = arglist_ [i] -> compare (*(e1. ArgList () [i]));
    if (res) return res;
  }

  // last chance, this might be an exprGroup or derived
  if ((c0 == COU_EXPRGROUP) ||
      (c0 == COU_EXPRQUAD)) {

    exprGroup *ne0 = dynamic_cast <exprGroup *> (this),
              *ne1 = dynamic_cast <exprGroup *> (&e1);

    int cg = ne0 -> compare (*ne1);

    if (cg) return cg; // exprGroup

    // last chance, the two are quadratic forms

    if (c0 == COU_EXPRQUAD) {

      exprQuad *ne0 = dynamic_cast <exprQuad *> (this),
   	       *ne1 = dynamic_cast <exprQuad *> (&e1);

      return ne0 -> compare (*ne1);
    }
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprOp::rank () {

  int maxrank = -1;

  for (expression **al = arglist_ + nargs_; 
       al-- > arglist_;) {
    int r = (*al) -> rank ();
    if (r > maxrank) maxrank = r;
  }

  return (maxrank);
}


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)
//
// For the base exprOp class we only do the first part (for argument
// list components only), and the calling class (Sum, Sub, Mul, Pow,
// and the like) will do the part for its own object

exprAux *exprOp::standardize (CouenneProblem *p, bool addAux) {

  exprVar *subst;

  for (int i = 0; i < nargs_; ++i)
    if ((subst = arglist_ [i] -> standardize (p))) {

      if ((subst -> Type () == VAR) ||
	  (subst -> Type () == AUX))
	arglist_ [i]    = new exprClone (subst);
      else arglist_ [i] = subst; // possibly a constant, should be nothing else
    }
  return NULL;
}


/// replace variable x with new (aux) w
void exprOp::replace (exprVar *x, exprVar *w) {

  expression **al = arglist_;
  int index = x -> Index ();

  for (int i = nargs_; i--; al++)

    switch ((*al) -> Type ()) {

    case AUX:
    case VAR:
      if ((*al) -> Index () == index) {
	delete *al;
	*al = new exprClone (w);
      }
      break;

    case UNARY:
    case N_ARY:
      (*al) -> replace (x, w);
      break;

    default:
      break;
    }
}


/// is this expression integer?
bool exprOp::isInteger () {

  for (int i = nargs_; i--;)

    if (!(arglist_ [i] -> isInteger ())) { 

      // this argument is not integer: check if constant and integer

      CouNumber lb, ub;
      arglist_ [i] -> getBounds (lb, ub);

      if ((fabs (lb - ub) > COUENNE_EPS) ||
	  !::isInteger (lb))
	return false;
    }

  return true;
}


/// fill in the set with all indices of variables appearing in the
/// expression
int exprOp::DepList (std::set <int> &deplist, 
		     enum dig_type type) {
  int tot = 0;

  //printf ("  exprop::deplist of %x: ", Original ()); print (); printf ("\n");

  for (int i = nargs_; i--;) {

    /*printf ("  ");
    arglist_ [i] -> print (); printf (": {");
    for (std::set <int>::iterator j=deplist.begin(); j != deplist.end(); ++j)
      printf ("%d ", *j);
      printf ("} -> {");*/

    tot += arglist_ [i] -> DepList (deplist, type);

    /*for (std::set <int>::iterator j=deplist.begin(); j != deplist.end(); ++j)
      printf ("%d ", *j);
      printf ("}\n");*/
  }

  return tot;
}

/// empty function to redirect variables to proper variable vector
void exprOp::realign (const CouenneProblem *p) {

  for (int i=0; i<nargs_; i++)
    arglist_ [i] -> realign (p);
}


// simplify n-ary expression f (g_1(x), g_2(x)... g_n(x))
expression *exprOp:: simplify () {

  //  Simplify arguments g_1(x), g_2(x)... g_n(x) first
  for (int i=0; i<nargs_; i++) {

    expression *subst;

    if ((subst = arglist_ [i] -> simplify ())) {

      delete arglist_ [i];
      arglist_ [i] = subst;
    }
  }

  return NULL;
}


//
// shrink argument list
//
// used by + and * (for now), accepts a constant resulting from
// applying an operator to the constants in the list of (pointers to)
// function arguments contained in el. The constant is inserted in the
// list if the result is not equal to null_element or if there are
// other non-constant terms in the arglist.
// 
// Example: f(x) + 3 + g(x) + 2 + 4 
//
// el    = {pf, NULL, pg, NULL, NULL}
// nargs = 5 
// c     = 3 + 2 + 4 = 9
// null_element = 0 (for sums)
// 
// where pf and pg are pointers to expression containing f and g,
// resp.
//
// Result: el and nargs are changed to
//
// el    = {pf, p9, pg}
// nargs = 3 
//
// Another example: f(x) + 2 + g(x) + (-4) + 2
// Result:
// el    = {pf, pg}
// nargs = 2
//
// Another example: f(x) * 3 * g(x) * 2 
//
// el    = {pf, NULL, pg, NULL}
// nargs = 4
// c     = 3 * 2 = 6 != null_element = 1 (for multiplications)
// Result:
// el    = {pf, p6, pg}
// nargs = 3
//

int exprOp::shrink_arglist (CouNumber c, CouNumber null_element) {

  int i=0, j=0;

  bool one_fun = false;

  // find first NULL spot (left by some constant)
  while ((i < nargs_) && (arglist_ [i])) 
    i++; 

  // no spots, leave
  if (i==nargs_) 
    return 0;

  // check if there is at least one non-constant expression
  for (int k=nargs_; k--;) 
    if (arglist_ [k]) {
      one_fun = true;
      break;
    }

  // add constant term if c is not null w.r.t. the operation or if it
  // would be an empty operand list otherwise
  if ((fabs (c - null_element) > COUENNE_EPS) || !one_fun)
    arglist_ [i++] = new exprConst (c);

  j = i;

  // now shift back all operands to compress argument list
  while (i < nargs_) {

    while ((i < nargs_) && !(arglist_ [i])) 
      i++;

    if (i < nargs_) 
      one_fun = true;

    while ((i < nargs_) && (arglist_ [i]))
      arglist_ [j++] = arglist_ [i++]; 
  }

  nargs_ = j;

  // only say shrinking simplified arg list if there is just the
  // constant
  return (nargs_ == 1);// && ((fabs (c - null_element) > COUENNE_EPS) || !one_fun));
}
