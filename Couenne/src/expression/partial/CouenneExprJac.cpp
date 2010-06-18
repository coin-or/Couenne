/* $Id$
 *
 * Name:    CouenneExprMatr.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Implementation, matrix expressions
 * 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneExprJac.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

ExprJac::ExprJac ():
  nnz_   (0),
  iRow_  (NULL),
  jCol_  (NULL),
  expr_  (NULL),
  nRows_ (0) {}


ExprJac::~ExprJac () {

  if (nnz_) {

    free (iRow_);
    free (jCol_);

    for (int i=0; i<nnz_; i++)
      free (expr_ [i]);

    free (expr_);
  }
}

/// code for refilling jacobian

#define reallocStep 100
static void reAlloc (int nCur, int &nMax, int *&r, int *&c, expression **&e) {

  if (nCur >= nMax) {

    nMax += reallocStep;

    r = (int         *) realloc (r, nMax * sizeof (int));
    c = (int         *) realloc (c, nMax * sizeof (int));
    e = (expression **) realloc (e, nMax * sizeof (expression *));
  }
}

ExprJac::ExprJac (CouenneProblem *p) {

  /// constraints: 
  /// 
  /// If they are of the form a <= x <= b, they should be ignored and
  /// be replaced by variable bound (simply ask the problem)
  ///
  /// All other constraints should be part of the jacobian

  nRows_ = 0;

  /// to be resized on demand
  int 
    cursize   = 0,
    nRealCons = 0;

  reAlloc (nnz_ = 0, cursize, iRow_, jCol_, expr_);

  // constraints ////////////////////////////////////////////////////////////

  for (int i = 0; i < p -> nCons (); i++) {

    CouenneConstraint *c = p -> Con (i);

    if (c -> Body () -> Type () == AUX ||
	c -> Body () -> Type () == VAR) 
      continue;

    // this is a constraint of the form f(x) <= 0. Find out the
    // variables (original or aux) it depends on, directly
    // (STOP_AT_AUX)

    nRows_ ++;

    std::set <int> deplist;

    c -> Body () -> DepList (deplist, STOP_AT_AUX);

    int nTerms = 0;

    for (std::set <int>::iterator k = deplist.begin (); k != deplist.end (); ++k) {

      expression 
	*J = c -> Body () -> differentiate (*k), // derivative of the
	// constraint's body
	// w.r.t. x_i
	*sJ = J -> simplify (),                 // a simplification
	*rJ = sJ ? sJ : J;                      // the real one

      if (sJ) 
	delete J; // the only remaining expression won't be wasted

      if ((rJ -> Type  () == CONST) &&
	  (rJ -> Value () == 0.)) 
	continue;

      // there is a nonzero entry!

      if (!nTerms)
	nRealCons++; // increase the counter of real constraints 

      reAlloc (nnz_ + 1, cursize, iRow_, jCol_, expr_);

      iRow_ [nnz_] = nRealCons;
      jCol_ [nnz_] = *k;
      expr_ [nnz_] = rJ;

      nnz_++;
      nTerms++;
    }
  }

  // auxiliaries ////////////////////////////////////////////////////////////

  /// Each should be considered a constraint

  for (int i = 0; i < p -> nVars (); i++) {

    exprVar *e = p -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    // this is a variable definition of the form y = f(x). Find out
    // the variables (original or aux) it depends on, directly
    // (STOP_AT_AUX)

    nRows_ ++;

    std::set <int> deplist;

    e -> DepList (deplist, STOP_AT_AUX);

    deplist.insert (e -> Index ());

    int nTerms = 0;

    for (std::set <int>::iterator k = deplist.begin (); k != deplist.end (); ++k) {

      expression 
	*J = (*k == e -> Index ()) ? 
  	     new exprConst (-1.) :
	     e -> Image () -> differentiate (*k), // derivative of the
	                                          // constraint's body
			  		          // w.r.t. x_i
	*sJ = J -> simplify (),                   // a simplification
	*rJ = sJ ? sJ : J;                        // the real one

      if (sJ) 
	delete J; // the only remaining expression won't be wasted

      if ((rJ -> Type  () == CONST) &&
	  (rJ -> Value () == 0.)) 
	continue; 

      // there is a nonzero entry!

      if (!nTerms)
	nRealCons++; // increase the counter of real constraints 

      reAlloc (nnz_ + 1, cursize, iRow_, jCol_, expr_);

      iRow_ [nnz_] = nRealCons;
      jCol_ [nnz_] = *k;
      expr_ [nnz_] = rJ;

      nnz_++;
      nTerms++;
    }
  }
}
