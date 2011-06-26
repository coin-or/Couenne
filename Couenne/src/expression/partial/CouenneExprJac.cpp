/* $Id$
 *
 * Name:    CouenneExprMatr.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Implementation, matrix expressions
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h> // ! must go

#include "CoinHelperFunctions.hpp"

#include "CouenneExprJac.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

//#define DEBUG

// empty constructor
ExprJac::ExprJac ():
  nnz_   (0),
  iRow_  (NULL),
  jCol_  (NULL),
  expr_  (NULL),
  nRows_ (0) {}


//destructor
ExprJac::~ExprJac () {

  if (nnz_) {

    free (iRow_);
    free (jCol_);

    for (int i=0; i<nnz_; i++)
      delete expr_ [i];

    free (expr_);
  }
}

// copy constructor
ExprJac::ExprJac  (const ExprJac &rhs)
{operator= (rhs);}


// clone
ExprJac *ExprJac::clone ()
{return new ExprJac (*this);}


// assignment
ExprJac &ExprJac::operator= (const ExprJac &rhs) {

  nnz_   = rhs. nnz_;
  nRows_ = rhs. nRows_;

  iRow_ = (nnz_ && rhs.iRow_) ? (int *) malloc (nnz_ * sizeof (int)) : NULL;
  jCol_ = (nnz_ && rhs.jCol_) ? (int *) malloc (nnz_ * sizeof (int)) : NULL;

  CoinCopyN (rhs.iRow_, nnz_, iRow_);
  CoinCopyN (rhs.jCol_, nnz_, jCol_);

  if (nnz_) {

    expr_ = (expression **) malloc (nnz_ * sizeof (expression *));

    for (int i=0; i<nnz_; i++)
      expr_ [i] = expr_ [i] -> clone ();

  } else expr_ = NULL;

  return *this;
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

// constructor
ExprJac::ExprJac (CouenneProblem *p):

  nnz_   (0),
  iRow_  (NULL),
  jCol_  (NULL),
  expr_  (NULL),
  nRows_ (0) {

  /// constraints: 
  /// 
  /// If they are variable constraints, i.e., of the form a <= x <= b,
  /// they should be ignored and be replaced by variable bound (simply
  /// ask the problem).
  ///
  /// All other constraints should be part of the jacobian

  /// to be resized on demand
  int 
    cursize   = 0,
    nRealCons = 0;

  reAlloc (nnz_, cursize, iRow_, jCol_, expr_);

  // constraints ////////////////////////////////////////////////////////////

  for (int i = 0; i < p -> nCons (); i++) {

    CouenneConstraint *c = p -> Con (i);

    if (c -> Body () -> Type () == AUX ||
	c -> Body () -> Type () == VAR) 
      continue;

    // This is a constraint of the form f(x) <= 0. Find out the
    // variables (original or aux) it depends on, directly
    // (STOP_AT_AUX)

    std::set <int> deplist;

    c -> Body () -> DepList (deplist, STOP_AT_AUX);

    int nTerms = 0;

    for (std::set <int>::iterator k = deplist.begin (); k != deplist.end (); ++k) {

      expression 
	*J = c -> Body () -> differentiate (*k), // derivative of the
                  	                         // constraint's body
	                                         // w.r.t. x_i

	*sJ = J -> simplify (),                  // a simplification
	*rJ = sJ ? sJ : J;                       // the real one

      if (sJ) 
	delete J; // the only remaining expression won't be wasted

      if ((rJ -> Type  () == CONST) &&
	  (rJ -> Value () == 0.)) 
	continue;

      // there is a nonzero entry!

      reAlloc (nnz_ + 1, cursize, iRow_, jCol_, expr_);

      rJ -> realign (p);

      iRow_ [nnz_] = nRealCons;
      jCol_ [nnz_] = *k;
      expr_ [nnz_] = rJ;

      nnz_++;
      nTerms++;
    }

    if (nTerms) {
      ++nRealCons; // increase the counter of real constraints 
      ++nRows_;    // and of rows
    }
  }

  // auxiliaries ////////////////////////////////////////////////////////////

  /// Each should be considered a constraint

  for (int i = 0; i < p -> nVars (); i++) {

    exprVar *e = p -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    // this is a variable definition of the form y </>/= f(x). Find
    // out the variables (original or aux) it depends on, directly
    // (STOP_AT_AUX)

    std::set <int> deplist;

    e -> Image () -> DepList (deplist, STOP_AT_AUX);

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

      rJ -> realign (p);

      // there is a nonzero entry!

      reAlloc (nnz_ + 1, cursize, iRow_, jCol_, expr_);

      iRow_ [nnz_] = nRealCons;
      jCol_ [nnz_] = *k;
      expr_ [nnz_] = rJ;

      ++nnz_;
      ++nTerms;
    }

    if (nTerms) {
      ++nRealCons; // increase the counter of real constraints 
      ++nRows_;
    }
  }

#ifdef DEBUG
  printf ("jacobian: %d nonzeros, %d rows\n", nnz_, nRows_);

  for (int i=0; i<nnz_; i++) {

    printf ("[%d,%d]: ", iRow_ [i], jCol_ [i]);

    fflush (stdout);
    expr_ [i] -> print (); 
    printf ("\n");
  }
#endif
}
