/* $Id$
 *
 * Name:    CouenneExprHess.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Hessian of the Lagrangian, implementation
 * 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneExprHess.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

ExprHess::ExprHess ():
  nnz_  (0),
  iRow_ (NULL),
  jCol_ (NULL),
  expr_ (NULL) {}


ExprHess::~ExprHess () {

  if (nnz_) {

    free (iRow_);
    free (jCol_);

    for (int i=0; i<nnz_; i++) {
      for (int j=0; j<numL_ [i]; j++)
	free (expr_ [i] [j]);
      free (expr_ [i]);
      free (lamI_);
    }

    free (numL_);
    free (expr_);
    free (lamI_);
  }
}

/// code for refilling jacobian

void HessElemFill (int i, int level,
		   std::set <int> &list, 
		   expression *expr,
		   int *row, int **lam, expression ***eee);

/// code for refilling jacobian
ExprHess::ExprHess (CouenneProblem *p) {

  /// for each j in (obj,con)
  ///   create j->deplist()

  std::set <int> *deplist = new std::set <int> [1 + p -> nVars () + p -> nCons ()];

  int nLevels = 0; // depth (<= #con + #var + 1)

  // fill in deplists for the objective...

  p -> Obj (0) -> Body () -> DepList (deplist [nLevels++], STOP_AT_AUX);

  // ... constraints...

  for (int i = 0; i < p -> nCons (); i++) {

    CouenneConstraint *c = p -> Con (i);

    if (c -> Body () -> Type () == AUX ||
	c -> Body () -> Type () == VAR) 
      continue;

    c -> Body () -> DepList (deplist [nLevels++], STOP_AT_AUX);
  }

  // and auxiliaries

  for (int i = 0; i < p -> nVars (); i++) {

    exprVar *e = p -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    e -> DepList (deplist [nLevels], STOP_AT_AUX);

    deplist [nLevels++].insert (e -> Index ());
  }

  /// for each variable i
  ///   create dense row
  ///   for each j in (obj,con)
  ///     for k in j->deplist(): k<=i
  ///       if d^2(j)/(di dj) nonzero
  ///         add term to list [i,j]
  ///   sparsify row

  int nVars = p -> nVars ();

  /// for each variable
  for (int i=0; i < nVars; i++) {

    // create dense row. These will become numL later
    int          *row = (int          *) malloc (nVars * sizeof (int));
    int         **lam = (int         **) malloc (nVars * sizeof (int *));
    expression ***eee = (expression ***) malloc (nVars * sizeof (expression **));

    CoinFillN (row, nVars, 0);
    for (int j=0; j<nVars; j++) {
      lam [j] = NULL;
      eee [j] = NULL;
    }

    // scan all levels

    int level = 0;

    /// fill term for objective

    HessElemFill (i, 0, deplist [0], p -> Obj (0) -> Body (), row, lam, eee);

    for (int j = 0; j < p -> nCons (); j++) {

      CouenneConstraint *c = p -> Con (j);
      if (c -> Body () -> Type () == AUX || c -> Body () -> Type () == VAR) continue;

      level++;

      HessElemFill (i, level, deplist [level], c -> Body (), row, lam, eee);
    }

    // and auxiliaries

    for (int j = 0; j < p -> nVars (); j++) {

      exprVar *e = p -> Var (j);

      if ((e -> Type () != AUX) || (e -> Multiplicity () <= 0))	continue;

      level++;

      HessElemFill (i, level, deplist [level], e -> Image (), row, lam, eee);
    }

    // sparsify row, eee, lam
    
    for (int j=0; j <= i; j++) 
      if (row [j]) {

	iRow_ [nnz_] = i;
	jCol_ [nnz_] = j;
	numL_ [nnz_] = row [j];
	lamI_ [nnz_] = (int         *) realloc (lam [j], row [j] * sizeof (int));
	expr_ [nnz_] = (expression **) realloc (eee [j], row [j] * sizeof (expression *));

	++nnz_;
      }
  }
}

#define reallocStep 100

void HessElemFill (int i, int level,
		   std::set <int> &list, 
		   expression *expr,
		   int *row, int **lam, expression ***eee) {

  if (list. find (i) == list.end ()) 
    return;

  for (int k=0; k <= i; k++) 
    if ((k==i) || (list. find (k) != list.end ())) {

      // objective depends on k and i. Is its second derivative, w.r.t. k and i, nonzero?

      expression 
	*d1  = expr -> differentiate (k),
	*sd1 = d1 -> simplify (),
	*rd1 = (sd1 ? sd1 : d1),
	*d2  = rd1 -> differentiate (i),
	*sd2 = d2 -> simplify (),
	*rd2 = (sd2 ? sd2 : d2);

      delete d1;
      if (sd1) delete sd1;
      if (sd2) delete d2;

      if ((rd2 -> Type () != CONST) ||
	  (rd2 -> Type () != 0.)) {

	// we have a nonzero

	if (!(row [k] % reallocStep)) {

	  lam [k] = (int         *) realloc (lam[k], (row[k] + reallocStep) * sizeof (int));
	  eee [k] = (expression **) realloc (eee[k], (row[k] + reallocStep) * sizeof (expression *));
	}

	lam [k] [row [k]] = level;
	eee [k] [row [k]] = rd2;
	row [k] ++;

      } else delete rd2;
    }
}
