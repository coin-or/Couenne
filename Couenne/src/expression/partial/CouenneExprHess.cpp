/* $Id$
 *
 * Name:    CouenneExprHess.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Hessian of the Lagrangian, implementation
 * 
 * This file is licensed under the Eclipse Public License (EPL) (EPL)
 */

#include <stdio.h> // ! must go

#include "CoinHelperFunctions.hpp"

#include "CouenneExprHess.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

//#define DEBUG

ExprHess::ExprHess ():

  nnz_  (0),
  iRow_ (NULL),
  jCol_ (NULL),
  numL_ (NULL),
  lamI_ (NULL),
  expr_ (NULL) {}


ExprHess::~ExprHess () {

  if (nnz_) {

    free (iRow_);
    free (jCol_);

    for (int i=0; i<nnz_; i++) {

      for (int j=0; j<numL_ [i]; j++)
	free (expr_ [i] [j]);

      free (expr_ [i]);
      free (lamI_ [i]);
    }

    free (numL_);
    free (lamI_);
    free (expr_);
  }
}

/// code for refilling jacobian
void HessElemFill (int i, int level,
		   std::set <int> &list, 
		   expression *expr,
		   int *row, int **lam, expression ***eee);

/// code for refilling arrays through realloc
static void reAlloc (int nCur, int &nMax, int *&r, int *&c, int *&n, int **&l, expression ***&e);


/// code for refilling jacobian
ExprHess::ExprHess (CouenneProblem *p):

  nnz_  (0),
  iRow_ (NULL),
  jCol_ (NULL),
  numL_ (NULL),
  lamI_ (NULL),
  expr_ (NULL) {

#ifdef DEBUG
  printf ("creating Hessian\n");
#endif

  /// for each j in (obj,con)
  ///   create j->deplist()

  std::set <int> *deplist = new std::set <int> [1 + p -> nVars () + p -> nCons ()];

  int nLevels = 0; // depth of this 3D structure (<= #con + #var + 1)

  // fill in dependence list for the objective...

  p -> Obj (0) -> Body () -> DepList (deplist [nLevels++], STOP_AT_AUX);

  // ... constraints...

  for (int i = 0; i < p -> nCons (); i++) {

    CouenneConstraint *c = p -> Con (i);

    if (c -> Body () -> Type () == AUX ||
	c -> Body () -> Type () == VAR ||
	c -> Body () -> Type () == CONST) 
      continue;

    c -> Body () -> DepList (deplist [nLevels++], STOP_AT_AUX);
  }

  // ... and auxiliaries

  for (int i = 0; i < p -> nVars (); i++) {

    exprVar *e = p -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    e -> Image () -> DepList (deplist [nLevels], STOP_AT_AUX);

    // insert auxiliary variable defined here as well
    // No, don't. No second order information will be there
    //deplist [nLevels].insert (e -> Index ());

    nLevels++;
  }

#ifdef DEBUG
  for (int i=0; i<nLevels; i++) {
    printf ("level %d: ", i);
    for (std::set <int>::iterator j = deplist [i]. begin () ; j != deplist [i].end (); ++j)
      printf ("%d ", *j);
    printf ("\n");
  }
#endif

  /// for each variable i
  ///   create dense row
  ///   for each j in (obj,con)
  ///     for k in j->deplist(): k<=i
  ///       if d^2(j)/(di dj) nonzero
  ///         add term to list [i,j]
  ///   sparsify row

  int 
    nVars   = p -> nVars (),
    curSize = 0; // used to expand array through realloc

  /// for each variable (filling a row of the hessian)
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

    level++;

    for (int j = 0; j < p -> nCons (); j++) {

      CouenneConstraint *c = p -> Con (j);
      if (c -> Body () -> Type () == AUX || 
	  c -> Body () -> Type () == VAR) continue;

      HessElemFill (i, level, deplist [level], c -> Body (), row, lam, eee);

      level++;
    }

    // and auxiliaries

    for (int j = 0; j < p -> nVars (); j++) {

      exprVar *e = p -> Var (j);

      if ((e -> Type         () != AUX) || 
	  (e -> Multiplicity () <= 0)) continue;

      HessElemFill (i, level, deplist [level], e -> Image (), row, lam, eee);
      //HessElemFill (i, level, deplist [level], e,             row, lam, eee); 
      // no need to pass e itself, as its second derivative is zero

      level++;
    }

    // sparsify row, eee, lam

    for (int j=0; j <= i; j++) 
      if (row [j]) {

	reAlloc (nnz_, curSize, iRow_, jCol_, numL_, lamI_, expr_);

	iRow_ [nnz_] = i;
	jCol_ [nnz_] = j;
	numL_ [nnz_] = row [j];
	lamI_ [nnz_] = (int         *) realloc (lam [j], row [j] * sizeof (int));
	expr_ [nnz_] = (expression **) realloc (eee [j], row [j] * sizeof (expression *));

	++nnz_;
      }

    /// upper triangular matrix has to be discarded
    for (int j=i+1; j < nVars; j++) 

      if (row [j]) { // lam and eee are non nil too

	free (lam [j]);

	for (int k=0; k<row[j]; k++)
	  if (eee [j] [k])
	    delete eee [j] [k];

	free (eee [j]);
      }

    free (row);
    free (lam);
    free (eee);
  }

#ifdef DEBUG
  printf ("hessian: %d nonzeros\n", nnz_);

  for (int i=0; i<nnz_; i++) {

    printf ("[%d,%d]: %d lambdas: ", 
	    iRow_ [i], jCol_ [i], numL_ [i]); 

    fflush (stdout);

    for (int j=0; j<numL_ [i]; j++) {
      printf ("(%d,", lamI_ [i][j]); 
      fflush (stdout);
      expr_ [i][j] -> print (); 
      fflush (stdout);
      printf (") ");
    }
    printf ("\n");
  }
#endif
}


#define reallocStep 100

void HessElemFill (int i, int level,
		   std::set <int> &list, 
		   expression *expr,
		   int *row, 
		   int **lam, 
		   expression ***eee) {

  if (list. find (i) == list.end ()) 
    return;

  for (int k=0; k <= i; k++) 
    if ((k==i) || (list. find (k) != list.end ())) {

      /// common case: expr is linear, just skip

      if (expr -> Linearity () <= LINEAR)
	continue;

      // objective depends on k and i. Is its second derivative, w.r.t. k and i, nonzero?

      expression 
	*d1  = expr -> differentiate (k),
	*sd1 = d1 -> simplify (),
	*rd1 = (sd1 ? sd1 : d1),
	*d2  = rd1 -> differentiate (i),
	*sd2 = d2 -> simplify (),
	*rd2 = (sd2 ? sd2 : d2);


#ifdef DEBUG
      printf (" rd2 [x_%d, x_%d]: ", k, i); fflush (stdout); 
      rd2 -> print (); printf (" -> "); 
      rd2 -> print (); printf ("\n");
#endif
      delete d1;
      if (sd1) delete sd1;
      if (sd2) delete d2;

      if ((rd2 -> Type  () != CONST) ||
	  (rd2 -> Value () != 0.)) {

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

static void reAlloc (int nCur, int &nMax, int *&r, int *&c, int *&n, int **&l, expression ***&e) {

  if (nCur >= nMax) {

    nMax += reallocStep;

    r = (int          *) realloc (r, nMax * sizeof (int));
    c = (int          *) realloc (c, nMax * sizeof (int));
    n = (int          *) realloc (n, nMax * sizeof (int));
    l = (int         **) realloc (l, nMax * sizeof (int *));
    e = (expression ***) realloc (e, nMax * sizeof (expression **));
  }
}

