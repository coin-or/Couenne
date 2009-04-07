/*
 * Name:    CouenneOS.cpp
 * Authors: 
 *          
 * Purpose: Creates a CouenneProblem object from an OSil instance
 *
 * (C) Carnegie-Mellon University, 2009. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"

#include "CouenneTypes.hpp"

#include "exprSum.hpp"
#include "exprMul.hpp"
#include "exprClone.hpp"
#include "exprGroup.hpp"


CouenneProblem *readProblem (Bonmin::BabSetupBase *base = NULL,
			     JnlstPtr jnlst = NULL
			     /* OSil pointer? */) {

  CouenneProblem *p = new CouenneProblem (NULL, base, jnlst);

  //p -> setProblemName (filename); // sets filename member, for later stats -- TODO

  // number of defined variables (aka common expressions)
  //ndefined_ = 0; // Don't know if OS has them: in AMPL they are
                 // defined by "var y := f(x)", with f(x) some
		 // expression

  int n_var = 0; // to be set to no. of variables

  // nonlinear in both objectives and constraints
  for (int i = 0; i < n_var; i++) 
    p -> addVariable (false, p -> domain ()); // true if integer

  // add objective function(s) 
  expression *expr = NULL; 
  // fill in the objective
  p -> addObjective (expr, "min");  // "max" for maximization

  // add constraints:

  /*
    addEQConstraint  (expr, new exprConst (ub));  // for equality
    addLEConstraint  (expr, new exprConst (ub));  // for <=
    addGEConstraint  (expr, new exprConst (lb));  // for >=
    addRNGConstraint (expr, new exprConst (lb), new exprConst (ub));  // for range
  */

  // create room for problem's variables and bounds
  CouNumber 
    *x  = (CouNumber *) malloc (n_var * sizeof (CouNumber)),
    *lb = (CouNumber *) malloc (n_var * sizeof (CouNumber)),
    *ub = (CouNumber *) malloc (n_var * sizeof (CouNumber));

  for (int i = n_var; i--;) {
    x  [i] =  0.;
    lb [i] = -COUENNE_INFINITY;
    ub [i] =  COUENNE_INFINITY;
  }

  // create domain point for Couenne
  p -> domain () -> push (n_var, x, lb, ub);

  free (x); free (lb); free (ub);

  // fill in lower and upper bounds ///////////////////////////////////////////////////////////////

  for (register int i=n_var; i--;) {
    p -> Lb (i) = - COUENNE_INFINITY;
    p -> Ub (i) =   COUENNE_INFINITY;
    p -> X  (i) = 0;
  }

  // initial x ////////////////////////////////////////////////////////////////////

  return p;
}
