/* $Id$
 *
 * Name:    CouenneOS.cpp
 * Authors: 
 *          
 * Purpose: Creates a CouenneProblem object from an OSil instance
 *
 * (C) Carnegie-Mellon University, 2009. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneOSInterface.hpp"

#include "CouenneProblem.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExprSum.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprGroup.hpp"

#include "BonTMINLP.hpp"

//#include "OSInstance.hpp"
//class OSInstance;

namespace Bonmin {

  class RegisteredOptions;
  class Bab;
  class TMINLP;
}

using Ipopt::SmartPtr;

using namespace Couenne;

void CouenneOSInterface::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions) {
	roptions->AddStringOption1("osilfile", "name of an osil file to read the problem instance from", "", "*", "name of osil file");
}

CouenneOSInterface::~CouenneOSInterface() {
	//delete osinstance;
	delete problem;
}

CouenneProblem* CouenneOSInterface::getCouenneProblem() {
	if (!osinstance) {
		// create osinstance from osilfile
	}

  problem = new CouenneProblem;

  //p -> setProblemName (filename); // sets filename member, for later stats -- TODO

  // number of defined variables (aka common expressions)
  //ndefined_ = 0; // Don't know if OS has them: in AMPL they are
                 // defined by "var y := f(x)", with f(x) some
		 // expression

  int n_var = 0; // to be set to no. of variables

  // nonlinear in both objectives and constraints
  for (int i = 0; i < n_var; i++) 
    problem -> addVariable (false, problem -> domain ()); // true if integer

  // add objective function(s) 
  expression *expr = NULL; 
  // fill in the objective
  problem -> addObjective (expr, "min");  // "max" for maximization

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
  problem -> domain () -> push (n_var, x, lb, ub);
  free (x); free (lb); free (ub);

  // suggested:
  // problem -> domain () -> push (n_var + problem -> nDefVars(), x, lb, ub, false);
  // //free (x); free (lb); free (ub);
  // saves three allocations (default last parameter is true, which copies x,l,b)

  // fill in lower and upper bounds ///////////////////////////////////////////////////////////////

  for (register int i=n_var; i--;) {
    problem -> Lb (i) = - COUENNE_INFINITY;
    problem -> Ub (i) =   COUENNE_INFINITY;
    problem -> X  (i) = 0;
  }

  // initial x ////////////////////////////////////////////////////////////////////

  return problem;
}

Ipopt::SmartPtr<Bonmin::TMINLP> CouenneOSInterface::getTMINLP() {
	return tminlp;
}

bool CouenneOSInterface::writeSolution(Bonmin::Bab& bab) {

	return false;
}
