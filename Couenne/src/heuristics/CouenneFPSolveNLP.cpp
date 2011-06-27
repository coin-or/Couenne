/* $Id$
 *
 * Name:    CouenneFPSolveNLP.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the NLP solution method for the Feasibility Pump 
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"

#include "CouenneTNLP.hpp"

using namespace Ipopt;
using namespace Couenne;

/// obtain continuous (if fractional) solution
CouNumber CouenneFeasPump::solveNLP (CouNumber *iSol, CouNumber *&nSol) {

  // Solve the continuous nonlinear programming problem
  //
  // min  f(x)
  // s.t. g(x) <= 0
  //
  // where g(x) are the original constraints and f(x) is one of the
  // following:
  //
  // 1) sum {i in Vars} (x_i - x_i^0)^2
  // 2) sum {i in I}    (x_i - x_i^0)^2
  // 3) sum {i in Vars} (P^i (x - x^0))^2
  // 4) sum {i in I}    (P^i (x - x^0))^2
  //
  // where is x^0 is the optimal solution of a MILP problem. P should
  // be a PSD matrix, but the Hessian is, in general, indefinite at
  // the IP point we are starting from. A cheap convexification
  // consists of computing the minimum eigenvalue lambda_min of H and,
  // if lambda_min < 0, replace H with
  //
  // H - lambda_min I
  //
  // Similarly to the MILP case, we have
  //
  // P = beta I + (1-beta) (H + lambda_min I) 
  //   = (beta + lambda_min (1 - beta)) I + (1-beta) H

  bool firstNLP = (nlp_ == NULL);

  if (firstNLP) // first call (in this call to FP). Create NLP
    nlp_ = new CouenneTNLP (problem_);

  problem_ -> domain () -> push (problem_ -> nVars (),
				 iSol,
				 NULL, // replaces problem_ -> domain () -> lb (),
				 NULL); // replaces problem_ -> domain () -> ub (),
				 //false); // to avoid overlapping with nsol within NLP

  // set new objective
  expression
    *oldObj = problem_ -> Obj (0) -> Body (),
    *newObj = updateNLPObj (iSol);

  newObj   -> realign (problem_);
  problem_ -> setObjective (0, newObj);
  nlp_     -> setObjective (newObj);

  if (problem_ -> Jnlst () -> ProduceOutput (J_STRONGWARNING, J_NLPHEURISTIC)) {
    printf ("now solving NLP:\n");
    problem_ -> print ();
    printf ("-----------------------\n");
  }

  // FIXME: probably the previous NLP optimum is a better starting point

  // compute H_2-closest NLP feasible solution
  nlp_ -> setInitSol (iSol);

  /////////////////////////////////////////////////////////

  ApplicationReturnStatus status = firstNLP ? 
    app_ -> OptimizeTNLP   (nlp_) :
    app_ -> ReOptimizeTNLP (nlp_);

  /////////////////////////////////////////////////////////

  if (nlp_ -> getSolution ()) // check if non-NULL

    if  (nSol)  CoinCopyN       (nlp_ -> getSolution (), problem_ -> nVars (), nSol);
    else nSol = CoinCopyOfArray (nlp_ -> getSolution (), problem_ -> nVars ());

  else problem_ -> Jnlst () -> Printf 
      (J_ERROR, J_NLPHEURISTIC, "FP: warning, NLP returns a NULL solution\n");

  delete newObj;

  CouNumber retval;

  problem_ -> setObjective (0, oldObj);

  if ((status != Solve_Succeeded) &&
      (status != Solved_To_Acceptable_Level))

    problem_ -> Jnlst () -> Printf 
      (J_ERROR, J_NLPHEURISTIC, "Feasibility Pump: Error solving NLP problem\n");

  retval = nlp_ -> getSolValue ();

  problem_ -> domain () -> pop ();

  return retval;
}
