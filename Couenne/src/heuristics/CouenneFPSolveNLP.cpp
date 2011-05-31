/* $Id$
 *
 * Name:    CouenneFPSolveNLP.cpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the NLP solution method for the Feasibility Pump 
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "IpIpoptApplication.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

#include "CouenneTNLP.hpp"

using namespace Couenne;

/// obtain continuous (if fractional) solution
CouNumber CouenneFeasPump::solveNLP (CouNumber *iSol, CouNumber *&nSol) {

  printf ("FP: solveNLP\n");

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
  // be a PSD matrix, but the Hessian is, in general, indefinite. A
  // cheap convexification consists of computing the minimum
  // eigenvalue lambda_min of H and, if lambda_min < 0, replace H with
  //
  // H - lambda_min I
  //
  // Similarly to the MILP case, we have
  //
  // P = beta I + (1-beta) (H + lambda_min I) 
  //   = (beta + lambda_min (1 - beta)) I + (1-beta) H

  bool firstNLP = (nlp_ == NULL);

  if (!nlp_) // first call (in this call to FP). Create NLP
    nlp_ = new CouenneTNLP (problem_);

  problem_ -> domain () -> push (problem_ -> nVars (),
				 iSol,
				 NULL, //problem_ -> domain () -> lb (),
				 NULL, //problem_ -> domain () -> ub (),
				 false);

  printf ("FP: created TNLP\n");

  // set new objective
  expression *newObj = updateNLPObj (iSol);
  problem_ -> setObjective (0, newObj);

  printf ("FP: created obj\n");

  // compute H_2-closest NLP feasible solution
  nlp_ -> setInitSol (iSol);

  /////////////////////////////////////////////////////////

  // shamelessly copied from hs071_main.cpp (it's Open Source too!)

  SmartPtr <IpoptApplication> app = IpoptApplicationFactory ();

  ApplicationReturnStatus status = app -> Initialize ();

  if (status != Solve_Succeeded)
    printf ("FP: Error in initialization\n");

  printf ("FP: optimize\n");

  status = firstNLP ? 
    app -> OptimizeTNLP   (nlp_) :
    app -> ReOptimizeTNLP (nlp_);

  problem_ -> domain () -> pop ();

  if (status != Solve_Succeeded)
    printf ("FP: Error solving problem\n");

  /////////////////////////////////////////////////////////

  if (nlp_ -> getSolution ()) // check if non-NULL
    nSol = CoinCopyOfArray (nlp_ -> getSolution (), problem_ -> nVars ());

  // integer solution with nlp cuts
  // until MINLP feasible

  delete newObj;

  return nlp_ -> getSolValue ();
}
