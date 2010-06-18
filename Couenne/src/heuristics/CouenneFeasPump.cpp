/* $Id$
 *
 * Name:    CouenneFeasPump.cpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the Feasibility Pump heuristic class
 * Created: August 5, 2009
 * 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneTNLP.hpp"

using namespace Couenne;

// Solve
int CouenneFeasPump::solution (double &objVal, double *newSolution) {

  if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
    return 0;

  // This FP works as follows:
  //
  // obtain current NLP solution xN or, if none available, current LP solution
  //
  // repeat {
  //
  //   1) compute MILP solution(s) xI that is H_1-closest to xN
  //   2) insert them in pool
  //   3) consider the most promising one in the whole pool
  //
  //   if it is MINLP-feasible (hack incumbent callback)
  //     update cutoff with min obj among MINLP-feasible
  //     run BT
  //   else
  //     apply linearization cuts to MILP
  //
  //   select subset of variables (cont or integer), fix them
  //
  //   compute NLP solution xN that is H_2-closest to xI
  //
  //   if MINLP-feasible
  //     update cutoff
  //
  // } until exit condition satisfied
  //
  // fix integer variables
  // resolve NLP

  CouNumber 
    *nSol = NULL, // solution of the nonlinear problem
    *iSol = NULL, // solution of the IP problem
    *best = NULL; // best solution found so far

  int
    niter  = 0, // current # of iterations
    retval = 0; // 1 if found a better solution

  /////////////////////////////////////////////////////////////////////////
  //                      _                   _
  //                     (_)                 | |
  //  _ __ ___     __ _   _    _ __          | |    ___     ___    _ __ 
  // | '_ ` _ \   / _` | | |  | '_ \         | |   / _ \   / _ \  | '_ `.
  // | | | | | | | (_| | | |  | | | |        | |  | (_) | | (_) | | |_) |
  // |_| |_| |_|  \__,_| |_|  |_| |_|        |_|   \___/   \___/  | .__/
  //						      	          | |
  /////////////////////////////////////////////////////////////// |_| /////

  milp_ = model_ -> solver () -> clone ();

  expression *origObj = problem_ -> Obj (0) -> Body ();

  // copy bounding box and current solution to the problem (better
  // linearization cuts)
  problem_ -> domain () -> push (milp_ -> getNumCols (),
				 milp_ -> getColSolution (),
				 milp_ -> getColLower (),
				 milp_ -> getColUpper ());
  do {

    // INTEGER PART /////////////////////////////////////////////////////////

    // Solve IP using nSol as the initial point to minimize weighted
    // l-1 distance from. If nSol==NULL, the MILP is created using the
    // original milp's LP solution.
    double z = solveMILP (nSol, iSol);

    if (problem_ -> checkNLP (iSol, z)) {

      // solution is MINLP feasible! 
      // Save and get out of the loop

      if (z < problem_ -> getCutOff ()) {

	retval = 1;
	objVal = z;
	best   = iSol;
	problem_ -> setCutOff (z);

	// if bound tightening clears the whole feasible set, stop
	if (!(problem_ -> btCore (NULL)))
	  break;
      }
    } else {

      // solution non MINLP feasible, it might get cut by
      // linearization cuts. If so, add a round of cuts and repeat.

      OsiCuts cs;
      // remaining three arguments at NULL by default
      couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); 

      if (cs.sizeRowCuts ()) { 

	// the (integer, NLP infeasible) solution could be separated

	milp_ -> applyCuts (cs);

	if (milpCuttingPlane_)
	  continue; // found linearization cut, now re-solve MILP (not quite a FP)
      }
    }

    // NONLINEAR PART ///////////////////////////////////////////////////////

    // solve the NLP to find a NLP (possibly non-MIP) feasible solution
    z = solveNLP (iSol, nSol); 

    // check if newly found NLP solution is also integer
    if (problem_ -> checkNLP (nSol, z) &&
	(z < problem_ -> getCutOff ())) {

      retval = 1;
      objVal = z;
      best   = nSol;
      problem_ -> setCutOff (z);

      // if bound tightening clears the whole feasible set, stop
      if (!(problem_ -> btCore (NULL)))
	break;
    }	

  } while ((niter++ < maxIter_) && 
	   (retval == 0));

  // OUT OF THE LOOP ////////////////////////////////////////////////////////

  // If MINLP solution found,
  //
  // 1) restore original objective 
  // 2) fix integer variables
  // 3) resolve NLP

  if (retval > 0) {

    if (!nlp_) // first call (in this call to FP). Create NLP
      nlp_ = new CouenneTNLP (problem_);

    problem_ -> setObjective (0, origObj);

    fixIntVariables (best);
    nlp_ -> setInitSol (best);

    ////////////////////////////////////////////////////////////////

    // shamelessly copied from hs071_main.cpp (it's Open Source too!)

    SmartPtr <IpoptApplication> app = IpoptApplicationFactory ();

    ApplicationReturnStatus status = app -> Initialize ();
    if (status != Solve_Succeeded) printf ("Error in initialization\n");

    status = app -> OptimizeTNLP (nlp_);
    if (status != Solve_Succeeded) printf ("Error solving problem\n");

    ////////////////////////////////////////////////////////////////

    double z = nlp_ -> solve (nSol);

    problem_ -> domain () -> pop (); // pushed in fixIntVariables

    // check if newly found NLP solution is also integer (unlikely...)
    if (problem_ -> checkNLP (nSol, z) &&
	(z < problem_ -> getCutOff ())) {

      problem_ -> setCutOff (z);
      objVal = z;
      best   = nSol;
    }	
  }

  if (retval > 0) 
    CoinCopyN (best, problem_ -> nVars (), newSolution);

  delete [] iSol;
  delete [] nSol; // best is either iSol or nSol, so don't delete [] it

  // release bounding box
  problem_ -> domain () -> pop ();

  // milp is deleted at every call since it changes not only in terms
  // of variable bounds but also in terms of linearization cuts added
  delete milp_;

  return retval;
}
