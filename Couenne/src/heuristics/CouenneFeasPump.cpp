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

#include "CouenneFeasPump.hpp"
#include "BonCouenneInterface.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "BonAuxInfos.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

// Solve //////////////////////////////////////////////////////// 
int CouenneFeasPump::solution (double & objectiveValue, double * newSolution) {

  if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
    return 0;

  // clones lower bounding problem (LP)
  //OsiSolverInterface * solver = model_ -> solver ();

  // OsiAuxInfo * auxInfo = solver -> getAuxiliaryInfo ();
  // Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (auxInfo);

  // if (babInfo) {

  //   babInfo -> setHasNlpSolution (false);

  //   // avoid doing this if we already found the problem to be infeasible
  //   if (babInfo -> infeasibleNode ())
  // 	return 0;
  // }

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
    *nSol = getContSolution (),
    *iSol = NULL,
    *best = NULL;

  int niter = 0;
  bool foundSol = false;
  CouNumber 
    cutoff = problem_ -> getCutOff (),
    bestZ  = COIN_DBL_MAX;

  /////////////////////////////////////////////////////////////////////////
  //                      _                   _                   
  //                     (_)                 | |                  
  //  _ __ ___     __ _   _    _ __          | |    ___     ___    _ __          
  // | '_ ` _ \   / _` | | |  | '_ \         | |   / _ \   / _ \  | '_ \         
  // | | | | | | | (_| | | |  | | | |        | |  | (_) | | (_) | | |_) |        
  // |_| |_| |_|  \__,_| |_|  |_| |_|        |_|   \___/   \___/  | .__/         
  //						      	          | |            
  /////////////////////////////////////////////////////////////// |_| /////

  do {

    // INTEGER PART /////////////////////////////////////////////////////////

    // Pointer to a solution in the pool with a "good solution" (if
    // nothing found in this call, return the most promising one)
    double z = getMILPSolution (iSol, nSol); 

    if (problem_ -> checkNLP (iSol, z)) {

      if (z < problem_ -> getCutOff ()) {

	foundSol = true;
	problem_ -> setCutOff (z);
	bestZ = z;
	best  = iSol;

	// if bound tightening clears the whole feasible set, stop
	if (!(problem_ -> btCore (NULL)))
	  return -1;
      }
    } else {

      OsiCuts cs;
      couenneCG_ -> genRowCuts (*nlp_, cs, 0, NULL); // remaining three argument at NULL by default
      milp_ -> applyCuts (cs);

      if (milpCuttingPlane_)
	continue; // found linearization cut, now re-solve MILP (not quite a FP)
    }

    // NONLINEAR PART ///////////////////////////////////////////////////////

    // set new objective
    expression *newObj = updateNLPObj (iSol);
    nlp_ -> setObj (0, newObj);

    // compute H_2-closest NLP feasible solution
    nlp_ -> setInitSol (iSol);
    z = nlp_ -> solve (nSol); // TODO: not necessary! keep cutting
    // integer solution with nlp cuts
    // until MINLP feasible

    // check if newly found NLP solution is also integer (unlikely...)
    if (problem_ -> checkNLP (nSol, z) &&
	(z < problem_ -> getCutOff ())) {

      foundSol = true;
      problem_ -> setCutOff (z);
      bestZ = z;
      best  = nSol;

      // if bound tightening clears the whole feasible set, stop
      if (!(problem_ -> btCore (NULL)))
	return -1;
    }	

  } while ((niter++ < maxIter_) && 
	   (!foundSol)          &&
	   (bestZ >= cutoff - COUENNE_EPS));

  // OUT OF THE LOOP: restore original objective and reoptimize

  // fix integer variables
  // resolve NLP

  if (foundSol) {

    nlp_ -> setObj (0, originalObjective_);

    fixIntVariables (best);
    nlp_ -> setInitSol (best);
    double z = nlp_ -> solve (nSol);

    // check if newly found NLP solution is also integer (unlikely...)
    if (problem_ -> checkNLP (nSol, z) &&
	(z < problem_ -> getCutOff ())) {

      problem_ -> setCutOff (z);
      bestZ = z;
      best  = nSol;

      // if bound tightening clears the whole feasible set, stop
      if (!(problem_ -> btCore (NULL)))
	return -1;
    }	
  }

  return foundSol;
}


/// find integer (possibly NLP-infeasible) point isol closest
/// (according to the l-1 norm of the hessian) to the current
/// NLP-feasible (but fractional) solution nsol
CouNumber CouenneFeasPump::getMILPSolution (CouNumber *iSol, CouNumber *nSol) {

  return 0.;
}

/// obtain continuous (if fractional) solution
CouNumber *CouenneFeasPump::getContSolution () {

  const double *solution = nlp_ -> getColSolution ();

  CouNumber *sol = CoinCopyOfArray (solution, problem_ -> nVars ());
  return sol;
}
