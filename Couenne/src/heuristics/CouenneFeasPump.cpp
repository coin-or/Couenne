/* $Id$
 *
 * Name:    CouenneFeasPump.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the Feasibility Pump heuristic class
 * Created: August 5, 2009
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneTNLP.hpp"

#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

// Solve
int CouenneFeasPump::solution (double &objVal, double *newSolution) {

  printf ("FP ====================================\n");

  if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
    return 0;

  // This FP works as follows:
  //
  // obtain current NLP solution xN or, if none available, current LP solution
  //
  // repeat {
  //
  //   1) compute MILP solution(s) {xI} H_1-closest to xN
  //   2) insert them in pool
  //   3) consider the most promising one xI in the whole pool
  //
  //   if xI is MINLP-feasible (hack incumbent callback)
  //     update cutoff with min obj among MINLP-feasible
  //     run BT
  //   [else // a little out of the FP scheme
  //     apply linearization cuts to MILP]
  //
  //   select subset of variables (cont or integer), fix them
  //
  //   compute NLP solution xN that is H_2-closest to xI
  //
  //   if MINLP feasible
  //     update cutoff with min obj among MINLP-feasible
  //     run BT
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

  expression *origObj = problem_ -> Obj (0) -> Body ();

  int objInd = problem_ -> Obj (0) -> Body () -> Index ();

  do {

    printf ("FP: main loop\n");

    // INTEGER PART /////////////////////////////////////////////////////////

    // Solve IP using nSol as the initial point to minimize weighted
    // l-1 distance from. If nSol==NULL, the MILP is created using the
    // original milp's LP solution.
    double z = solveMILP (nSol, iSol);

    bool isChecked = false;
#ifdef FM_CHECKNLP2
    isChecked = problem_->checkNLP2(iSol, 0, false, // do not care about obj 
				    true, // stopAtFirstViol
				    true, // checkALL
				    problem_->getFeasTol());
    if(isChecked) {
      z = problem_->getRecordBestSol()->getModSolVal();
    }
#else /* not FM_CHECKNLP2 */
    isChecked = problem_ -> checkNLP (iSol, z, true);
#endif  /* not FM_CHECKNLP2 */

    if (isChecked) {

      printf ("FP: solution found is MINLP-feasible! Comparing now %g with cutoff %g\n", 
	      z, problem_ -> getCutOff ());

      // solution is MINLP feasible! 
      // Save 

      if (z < problem_ -> getCutOff ()) {

	printf ("FP: (and it is better than the cutoff)\n");

	retval = 1;
	objVal = z;

#ifdef FM_CHECKNLP2
#ifdef FM_TRACE_OPTSOL
	problem_->getRecordBestSol()->update();
	best = problem_->getRecordBestSol()->getSol();
	objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
	best = problem_->getRecordBestSol()->getModSol();
	objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#else /* not FM_CHECKNLP2 */
#ifdef FM_TRACE_OPTSOL
	problem_->getRecordBestSol()->update(iSol, problem_->nVars(),
					     z, problem_->getFeasTol());
	best = problem_->getRecordBestSol()->getSol();
	objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
	best   = iSol;
	objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */

	problem_ -> setCutOff (objVal);

	t_chg_bounds *chg_bds = NULL;

	if (objInd >= 0) {
	
	  chg_bds = new t_chg_bounds [problem_ -> nVars ()];
	  chg_bds [objInd].setUpper (t_chg_bounds::CHANGED); 
	}

	// if bound tightening clears the whole feasible set, stop
	bool is_still_feas = problem_ -> btCore (chg_bds);

	if (chg_bds) 
	  delete [] chg_bds;

	if (!is_still_feas)
	  break;
      }
    } else {

      printf ("FP: MINLP-infeasible, try some cuts\n");

      // solution non MINLP feasible, it might get cut by
      // linearization cuts. If so, add a round of cuts and repeat.

      OsiCuts cs;
      // remaining three arguments at NULL by default
      couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); 

      if (cs.sizeRowCuts ()) { 

	printf ("FP: found cuts\n");
    
	// the (integer, NLP infeasible) solution could be separated

	milp_ -> applyCuts (cs);

	if (milpCuttingPlane_)
	  continue; // found linearization cut, now re-solve MILP (not quite a FP)
      }
    }

    // NONLINEAR PART ///////////////////////////////////////////////////////

    // fix some variables and solve the NLP to find a NLP (possibly
    // non-MIP) feasible solution

    z = solveNLP (iSol, nSol); 

    // check if newly found NLP solution is also integer

    isChecked = false;
#ifdef FM_CHECKNLP2
    isChecked = problem_->checkNLP2(nSol, 0, false, // do not care about obj 
				    true, // stopAtFirstViol
				    true, // checkALL
				    problem_->getFeasTol());
    if(isChecked) {
      z = problem_->getRecordBestSol()->getModSolVal();
    }
#else /* not FM_CHECKNLP2 */
    isChecked = problem_ -> checkNLP (nSol, z, true);
#endif  /* not FM_CHECKNLP2 */

    if (isChecked &&
	(z < problem_ -> getCutOff ())) {

      retval = 1;

#ifdef FM_CHECKNLP2
#ifdef FM_TRACE_OPTSOL
      problem_->getRecordBestSol()->update();
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
      best = problem_->getRecordBestSol()->getModSol();
      objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#else /* not FM_CHECKNLP2 */
#ifdef FM_TRACE_OPTSOL
      problem_->getRecordBestSol()->update(nSol, problem_->nVars(),
					   z, problem_->getFeasTol());
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
      best   = nSol;
      objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */
      
      problem_ -> setCutOff (objVal);

      t_chg_bounds *chg_bds = NULL;

      if (objInd >= 0) {
	
	chg_bds = new t_chg_bounds [problem_ -> nVars ()];
	chg_bds [objInd].setUpper (t_chg_bounds::CHANGED); 
      }

      // if bound tightening clears the whole feasible set, stop
      bool is_still_feas = problem_ -> btCore (chg_bds);

      if (chg_bds) 
	delete [] chg_bds;

      if (!is_still_feas)
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

    printf ("FP: final NLP\n");

    if (!nlp_) // first call (in this call to FP). Create NLP
      nlp_ = new CouenneTNLP (problem_);

    problem_ -> setObjective (0, origObj);

    fixIntVariables (best);
    nlp_ -> setInitSol (best);

    ////////////////////////////////////////////////////////////////

    // shamelessly copied from hs071_main.cpp (it's Open Source too!)

    SmartPtr <IpoptApplication> app = IpoptApplicationFactory ();

    ApplicationReturnStatus status = app -> Initialize ();
    if (status != Solve_Succeeded) printf ("FP: error in initialization\n");

    status = app -> OptimizeTNLP (nlp_);
    if (status != Solve_Succeeded) printf ("FP: error solving problem\n");

    ////////////////////////////////////////////////////////////////

    double z = nlp_ -> getSolValue ();

    problem_ -> domain () -> pop (); // pushed in fixIntVariables

    // check if newly found NLP solution is also integer (unlikely...)
    bool isChecked = false;
#ifdef FM_CHECKNLP2
    isChecked = problem_->checkNLP2(nSol, 0, false, // do not care about obj 
				    true, // stopAtFirstViol
				    true, // checkALL
				    problem_->getFeasTol());
    if(isChecked) {
      z = problem_->getRecordBestSol()->getModSolVal();
    }
#else /* not FM_CHECKNLP2 */
    isChecked = problem_ -> checkNLP (nSol, z, true);
#endif  /* not FM_CHECKNLP2 */

    if (isChecked &&
	(z < problem_ -> getCutOff ())) {

#ifdef FM_CHECKNLP2
#ifdef FM_TRACE_OPTSOL
      problem_->getRecordBestSol()->update();
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
      best = problem_->getRecordBestSol()->getModSol();
      objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#else /* not FM_CHECKNLP2 */
#ifdef FM_TRACE_OPTSOL
      problem_->getRecordBestSol()->update(nSol, problem_->nVars(),
					   z, problem_->getFeasTol());
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
      best   = nSol;
      objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */

      problem_ -> setCutOff (objVal);
    }	
  }

  if (retval > 0) 
    CoinCopyN (best, problem_ -> nVars (), newSolution);

  delete [] iSol;
  delete [] nSol; // best is either iSol or nSol, so don't delete [] it

  // release bounding box
  problem_ -> domain () -> pop (); // pushed in first call to solveMILP

  // milp is deleted at every call since it changes not only in terms
  // of variable bounds but also in terms of linearization cuts added
  delete milp_;

  return retval;
}
