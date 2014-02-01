/* $Id$
 *
 * Name:    CouenneFPphaseNLP.cpp
 * Authors: Pietro Belotti
 * Purpose: NLP part of the loop
 * Created: February 1, 2014
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneExprAux.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneTNLP.hpp"
#include "CouenneFPpool.hpp"
#include "CouenneRecordBestSol.hpp"

#ifdef COIN_HAS_SCIP
/* general SCIP includes */
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#endif

using namespace Ipopt;
using namespace Couenne;

void printDist   (CouenneProblem *p, const double *iSol, double *nSol);
void printCmpSol (CouenneProblem *p, const double *iSol, double *nSol, int direction);

// Solve
int CouenneFeasPump::solution (double &objVal, double *newSolution) {

  const int depth = (model_ && (model_ -> currentNode ())) ? model_ -> currentNode () -> depth () : 0;

  double time0 = CoinCpuTime();

  if ((nCalls_ == 0) ||                                   // check upper limit on number of runs
      //(problem_ -> nIntVars () <= 0) ||                 // feas pump on NLP? Why not?
      (CoinCpuTime () >= problem_ -> getMaxCpuTime ()) || // don't start if time is out
      ((numberSolvePerLevel_ >= 0) &&                     // stop FP after a certain level
       (CoinDrand48 () > 1. / CoinMax                     // decided randomly and inversely proportional
	(1., (double) ((CoinMax (0, depth - numberSolvePerLevel_ + 1)) *   // to BB tree depth
   		                   (depth - numberSolvePerLevel_ + 1))))))
    return 0;

  if (nCalls_ > 0) // only decrease it if positive. If -1, no limit on # calls
    --nCalls_;

  problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "==================================================== FP: BEGIN\n");
  problem_ -> Jnlst () -> Printf (J_ERROR,   J_NLPHEURISTIC, "[FeasPump] Initializing\n");

  // Solve NLP once at the beginning ////////////////////////
  //
  // First, it possibly provides an NLP feasible point (the LP
  // provides a neither NLP nor integer point)
  //
  // Second, we need it to compile the Lagrangian Hessian used as a
  // distance

  if (!nlp_) // first call (in this run of FP). Create NLP
    nlp_ = new CouenneTNLP (problem_);

  problem_ -> domain () -> push (*(problem_ -> domain () -> current ()));

  // Initial Bound Tightening: since Cbc prefers (reasonably) to call
  // a heuristic before any cut separator and because BT is instead
  // quite fast, let's see if we can reduce bounds to get a tighter
  // MILP relaxation.

  t_chg_bounds *chg_bds = new t_chg_bounds [problem_ -> nVars ()];

  for (int i=problem_ -> nVars (); i--;) {
    chg_bds [i]. setLower (t_chg_bounds::CHANGED);
    chg_bds [i]. setUpper (t_chg_bounds::CHANGED);
  }

  bool is_still_feas = problem_ -> btCore (chg_bds);

  delete [] chg_bds;

  // Now that bounds are possibly reduced, put current point within
  // bounds

  double
    *lb  = problem_ -> domain () -> lb (),
    *ub  = problem_ -> domain () -> ub (),
    *sol = problem_ -> domain () -> x  ();

  for (int i=problem_ -> nVars (); i--;)

    if      (sol [i] < lb [i]) sol [i] = lb [i];
    else if (sol [i] > ub [i]) sol [i] = ub [i];

  // fix integer coordinates of current (MINLP feasible!) solution
  // and set it as initial (obviously NLP feasible) solution

  nlp_ -> setInitSol (sol);

  ////////////////////////////////////////////////////////////////

  if ((multHessNLP_  != 0.) || 
      (multHessMILP_ != 0.))
    nlp_ -> getSaveOptHessian () = true;

  if (problem_ -> Jnlst () -> ProduceOutput (J_ALL, J_NLPHEURISTIC)) {
    printf ("initial: solving NLP:---------\n");
    problem_ -> print ();
    printf ("---------------------\n");
  }

  problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: Initial NLP... "); fflush (stdout);

  // Solve with original objective function
  ApplicationReturnStatus status = app_ -> OptimizeTNLP (nlp_);

  problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "done\n"); fflush (stdout);

  if (nlp_ -> getSolution () && (problem_ -> Jnlst () -> ProduceOutput (J_ALL, J_NLPHEURISTIC))) { // check if non-NULL
    printf ("######################## NLP solution (init):\n");
    for (int i=0; i< problem_ -> nVars ();) {
      printf ("%+e ", nlp_ -> getSolution () [i]);
      if (!(++i % 15)) printf ("\n");
    }
    if ((problem_ -> nVars () - 1) % 5) printf ("\n");
  }

  ////////////////////////////////////////////////////////////////

  problem_ -> domain () -> pop ();

  if ((status != Solve_Succeeded) && (status != Solved_To_Acceptable_Level))
    problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "Feasibility Pump: error in initial NLP problem\n");

  if ((multHessNLP_  != 0.) ||
      (multHessMILP_ != 0.))
    nlp_ -> getSaveOptHessian () = false;

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
    niter    = 0,   // iteration counter
    nsuciter = 0,   // counter for consecutive successful applications of one MILP solving method
    retval   = 0,   // 1 if found a better solution
    nSep     = 0,   // If separation short circuit, max # of consecutive separations
    objInd   = problem_ -> Obj (0) -> Body () -> Index ();

  /////////////////////////////////////////////////////////////////////////
  //                      _                   _
  //                     (_)                 | |
  //  _ __ ___     __ _   _    _ __          | |    ___     ___    _ __ 
  // | '_ ` _ \   / _` | | |  | '_ \         | |   / _ \   / _ \  | '_ `.
  // | | | | | | | (_| | | |  | | | |        | |  | (_) | | (_) | | |_) |
  // |_| |_| |_|  \__,_| |_|  |_| |_|        |_|   \___/   \___/  | .__/
  //						      	          | |
  /////////////////////////////////////////////////////////////// |_| /////

  // copy bounding box and current solution to the problem (better
  // linearization cuts). Note that this push is pop()'d at the end
  // of the main routine

  problem_ -> domain () -> push (model_ -> solver ());

  expression *originalObjective = problem_ -> Obj (0) -> Body ();

  double
    save_mDN = multDistNLP_,
    save_mHN = multHessNLP_,
    save_mON = multObjFNLP_,
    save_mDM = multDistMILP_,
    save_mHM = multHessMILP_,
    save_mOM = multObjFMILP_;

  do {

    if (niter) {
      multDistNLP_  *= fabs (save_mDN); // update within loop
      multHessNLP_  *= fabs (save_mHN);
      multObjFNLP_  *= fabs (save_mON);
      multDistMILP_ *= fabs (save_mDM);
      multHessMILP_ *= fabs (save_mHM);
      multObjFMILP_ *= fabs (save_mOM);
    }

    if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
      break;

    problem_ -> Jnlst () -> Printf (J_ERROR, J_NLPHEURISTIC, "[FeasPump] Iteration %d [%gs]\n", niter, CoinCpuTime() - time0);

    // INTEGER PART /////////////////////////////////////////////////////////

    // Solve IP using nSol as the initial point to minimize weighted
    // l-1 distance from. If nSol==NULL, the MILP is created using the
    // original milp's LP solution.
            
    double z = solveMILP (nSol, iSol, niter, &nsuciter);

    // if no MILP solution was found, bail out

    bool try_again = false;

    if (!iSol || z >= COIN_DBL_MAX/2) {

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: could not find IP solution\n");

      // find non-tabu solution in the solution pool
      while (!pool_ -> Set (). empty ()) {

	// EXTRACT the closest (to nSol) IP solution from the pool
	pool_ -> findClosestAndReplace (iSol, iSol, problem_ -> nVars ());

	CouenneFPsolution newSol (problem_, iSol);

	// we found a solution that is not in the tabu list
	if (tabuPool_ . find (newSol) == tabuPool_ . end ()) {

	  try_again = true;
	  break;
	}
      }

      if (!try_again) { // try moving around current solution 

	// SCIP could not find a MILP solution and we're somewhat
	// locked. Round current solution and feed it to the NLP. This
	// is better than just bailing out.

	int n = problem_ -> nVars ();

	if (!iSol)
	  iSol = new double [n];

	for (int i=0; i<n; i++)
	  iSol [i] = (problem_ -> Var (i) -> isInteger ()) ? 
	    COUENNE_round (nSol [i]) : 
	    nSol [i];
      }
 
      // if (!try_again) { // nothing to do, bail out	
      // 	problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: could not find from pool either, bailing out\n");
      // 	break;
      // }
    }

    bool isChecked = false;

    // If a solution was found, but is in the tabu list, two choices:
    //
    // 1) the pool is empty: do a round of cuts and repeat;
    //
    // 2) the pool is nonempty: extract the best solution from the
    //    pool and use it instead

    CouenneFPsolution checkedSol (problem_, iSol, false); // false is for not allocating space for this

    if (tabuPool_. find (checkedSol) == tabuPool_ . end ())

      tabuPool_. insert (CouenneFPsolution (problem_, iSol)); // only insertion to tabu pool: we check its feasibility now

    else {

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: found solution is tabu\n");

      // Current solution was visited before. Replace it with another
      // MILP solution from the pool, if any.

      if         (tabuMgt_ == FP_TABU_NONE) break;

      else   if ((tabuMgt_ == FP_TABU_POOL) && !(pool_ -> Set (). empty ())) {

	// try to find non-tabu solution in the solution pool
	do {

            // retrieve the top solution from the pool
            pool_ -> findClosestAndReplace (iSol, nSol, problem_ -> nVars ());

	    CouenneFPsolution newSol (problem_, iSol);

            // we found a solution that is not in the tabu list
            if (tabuPool_. find (newSol) == tabuPool_ . end ())
	      break;

            // the pool is empty -> bail out
            if (pool_ -> Set ().empty ())
	      {
		delete[] iSol;
		iSol = NULL;
	      }

	} while( !pool_ -> Set ().empty() );

      } else if (((tabuMgt_ == FP_TABU_CUT)   ||  
		  ((pool_ -> Set (). empty ()) && iSol))) {

	OsiCuts cs;

	problem_   -> domain () -> push (problem_ -> nVars (), iSol, milp_ -> getColLower (), milp_ -> getColUpper (), false); // don't copy vectors
	couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); // remaining three arguments NULL by default
	problem_   -> domain () -> pop ();

	if (cs.sizeRowCuts ()) {

	  milp_ -> applyCuts (cs);

	  if (nSep++ < nSepRounds_)
	    continue;

	} else break; // nothing left to do, just bail out

      } else if ((tabuMgt_ == FP_TABU_PERTURB) && iSol) {

	// perturb solution	

	const CouNumber 
	  *lb = milp_ -> getColLower (),
	  *ub = milp_ -> getColUpper ();

	double 
	  downMoves = 0.,
	  upMoves   = 0.;

	int startIndex = (int) floor (CoinDrand48 () * problem_ -> nOrigVars ());

	for (int ii = problem_ -> nOrigVars (); ii--; lb++, ub++) {

	  if (problem_ -> Var (ii) -> Multiplicity () <= 0)
	    continue;

	  // make perturbation start from pseudo-random points

	  int i = (startIndex + ii) % problem_ -> nOrigVars ();

	  if (problem_ -> Var (i) -> isInteger ()) {

	    double
	      rnd  = CoinDrand48 (), 
	      down = 0.,
	      up   = 1.;

	    // if there is room on the left (resp. right) of the
	    // solution, consider moving down (resp. up). Make moves
	    // less likely as we don't want to perturb too many
	    // variables

#define RND_DECR_EXPONENT .5

	    if (iSol [i] >= lb [i] + 1.) down =      1. / pow (1. + (downMoves += 1.), RND_DECR_EXPONENT);
	    if (iSol [i] <= ub [i] - 1.) up   = 1. - 1. / pow (1. + (upMoves   += 1.), RND_DECR_EXPONENT);

	    if      (rnd < down) iSol [i] -= 1.;
	    else if (rnd > up)   iSol [i] += 1.;
	  }
	}
      }
    } 

    problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: checking IP solution for feasibility\n");

    isChecked = problem_ -> checkNLP0 (iSol, z, true, 
				       false, // don't care about obj
				       true,  // stop at first violation
				       true); // checkAll

    // Possible improvement: if IP-infeasible or if z does not improve
    // the best solution found so far, then try the following:
    //
    // 1) Fix integer variables to IP solution's corresponding components
    // 2) Solve restriction with original obj
    // 3) While not found (MI)NLP solution or
    //          solution MINLP infeasible or
    //          z not better
    // 4)   Get solution x* from pool
    // 5)   Solve with original objective
    // 6) Done
    // 7) If found solution, set it, otherwise keep previously found (IP-feasible) one

    if ((!isChecked ||                     // not MINLP feasible
	 z > problem_ -> getCutOff ())) {  // not improving

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: infeasible/non-improving (feas==%d, z=%g, cutoff=%g), looping on pool\n", isChecked, z, problem_->getCutOff ());

      bool try_again;

      do {

	try_again = false;

	if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	  break;

	// Check if fixing integers to those in iSol yields an
	// infeasible problem. If so, don't optimize
	if (fixIntVariables (iSol)) {

	  nlp_ -> setInitSol (iSol);

	  if (problem_ -> Jnlst () -> ProduceOutput (J_ALL, J_NLPHEURISTIC)) {
	    printf ("----------------------- Solving NLP:\n");
	    problem_ -> print ();
	    printf ("-----------------------\n");
	  }

	  status = app_ -> OptimizeTNLP (nlp_);

	  if (nlp_ -> getSolution ()) { // check if non-NULL
	    if  (nSol)  CoinCopyN       (nlp_ -> getSolution (), problem_ -> nVars (), nSol);
	    else nSol = CoinCopyOfArray (nlp_ -> getSolution (), problem_ -> nVars ());
	  }

	  if (nlp_ -> getSolution () && (problem_ -> Jnlst () -> ProduceOutput (J_ALL, J_NLPHEURISTIC))) { // check if non-NULL
	    printf ("######################## NLP solution (loop through pool):\n");
	      for (int i=0; i< problem_ -> nVars ();) {
		printf ("%+e ", nlp_ -> getSolution () [i]);
		if (!(++i % 15)) printf ("\n");
	      }
	  }

	  z = nlp_ -> getSolValue ();

	  if (z < problem_ -> getCutOff ()) // don't waste time if not improving

	    isChecked = problem_ -> checkNLP0 (nSol, z, true, 
					       false, // don't care about obj
					       true,  // stop at first violation
					       true); // checkAll

	  if (isChecked && (z < problem_ -> getCutOff ())) 
	    break;
	} 

	// find non-tabu solution in the solution pool
	while (!pool_ -> Set (). empty ()) {

	  // EXTRACT the closest (to nSol) IP solution from the pool
	  pool_ -> findClosestAndReplace (iSol, nSol, problem_ -> nVars ());

	  CouenneFPsolution newSol (problem_, iSol);

	  // we found a solution that is not in the tabu list
	  if (tabuPool_ . find (newSol) == tabuPool_ . end ()) {
	    try_again = true;
	    break;
	  }
	} 

      } while (try_again);
    }

    // Whatever the previous block yielded, we are now going to check
    // if we have a new solution, and if so we save it.

    if (isChecked) {

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: IP solution is MINLP feasible\n");

      // solution is MINLP feasible! Save it.

      retval = 1;
      objVal = z;

      // Found a MINLP-feasible solution, but to keep diversity do NOT
      // use the best available. Just use this.
      //
      // #ifdef FM_CHECKNLP2
      // #  ifdef FM_TRACE_OPTSOL
      //       problem_->getRecordBestSol()->update();
      //       best = problem_->getRecordBestSol()->getSol();
      //       objVal = problem_->getRecordBestSol()->getVal();
      // #  else /* not FM_TRACE_OPTSOL */
      //       best = problem_->getRecordBestSol()->getModSol(problem_->nVars());
      //       objVal = z;
      // #  endif /* not FM_TRACE_OPTSOL */
      // #else /* not FM_CHECKNLP2 */
      // #  ifdef FM_TRACE_OPTSOL
      //       problem_->getRecordBestSol()->update(iSol, problem_->nVars(),
      // 					   z, problem_->getFeasTol());
      //       best = problem_->getRecordBestSol()->getSol();
      //       objVal = problem_->getRecordBestSol()->getVal();
      // #  else /* not FM_TRACE_OPTSOL */

      best = iSol;

      //       objVal = z;
      // #  ifdef FM_TRACE_OPTSOL
      //       problem_->getRecordBestSol()->update();
      //       best = problem_->getRecordBestSol()->getSol();
      //       objVal = problem_->getRecordBestSol()->getVal();
      // #  endif /* not FM_TRACE_OPTSOL */
      // #endif /* not FM_CHECKNLP2 */

      if (z < problem_ -> getCutOff ()) {

	problem_ -> setCutOff (objVal);

	t_chg_bounds *chg_bds = NULL;

	if (objInd >= 0) {
	  chg_bds = new t_chg_bounds [problem_ -> nVars ()];
	  chg_bds [objInd].setUpper (t_chg_bounds::CHANGED); 
	}

	// If, with new cutoff, bound tightening clears the whole
	// feasible set, stop
	bool is_still_feas = problem_ -> btCore (chg_bds);

	if (chg_bds)
	  delete [] chg_bds;

	// don't tighten MILP if BT says it's infeasible

	if (!is_still_feas)
	  break;

	// Update lb/ub on milp and nlp here
	const CouNumber 
	  *plb = problem_ -> Lb (),
	  *pub = problem_ -> Ub (),
	  *mlb = milp_    -> getColLower (),
	  *mub = milp_    -> getColUpper ();

	for (int i=problem_ -> nVars (), j=0; i--; ++j, ++plb, ++pub) {

	  bool neglect = problem_ -> Var (j) -> Multiplicity () <= 0;

	  if (*plb > *mlb++) milp_ -> setColLower (j, neglect ? 0. : *plb);
	  if (*pub < *mub++) milp_ -> setColUpper (j, neglect ? 0. : *pub);
	}
      }

      break;

    } else {

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: IP solution NOT MINLP feasible\n");

      if (milpCuttingPlane_ == FP_CUT_EXTERNAL || 
	  milpCuttingPlane_ == FP_CUT_POST) {

	// Solution is IP- but not MINLP feasible: it might get cut by
	// linearization cuts. If so, add a round of cuts and repeat.

	OsiCuts cs;

	problem_   -> domain () -> push (milp_);
	couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); // remaining three arguments NULL by default
	problem_   -> domain () -> pop ();

	if (cs.sizeRowCuts ()) { 

	  // the (integer, NLP infeasible) solution could be separated

	  milp_ -> applyCuts (cs);

	  // found linearization cut, now re-solve MILP (not quite a FP)
	  if (milpCuttingPlane_ == FP_CUT_EXTERNAL && 
	      nSep++ < nSepRounds_)
	    continue;
	}
      }
    }

    //
    // reset number of separation during non-separating iteration
    //

    nSep = 0;

    // NONLINEAR PART ///////////////////////////////////////////////////////

    // fix some variables and solve the NLP to find a NLP (possibly
    // non-MIP) feasible solution

    z = solveNLP (iSol, nSol); 

    if ((nSol && iSol) &&
	(problem_ -> Jnlst () -> ProduceOutput (J_WARNING, J_NLPHEURISTIC))) {

      double dist = 0.;
      int nNonint = 0;

      for (int i = 0; i < problem_ -> nVars (); ++i) {

	exprVar *e = problem_ -> Var (i);

	if (e -> Multiplicity () <= 0)
	  continue;

	if (e -> isInteger () &&
	    (fabs (iSol [i] - ceil (iSol [i] - .5)) > COUENNE_EPS))
	  ++nNonint;

	dist += 
	  (iSol [i] - nSol [i]) * 
	  (iSol [i] - nSol [i]);
      }

      printf ("FP: after NLP, distance %g, %d nonintegers\n", sqrt (dist), nNonint);
    }

    if (problem_ -> Jnlst () -> ProduceOutput (J_WARNING, J_NLPHEURISTIC)) {

      printDist   (problem_, iSol, nSol);
      printCmpSol (problem_, iSol, nSol, 0);
    }

    if (z > COIN_DBL_MAX/2) // something went wrong in the NLP, bail out
      break;

    // check if newly found NLP solution is also integer

    isChecked = false;

    if (nSol) {

      isChecked = problem_->checkNLP0 (nSol, z, true,
				       false, // do not care about obj 
				       true, // stopAtFirstViol
				       true); // checkALL
    }

    if (nSol &&	isChecked) {

      problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: NLP solution is MINLP feasible\n");

      retval = 1;
      objVal = z;

#ifdef FM_TRACE_OPTSOL // - if ----------------------------
#ifdef FM_CHECKNLP2
      problem_->getRecordBestSol()->update();
#else
      problem_->getRecordBestSol()->update(nSol, problem_->nVars(), z, problem_->getFeasTol());
#endif
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else                  // - else --------------------------
#ifdef FM_CHECKNLP2
      best = problem_->getRecordBestSol()->getModSol(problem_ -> nVars ());
#else
      best   = nSol;
#endif
      objVal = z;
#endif                 // - endif -------------------------

      if (z < problem_ -> getCutOff ()) {

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
    }

    // update weights for hessian, distance, and objective, for both
    // MILP and NLP parts

    //multHessNLP_ *= save_multHessNLP_; // exponential reduction

  } while ((niter++ < maxIter_) && 
	   (retval == 0));

  problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: out of the FP loop (feas=%d\n)", retval);

  // OUT OF THE LOOP ////////////////////////////////////////////////////////

  // If MINLP solution found,
  //
  // 1) restore original objective 
  // 2) fix integer variables
  // 3) resolve NLP

  if (nlp_)
    nlp_ -> setObjective (originalObjective);

  if (retval > 0) {

    problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: fix integer variables and rerun NLP\n");

    if (!nlp_) // first call (in this run of FP). Create NLP
      nlp_ = new CouenneTNLP (problem_);

    problem_ -> domain () -> push (*(problem_ -> domain () -> current ()));

    // fix integer coordinates of current (MINLP feasible!) solution
    // and, if feasible, set it as initial (obviously NLP feasible)
    // solution

    if (fixIntVariables (best)) {

      nlp_ -> setInitSol (best);

      ////////////////////////////////////////////////////////////////

      //app_ -> Options () -> SetStringValue ("fixed_variable_treatment", "make_parameter");

      // Solve with original objective function
      status = app_ -> OptimizeTNLP (nlp_);

      ////////////////////////////////////////////////////////////////

      if ((status != Solve_Succeeded) && 
	  (status != Solved_To_Acceptable_Level))
 
	problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, 
					"Feasibility Pump: error in final NLP problem (due to fixing integer variables?)\n");

      if (nlp_ -> getSolution () && (problem_ -> Jnlst () -> ProduceOutput (J_ALL, J_NLPHEURISTIC))) { // check if non-NULL
	printf ("######################## NLP solution (post):\n");
	  for (int i=0; i< problem_ -> nVars ();) {
	    printf ("%+e ", nlp_ -> getSolution () [i]);
	    if (!(++i % 15)) printf ("\n");
	  }
      }

      // if found a solution with the last NLP, check & save it

      double z = nlp_ -> getSolValue ();

      // check if newly found NLP solution is also integer (unlikely...)
      bool isChecked = false;

      if (nSol) {

	problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: found nlp solution, check it\n");

	isChecked = problem_ -> checkNLP0 (nSol, z, true,
					   false,
					   true,
					   true);
      }

      if (nSol &&
	  isChecked &&
	  (z < problem_ -> getCutOff ())) {

	problem_ -> Jnlst () -> Printf (J_WARNING, J_NLPHEURISTIC, "FP: feasible solution is improving\n");

#ifdef FM_TRACE_OPTSOL

#ifdef FM_CHECKNLP2
	problem_->getRecordBestSol()->update();
#else
	problem_->getRecordBestSol()->update(nSol, problem_->nVars(), z, problem_->getFeasTol());
#endif
	best = problem_->getRecordBestSol()->getSol();
	objVal = problem_->getRecordBestSol()->getVal();
#else

#ifdef FM_CHECKNLP2
	best = problem_->getRecordBestSol()->getModSol(problem_ -> nVars ());
#else
	best   = nSol;
#endif
	objVal = z;
#endif

	problem_ -> setCutOff (objVal);
      }
    }

    problem_ -> domain () -> pop ();

    if (problem_ -> Jnlst () -> ProduceOutput (J_WARNING, J_NLPHEURISTIC)) {
      printf ("FP: returning MINLP feasible solution:\n");
      printDist (problem_, best, nSol);
    }

    CoinCopyN (best, problem_ -> nVars (), newSolution);
  }

  // post-call fading update of coefficients

  multDistNLP_  = save_mDN * fadeMult_;
  multHessNLP_  = save_mHN * fadeMult_;
  multObjFNLP_  = save_mON * fadeMult_;
  multDistMILP_ = save_mDM * fadeMult_;
  multHessMILP_ = save_mHM * fadeMult_;
  multObjFMILP_ = save_mOM * fadeMult_;

  if (iSol) delete [] iSol;
  if (nSol) delete [] nSol;

  // release bounding box
  problem_ -> domain () -> pop ();

  // deleted at every call from Cbc, since it changes not only in
  // terms of variable bounds but also in terms of linearization cuts
  // added

  delete milp_;
  delete postlp_;
  milp_ = postlp_ = NULL;

  problem_ -> Jnlst () -> Printf 
    (J_WARNING, J_NLPHEURISTIC, "FP: done ===================\n");

  if (!retval)
    problem_ -> Jnlst () -> Printf 
      (J_ERROR, J_NLPHEURISTIC, "FP: No solution found\n");

  if (retval && (-2 == nCalls_)) {
    problem_ -> bonBase () -> options () -> SetNumericValue ("time_limit", 0.001, "couenne.");
    problem_ -> bonBase () -> setDoubleParameter (Bonmin::BabSetupBase::MaxTime, 0.001);
  }

  return retval;
}

// gather data on single solution for later printout
void compDistSingle (CouenneProblem *p,
		     int n,
		     const double *v,
		     double &norm,
		     int &nInfI,
		     int &nInfN,
		     double &infI,
		     double &infN) {

  p -> domain () -> push (n, v, NULL, NULL); // don't care about bounds as not used

  norm = infI = infN = 0.;
  nInfI = nInfN = 0;

  while (n--) {
    norm += (*v * *v);
    ++v;
  }
  v -= n;

  norm = sqrt (norm);

  for (std::vector <exprVar *>::iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end ();
       ++i) {

    if ((*i) -> Multiplicity () <= 0)
      continue;

    CouNumber 
      vval = (**i) ();

    if ((*i) -> isInteger ()) {

      double inf = CoinMax (vval - floor (vval + COUENNE_EPS),
			    ceil (vval - COUENNE_EPS) - vval);

      if (inf > COUENNE_EPS) {
	++nInfI;
	infI += inf;
      }
    }

    if ((*i) -> Type () == AUX) {

      double 
	diff = 0.,
	fval = (*((*i) -> Image ())) ();

      if      ((*i) -> sign () != expression::AUX_GEQ) diff = CoinMax (diff, vval - fval);
      else if ((*i) -> sign () != expression::AUX_LEQ) diff = CoinMax (diff, fval - vval);

      if (diff > COUENNE_EPS) {
	++nInfN;
	infN += diff;
      }
    }
  }

  p -> domain () -> pop ();
}


// print solutions and distances
void printDist (CouenneProblem *p, const double *iSol, double *nSol) {

  int nInfII = -1, nInfNI = -1, nInfIN = -1, nInfNN = -1;

  double 
    dist  = -1., 
    normI = -1., normN = -1., 
    infII = -1., infNI = -1.,
    infIN = -1., infNN = -1.;

  if (iSol) compDistSingle (p, p -> nVars (), iSol, normI, nInfII, nInfNI, infII, infNI);
  if (nSol) compDistSingle (p, p -> nVars (), nSol, normN, nInfIN, nInfNN, infIN, infNN);
  
  if (iSol && nSol) {

    dist = 0.;

    for (int i = p -> nVars (); i--;) {

      if (p -> Var (i) -> Multiplicity () <= 0)
	continue;

      dist += 
	(iSol [i] - nSol [i]) * 
	(iSol [i] - nSol [i]);
    }

    dist = sqrt (dist);
  }

  printf ("FP: ");
  printf ("IP norm i:%e n:%e dist %e #inf i:%4d n:%4d max inf i:%e n:%e ", normI, normN, dist, nInfII, nInfNI, infII, infNI);
  printf ("NL #inf i:%4d n:%4d max inf i:%e n:%e\n", nInfIN, nInfNN, infIN, infNN);
}


#define WRAP 3

void printCmpSol (CouenneProblem *p, const double *iSol, double *nSol, int direction) {

  int n = p -> nVars ();

  printf ("i:%p n:%p\nFP # ", 
	  (void *) iSol, (void *) nSol);

  double 
    distance = 0.,
    diff;

  char c =
    direction < 0 ? '<' : 
    direction > 0 ? '>' : '-';

  for (int i=0; i<n; i++) {

    if (p -> Var (i) -> Multiplicity () <= 0)
      continue;

    if (i && !(i % WRAP))
      printf ("\nFP # ");

    double 
      iS = iSol ? iSol [i] : 12345., // flag values, to be seen if something is wrong
      nS = nSol ? nSol [i] : 54321.;

    printf ("[%4d %+e -%c- %+e (%e)] ", 
	    i, iS, c, nS,
	    (iSol && nSol) ? fabs (iS - nS) : 0.);

    if (iSol && nSol) {
      diff = iS - nS;
      distance += (diff*diff);
    }
  }

  if (iSol && nSol) {

    distance = sqrt (distance);
    printf ("\n### distance: %e\n", distance);
  }
} 
