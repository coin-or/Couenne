/* $Id$
 *
 * Name:    CouenneFPphaseMILP.cpp
 * Authors: Pietro Belotti
 * Purpose: MILP part of the loop
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
int CouenneFeasPump::milpPhase (double *nSol, double *iSol, int niter, int *nsuciter) {

  const int depth = (model_ && (model_ -> currentNode ())) ? model_ -> currentNode () -> depth () : 0;

  double time0 = CoinCpuTime();

  // INTEGER PART /////////////////////////////////////////////////////////

  // Solve IP using nSol as the initial point to minimize weighted
  // l-1 distance from. If nSol==NULL, the MILP is created using the
  // original milp's LP solution.
            
  double z = solveMILP (nSol, iSol, niter, nsuciter);

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
}
