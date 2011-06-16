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

//#include "CouenneExprVar.hpp"
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

using namespace Couenne;

/// When the current IP (non-NLP) point is not MINLP feasible, linear
/// cuts are added and the IP is re-solved not more than this times
const int numConsCutPasses = 5;

void printDist (CouenneProblem *p, double *iSol, double *nSol);

// Solve
int CouenneFeasPump::solution (double &objVal, double *newSolution) {

  const int depth = (model_ -> currentNode ()) ? model_ -> currentNode () -> depth () : 0;

  if ((problem_ -> nIntVars () <= 0) ||                   // feas pump on NLP? Not yet...
      (CoinCpuTime () > problem_ -> getMaxCpuTime ()) ||  // don't start if time is out
      ((numberSolvePerLevel_ >= 0) &&                     // stop FP after a certain level
       (CoinDrand48 () > 1. / CoinMax 
       (1., (double) ((depth - numberSolvePerLevel_) * 
		      (depth - numberSolvePerLevel_))))))
    return 0;

  problem_ -> Jnlst () -> Printf 
    (J_ERROR, J_NLPHEURISTIC, "FP: start ===================\n");

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
    niter  = 0,   // iteration counter
    nsuciter = 0,// counter for consecutive successful application of one MILP solving method
    retval = 0,   // 1 if found a better solution
    nSep   = 0,   // If separation short circuit, max # of consecutive separations
    objInd = problem_ -> Obj (0) -> Body () -> Index ();

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

  do {

    // INTEGER PART /////////////////////////////////////////////////////////

    // Solve IP using nSol as the initial point to minimize weighted
    // l-1 distance from. If nSol==NULL, the MILP is created using the
    // original milp's LP solution.

     double z = solveMILP (nSol, iSol, niter, &nsuciter);

    // placeholder for how to use pool

    // if no MILP solution was found, bail out

    if (!iSol || z >= COIN_DBL_MAX/2) 
      break;

    bool isChecked = false;

    // if a solution was found, but it is in the tabu list, then we
    // have two choices:
    //
    // 1) the pool is empty: do a round of cuts and repeat;
    //
    // 2) the pool is nonempty: extract the best solution from the
    //    pool and use it instead

    CouenneFPsolution checkedSol (problem_, iSol, false); // false is for not allocating space for this

    if (tabuPool_. find (checkedSol) != tabuPool_ . end ()) {

      if (pool_ -> Queue (). empty ()) {

	OsiCuts cs;

	problem_   -> domain () -> push (milp_);
	couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); // remaining three arguments NULL by default
	problem_   -> domain () -> pop ();

	if (cs.sizeRowCuts ()) {

	  milp_ -> applyCuts (cs);

	  if (nSep++ < numConsCutPasses)
	    continue;

	} else break; // nothing left to do, just bail out

      } else {

         // try to find non-tabu solution in the solution pool
         do
         {
            // retrieve the top solution from the pool
            const CouenneFPsolution &newSol = pool_ -> Queue (). top ();
            CoinCopyN (newSol. x (), problem_ -> nVars (), iSol);
            pool_ -> Queue (). pop ();

            // we found a solution that is not in the tabu list
            if (tabuPool_. find(newSol)  == tabuPool_ . end ())
               break;

            // the pool is empty -> bail out
            if (pool_ -> Queue ().empty())
            {
               delete[] iSol;
               iSol = NULL;
            }

         } while( !pool_ -> Queue ().empty() );
      }

    } else tabuPool_. insert (CouenneFPsolution (problem_, iSol));

#ifdef FM_CHECKNLP2
    isChecked = problem_ -> checkNLP2 (iSol, 0, false, // do not care about obj 
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

      // solution is MINLP feasible! Save it.

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

      if (z < problem_ -> getCutOff ()) {

	problem_ -> setCutOff (objVal);

	t_chg_bounds *chg_bds = NULL;

	if (objInd >= 0) {
	
	  chg_bds = new t_chg_bounds [problem_ -> nVars ()];
	  chg_bds [objInd].setUpper (t_chg_bounds::CHANGED); 
	}

	// if bound tightening clears the whole feasible set, stop
	bool is_still_feas = problem_ -> btCore (chg_bds);

	// don't tighten MILP if BT says it's infeasible

	if (chg_bds) 
	  delete [] chg_bds;

	if (!is_still_feas)
	  break;

	// Update lb/ub on milp and nlp here
	const CouNumber 
	  *plb = problem_ -> Lb (),
	  *pub = problem_ -> Ub (),
	  *mlb = milp_    -> getColLower (),
	  *mub = milp_    -> getColUpper ();

	for (int i=problem_ -> nVars (), j=0; i--; ++j, ++plb, ++pub) {
	    
	  if (*plb > *mlb++) milp_ -> setColLower (j, *plb);
	  if (*pub < *mub++) milp_ -> setColUpper (j, *pub);
	}
      }

      break;

    } else {

      // solution is not MINLP feasible, it might get cut by
      // linearization cuts. If so, add a round of cuts and repeat.

      OsiCuts cs;

      problem_   -> domain () -> push (milp_);
      couenneCG_ -> genRowCuts (*milp_, cs, 0, NULL); // remaining three arguments NULL by default
      problem_   -> domain () -> pop ();

      if (cs.sizeRowCuts ()) { 

	// the (integer, NLP infeasible) solution could be separated

	milp_ -> applyCuts (cs);

	// found linearization cut, now re-solve MILP (not quite a FP)
	if (milpCuttingPlane_ && (nSep++ < numConsCutPasses))
	  continue;
      }
    }

    // reset number of separation during non-separating iteration
    nSep = 0;

    // NONLINEAR PART ///////////////////////////////////////////////////////

    // fix some variables and solve the NLP to find a NLP (possibly
    // non-MIP) feasible solution

    z = solveNLP (iSol, nSol); 

    if (problem_ -> Jnlst () -> ProduceOutput (J_ERROR, J_NLPHEURISTIC)) {

      printDist   (problem_,             iSol, nSol);
      printCmpSol (problem_ -> nVars (), iSol, nSol, 0);
    }

    if (z > COIN_DBL_MAX/2) // something went wrong in the NLP, bail out
      break;

    // check if newly found NLP solution is also integer

    isChecked = false;

    if (nSol) {

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
    }

    if (nSol &&	isChecked) {
      
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
      problem_->getRecordBestSol()->update(nSol, problem_->nVars(),
					   z, problem_->getFeasTol());
      best = problem_->getRecordBestSol()->getSol();
      objVal = problem_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
      best   = nSol;
      objVal = z;
#endif /* not FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */
      
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

  } while ((niter++ < maxIter_) && 
	   (retval == 0));

  // OUT OF THE LOOP ////////////////////////////////////////////////////////

  // If MINLP solution found,
  //
  // 1) restore original objective 
  // 2) fix integer variables
  // 3) resolve NLP

  if (retval > 0) {

    if (!nlp_) // first call (in this run of FP). Create NLP
      nlp_ = new CouenneTNLP (problem_);

    problem_ -> domain () -> push (problem_ -> nVars (),
				   problem_ -> domain () -> x  (),
				   problem_ -> domain () -> lb (),
				   problem_ -> domain () -> ub ());

    fixIntVariables (best);
    nlp_ -> setInitSol (best);

    ////////////////////////////////////////////////////////////////

    // Solve with original objective function
    ApplicationReturnStatus status = app_ -> OptimizeTNLP (nlp_);

    ////////////////////////////////////////////////////////////////

    problem_ -> domain () -> pop ();

    if ((status != Solve_Succeeded) && 
	(status != Solved_To_Acceptable_Level))
 
      problem_ -> Jnlst () -> Printf (J_ERROR, J_NLPHEURISTIC, 
				      "Feasibility Pump: error in final NLP problem\n");

    // if found a solution with the last NLP, check & save it

    double z = nlp_ -> getSolValue ();

    // check if newly found NLP solution is also integer (unlikely...)
    bool isChecked = false;

    if (nSol) {
      
#ifdef FM_CHECKNLP2
      isChecked = problem_->checkNLP2(nSol, 0, false, // do not care about obj 
				      true, // stopAtFirstViol
				      true, // checkALL
				      problem_->getFeasTol());
      if (isChecked) {
	z = problem_->getRecordBestSol()->getModSolVal();
      }

#else /* not FM_CHECKNLP2 */
    isChecked = problem_ -> checkNLP (nSol, z, true);
#endif  /* not FM_CHECKNLP2 */
    }

    if (nSol &&
	isChecked &&
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

  if (retval > 0) {

    if (problem_ -> Jnlst () -> ProduceOutput (J_ERROR, J_NLPHEURISTIC)) {

      printf ("FP: returning MINLP feasible solution:\n");
      printDist (problem_, best, nSol);
    }

    CoinCopyN (best, problem_ -> nVars (), newSolution);
  }

  if (iSol) delete [] iSol;
  if (nSol) delete [] nSol;

  // release bounding box
  problem_ -> domain () -> pop ();

  // deleted at every call from Cbc, since it changes not only in
  // terms of variable bounds but also in of linearization cuts added

  delete milp_;
  milp_ = NULL;

  problem_ -> Jnlst () -> Printf 
    (J_ERROR, J_NLPHEURISTIC, "FP: done ===================\n");

  return retval;
}

void compDistSingle (CouenneProblem *p,
		     int n, double *v, 
		     double &norm, 
		     int &nInfI, 
		     int &nInfN, 
		     double &infI, 
		     double &infN) {

  p -> domain () -> push (n, v, NULL, NULL);

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

    CouNumber 
      vval = (**i) ();

    if ((*i) -> Multiplicity () <= 0)
      continue;

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

void printDist (CouenneProblem *p, double *iSol, double *nSol) {

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

    for (int i = p -> nVars (); i--;)
      dist += 
	(iSol [i] - nSol [i]) * 
	(iSol [i] - nSol [i]);

    dist = sqrt (dist);
  }

  printf ("FP: ");

  printf ("IP norm i:%e n:%e dist %e #inf i:%4d n:%4d max inf i:%e n:%e ", 
	  normI, normN, dist, nInfII, nInfNI, infII, infNI);

  printf ("NL #inf i:%4d n:%4d max inf i:%e n:%e\n", 
	  nInfIN, nInfNN, infIN, infNN);
}


#define WRAP 3

void printCmpSol (int n, double *iSol, double *nSol, int direction) {
  
  printf ("i:%p n:%p\nFP # ", 
	  (void *) iSol, (void *) nSol);

  double 
    distance = 0.,
    diff;

  char c =
    direction < 0 ? '<' : 
    direction > 0 ? '>' : '-';

  for (int i=0; i<n; i++) {

    if (i && !(i % WRAP))
      printf ("\nFP # ");

    double 
      iS = iSol ? iSol [i] : 12345.,
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
