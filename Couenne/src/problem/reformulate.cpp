/* $Id$
 *
 * Name:    reformulate.cpp
 * Author:  Pietro Belotti
 * Purpose: apply reformulation
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "BonBabSetupBase.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExprVar.hpp"

#include "CouenneProblem.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneLQelems.hpp"

#include "CouenneRecordBestSol.hpp"

#define THRESHOLD_OUTPUT_REFORMULATE 1000

using namespace Couenne;

/// preprocess problem in order to extract linear relaxations etc.
void CouenneProblem::reformulate (CouenneCutGenerator *cg) {
  
  double now = CoinCpuTime ();

  if (nVars () > THRESHOLD_OUTPUT_REFORMULATE) {
    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "Reformulating problem: "); 
    fflush (stdout);
  }

  if (domain_.current () == NULL) {

    // create room for problem's variables and bounds, if no domain exists
    CouNumber 
      *x  = (CouNumber *) malloc (nVars() * sizeof (CouNumber)),
      *lb = (CouNumber *) malloc (nVars() * sizeof (CouNumber)),
      *ub = (CouNumber *) malloc (nVars() * sizeof (CouNumber));

    for (int i = nVars(); i--;) {
      x  [i] =  0.;
      lb [i] = -COUENNE_INFINITY;
      ub [i] =  COUENNE_INFINITY;
    }

    domain_.push (nVars (), x, lb, ub);
  }

  // link initial variables to problem's domain
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)
    (*i) -> linkDomain (&domain_);

  if (jnlst_ -> ProduceOutput(Ipopt::J_SUMMARY, J_PROBLEM))
    print (std::cout);

  // save -- for statistic purposes -- number of original
  // constraints. Some of them will be deleted as definition of
  // auxiliary variables.
  nOrigCons_    = constraints_. size ();
  nOrigIntVars_ = nIntVars ();

  jnlst_->Printf (Ipopt::J_ERROR, J_PROBLEM,
		  "Problem size before reformulation: %d variables (%d integer), %d constraints.\n",
		  nOrigVars_, nOrigIntVars_, nOrigCons_);

  // reformulation
  if (!standardize ()) { // problem is infeasible if standardize returns false

    jnlst_->Printf(Ipopt::J_ERROR, J_COUENNE,
		   "Problem infeasible after reformulation\n");
    // fake infeasible bounds for Couenne to bail out
    for (int i = nVars (); i--;)
      Ub (i) = - (Lb (i) = 1.);

    return;
  }

  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();

  // give a value to all auxiliary variables. Do it now to be able to
  // recognize complementarity constraints in fillDependence()
  initAuxs ();

  bool *isInt = new bool[nVars()];
  for(int i=0; i<nVars(); i++) {
    if(variables_[i]->isInteger()) {
      isInt[i] = true;
    }
    else {
      isInt[i] = false;
    }
  }
  recBSol->setInitIsInt(isInt, nVars());
  recBSol->setInitDomLb(domain()->lb(), nVars());
  recBSol->setInitDomUb(domain()->ub(), nVars());
  delete[] isInt;

  CouNumber cutoff;
  // check for initial solution given to Couenne. If feasible, set cutoff

#ifdef FM_CHECKNLP2
  cutoff = X (objectives_ [0] -> Body () -> Index ());
  if(checkNLP2(X(), cutoff, false, // do not care about obj value
	       true, // stop at first viol 
	       false, // checkAll
	       getFeasTol())) {
    
    jnlst_ -> Printf (Ipopt::J_ERROR, J_PROBLEM,
		      "Couenne: initial solution (value %g) is MINLP feasible\n",
		      cutoff);

#ifdef FM_TRACE_OPTSOL
    getRecordBestSol()->update();
    setCutOff(getRecordBestSol()->getVal());
#else /* not FM_TRACE_OPTSOL */

#ifdef FM_UP_BND
    setCutOff(getRecordBestSol()->getModSolVal());
#else
    setCutOff(getRecordBestSol()->getModSolVal(), 
	      getRecordBestSol()->getModSol(nVars()));
#endif
#endif /* not FM_TRACE_OPTSOL */

  }
#else /* not FM_CHECKNLP2 */
  if (checkNLP (X (), cutoff = X (objectives_ [0] -> Body () -> Index ()), true)) {
    jnlst_ -> Printf (Ipopt::J_ERROR, J_PROBLEM,
		      "Couenne: initial solution (value %g) is MINLP feasible\n",
		      cutoff);

#ifdef FM_TRACE_OPTSOL
    getRecordBestSol()->update(X(), nVars(), cutoff, getFeasTol());
    setCutOff(getRecordBestSol()->getVal());
#else /* not FM_TRACE_OPTSOL */

#ifdef FM_UP_BND
    setCutOff (cutoff);
#else
    setCutOff (cutoff, X ());    
#endif
#endif /* not FM_TRACE_OPTSOL */

  }
#endif /* not FM_CHECKNLP2 */
 
  // fill dependence_ structure
  fillDependence (bonBase_, cg);

  // quadratic handling
  fillQuadIndices ();

  // if ((now = (CoinCpuTime () - now)) > 10.)
  //   jnlst_->Printf(Ipopt::J_ERROR, J_PROBLEM,
  //   "reformulation time %.3fs\n", now);

  jnlst_->Printf (Ipopt::J_WARNING, J_PROBLEM, "Initializing auxiliaries\n");

  // give a value to all auxiliary variables
  initAuxs ();

  int nActualVars = nIntVars_ = 0;

  // check how many integer variables we have now (including aux)
  for (int i=0; i<nVars(); i++)
    if (variables_ [i] -> Multiplicity () > 0) {

      nActualVars++;
      if (variables_ [i] -> isDefinedInteger ())
	nIntVars_++;
    }

  jnlst_->Printf(Ipopt::J_ERROR, J_PROBLEM,
		  "Problem size after  reformulation: %d variables (%d integer), %d constraints.\n",
		  nActualVars, nIntVars_, nCons());

  // check if optimal solution is available (for debug purposes)
  readOptimum ();

  if (bonBase_) {

    CouNumber 
      art_cutoff =  COIN_DBL_MAX,
      art_lower  = -COIN_DBL_MAX;

    bonBase_ -> options() -> GetNumericValue ("art_cutoff", art_cutoff, "couenne.");
    bonBase_ -> options() -> GetNumericValue ("art_lower",  art_lower,  "couenne.");

    if (art_cutoff <  1.e50) setCutOff (art_cutoff);
    if (art_lower  > -1.e50) {
      int indobj = objectives_ [0] -> Body () -> Index ();
      if (indobj >= 0)
	domain_.lb (indobj) = art_lower;
    }
  }

  if (jnlst_->ProduceOutput(Ipopt::J_DETAILED, J_PROBLEM)) {
    // We should route that also through the journalist
    print (std::cout);
  }

  createUnusedOriginals ();

  if (nOrigVars_  > THRESHOLD_OUTPUT_REFORMULATE)
    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "%.1f seconds\n", CoinCpuTime () - now); 
  else if (nVars () > 2*THRESHOLD_OUTPUT_REFORMULATE) 
    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "Reformulation: %.1f seconds\n", CoinCpuTime () - now); 

  if (orbitalBranching_) {

    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "Setting up symmetry groups\n"); 
    setupSymmetry ();
  }

  //writeAMPL ("extended-aw.mod", true);
  //writeAMPL ("original.mod", false);
}
