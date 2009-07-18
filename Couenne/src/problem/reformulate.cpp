/* $Id$
 *
 * Name:    reformulate.cpp
 * Author:  Pietro Belotti
 * Purpose: apply reformulation
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CouenneTypes.hpp"

#include "exprVar.hpp"

#include "CouenneProblem.hpp"
#include "depGraph.hpp"
#include "lqelems.hpp"


/// preprocess problem in order to extract linear relaxations etc.
void CouenneProblem::reformulate (CouenneCutGenerator *cg) {

  double now = CoinCpuTime ();

  if (domain_.current () == NULL) {

    // create room for problem's variables and bounds, if no domain exists
    CouNumber 
      *x  = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber)),
      *lb = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber)),
      *ub = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber));

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

    jnlst_->Printf(Ipopt::J_ERROR, J_PROBLEM,
		   "Couenne: problem infeasible after reformulation\n");
    // fake infeasible bounds for Couenne to bail out
    for (int i = nVars (); i--;)
      Ub (i) = - (Lb (i) = 1.);
  }

  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();

  // give a value to all auxiliary variables. Do it now to be able to
  // recognize complementarity constraints in fillDependence()
  initAuxs ();

  // fill dependence_ structure
  fillDependence (bonBase_, cg);

  // quadratic handling
  fillQuadIndices ();

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_->Printf(Ipopt::J_ERROR, J_PROBLEM,
    "Couenne: reformulation time %.3fs\n", now);

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

  //writeAMPL ("extended-aw.mod", true);
  //writeAMPL ("original.mod", false);
}
