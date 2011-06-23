/* $Id$
 *
 * Name:    CouenneVarObject.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Carnegie-Mellon University, 2008-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "BonBabSetupBase.hpp"

#include "CouenneProblem.hpp"
#include "CouenneVarObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneComplObject.hpp"
#include "CouenneProblemElem.hpp"

using namespace Ipopt;
using namespace Couenne;

/// Constructor with information for branching point selection strategy
CouenneVarObject::CouenneVarObject (CouenneCutGenerator *c,
				    CouenneProblem *p,
				    exprVar *ref, 
				    Bonmin::BabSetupBase *base, 
				    JnlstPtr jnlst,
				    int varSelection):

  // Do not set variable (expression).
  // This way, no expression-dependent strategy is chosen
  CouenneObject (c, p, ref, base, jnlst),
  varSelection_ (varSelection) {

  if (jnlst_ -> ProduceOutput (J_SUMMARY, J_BRANCHING)) {
    printf ("created Variable Object: "); 
    reference_ -> print (); 
    printf (" with %s strategy [clamp=%g, alpha=%g]\n", 
	    (strategy_ == LP_CLAMPED)   ? "lp-clamped" : 
	    (strategy_ == LP_CENTRAL)   ? "lp-central" : 
	    (strategy_ == BALANCED)     ? "balanced"   : 
	    (strategy_ == MIN_AREA)     ? "min-area"   : 
	    (strategy_ == MID_INTERVAL) ? "mid-point"  : 
	    (strategy_ == NO_BRANCH)    ? "no-branching (null infeasibility)" : 
	                                  "no strategy",
	    lp_clamp_, alpha_);
  }
}


/// apply the branching rule 
OsiBranchingObject *CouenneVarObject::createBranch (OsiSolverInterface *si, 
						    const OsiBranchingInformation *info, 
						    int way) const {

  // Before anything, use the LP point if so instructed

  problem_ -> domain () -> push 
    (problem_ -> nVars (),
     info -> solution_,
     info -> lower_,
     info -> upper_, false); // don't have to alloc+copy

  OsiBranchingObject *obj = NULL;

  // For some obscure reason, this only seems to work if
  // strong/reliability branching is not used. I suppose it has to do
  // with omitted pseudocost multiplier update, or some other hidden
  // part of the code. 

  if ((varSelection_ == Bonmin::BabSetupBase::OSI_SIMPLE) && 
      ((strategy_ == CouenneObject::LP_CLAMPED) ||
       (strategy_ == CouenneObject::LP_CENTRAL) ||
       (strategy_ == CouenneObject::MID_INTERVAL))) {

    int indVar = reference_ -> Index ();

    CouNumber
      brpt  = info -> solution_ [indVar],
      l     = info -> lower_    [indVar],
      u     = info -> upper_    [indVar];
      
    // these vanilla assignments would drive Couenne crazy. If any of
    // the bounds is (in absolute value) above 1e20, branching values
    // can be detrimental for bound tightening and convexification
    // cuts, for instance.

    // Actually, even smaller values (such as 1e10) can trigger bad BB
    // tree development (see wall in test files, not globallib/wall)

#define LARGE_VALUE 1e8

    if ((l < - LARGE_VALUE) && 
	(u >   LARGE_VALUE) &&
	(fabs (brpt) > LARGE_VALUE / 10))
      brpt = 0.;

    if (l < - COUENNE_INFINITY) l = -1. - 2 * fabs (brpt);
    if (u >   COUENNE_INFINITY) u = +1. + 2 * fabs (brpt);

    CouNumber width = lp_clamp_ * (u-l);

    switch (strategy_) {

    case CouenneObject::LP_CLAMPED:   brpt = CoinMax (l + width, CoinMin (brpt, u - width));        break;
    case CouenneObject::LP_CENTRAL:   if ((brpt < l + width) || 
					  (brpt > u - width)) brpt = .5 * (l+u);                    break;
    case CouenneObject::MID_INTERVAL: brpt = midInterval (brpt, 
    							  info -> lower_ [indVar], 
    							  info -> upper_ [indVar], info);           break;
    default: assert (false); // this will never be used
    }

    obj = new CouenneBranchingObject (si, this, jnlst_, cutGen_, problem_, reference_, 
							  TWO_LEFT, brpt, doFBBT_, doConvCuts_);
  } else {

    // now deal with the more complicated branching selections

    // The infeasibility on an (independent) variable x_i is given by
    // something more elaborate than |w-f(x)|, that is, a function of
    // all infeasibilities of all expressions which depend on x_i.

    int bestWay;
    const CouenneObject *criticalObject = NULL; // should create the branchingObject

    CouNumber bestPt = computeBranchingPoint (info, bestWay, criticalObject);

    ///////////////////////////////////////////

    int indVar = reference_ -> Index ();

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, ":::: creating branching on x_%d @%g [%g,%g]\n", 
		      indVar, 
		      info -> solution_ [indVar],
		      info -> lower_    [indVar],
		      info -> upper_    [indVar]);

    obj = criticalObject ? 
      criticalObject -> createBranch (si, info, way) :
      new CouenneBranchingObject (si, this, jnlst_, cutGen_, problem_, reference_, 
				  way, bestPt, doFBBT_, doConvCuts_);
  }

  problem_ -> domain () -> pop ();

  return obj;
}


/// compute branching point (used in createBranch ())
CouNumber CouenneVarObject::computeBranchingPoint (const OsiBranchingInformation *info,
						   int& bestWay,
						   const CouenneObject *&criticalObject) const {
  criticalObject = NULL;

  if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
    printf ( "---------- computeBRPT for "); 
    reference_ -> print (); 
    printf (" [%g,%g]\n", 
	    info -> lower_ [reference_ -> Index ()],
	    info -> upper_ [reference_ -> Index ()]);
  }

  expression *brVar = NULL; // branching variable

  CouNumber
    brdistDn = 0.,
    brdistUp = 0.,
    bestPt   = 0.,
   *brPts    = NULL, // branching point(s)
   *brDist   = NULL, // distances from LP point to each new convexification
    maxdist  = - COIN_DBL_MAX;

  bool chosen = false;

  bestWay = TWO_LEFT;

  int 
    whichWay = TWO_LEFT,
    index    = reference_ -> Index ();

  std::set <int> deplist = problem_ -> Dependence () [index];

  for (std::set <int>::iterator i = deplist.begin (); i != deplist.end (); ++i) {

    const CouenneObject *obj = problem_ -> Objects () [*i];

    CouNumber improv = 0.;

    assert (obj -> Reference ());

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      printf ("  dependence: "); 
      obj -> Reference () -> print (); 
      if (reference_ -> Image ()) {printf (" := "); obj -> Reference () -> Image () -> print ();}
      printf ("\n");
    }

    if (obj -> Reference ()) {
      if (obj -> Reference () -> Image ())

	improv = obj -> Reference () -> Image ()
	  -> selectBranch (obj, info,                       // input parameters
			   brVar, brPts, brDist, whichWay); // result: who, where, distances, direction
      else {

	brVar = obj -> Reference ();
	brPts  = (double *) realloc (brPts,      sizeof (double)); 
	brDist = (double *) realloc (brDist, 2 * sizeof (double)); 

	double point = info -> solution_ [obj -> Reference () -> Index ()];

	*brPts = point;
	improv = 0.;

	if (point > floor (point)) {improv =                  brDist [0] = point - floor (point);}
	if (point < ceil  (point)) {improv = CoinMin (improv, brDist [1] = ceil (point) - point);}

	point -= floor (point);
	whichWay = (point < 0.45) ? TWO_LEFT : (point > 0.55) ? TWO_RIGHT : TWO_RAND;
      }
    }

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      printf ("  --> Branching on "); 
      if (brVar) {
	brVar -> print (); 
	if (brPts) 
	  printf (" at %g, improv %g <%g>, indices = %d,%d\n", 
		  *brPts, improv, maxdist, index, brVar -> Index ());
      }
    }

    if (brVar &&
	(brVar -> Index () == index) &&    // it's us!
	(fabs (improv) > maxdist) &&       // this branching seems to induce a higher improvement
	(fabs (*brPts) < COU_MAX_COEFF)) { // and branching pt is limited

      criticalObject = (problem_ -> Objects () [*i]); // set this object as the branch creator

      brdistDn = brDist [0];
      brdistUp = brDist [1];
      chosen = true;
      bestPt = *brPts;
      maxdist = improv;
      bestWay = whichWay;
    }
  }

  // no hits on this VarObject's variable, that is, this variable was
  // never chosen 

  if (!chosen) {

    bestPt = info -> solution_ [index];
    brVar  = reference_;

    CouNumber 
      l     = info -> lower_ [index], 
      u     = info -> upper_ [index],
      width = lp_clamp_ * (u-l);

    switch (strategy_) {

    case CouenneObject::LP_CLAMPED:
      bestPt = CoinMax (l + width, CoinMin (bestPt, u - width));
      break;

    case CouenneObject::LP_CENTRAL: 
      if ((bestPt < l + width) || (bestPt > u - width))
	bestPt = (l+u)/2;
      break;

    case CouenneObject::MID_INTERVAL: 

    default: 
      // all other cases (balanced, min-area)
      bestPt = midInterval (bestPt, l, u, info);

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
	if (CoinMin (fabs (bestPt - l), fabs (bestPt - u)) < 1e-3) {
	  printf ("computed failsafe %g [%g,%g] for ", bestPt, l,u);
	  reference_ -> print (); printf ("\n");
	}
      }

      break;
    }

    if ((l < - large_bound) && 
	(u >   large_bound) &&
	(fabs (bestPt) > large_bound))
      bestPt = 0.;

    brPts  = (double *) realloc (brPts, sizeof (double));
    *brPts = bestPt;

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      printf ("  ::: failsafe:  %g [%g,%g] for ", 
	      bestPt, info -> lower_ [index], info -> upper_ [index]); 
      reference_ -> print ();
      printf ("\n");
    }

  } else {

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      if (CoinMin (fabs (bestPt - info -> lower_ [index]), 
		   fabs (bestPt - info -> upper_ [index])) < 1e-3) {
	printf ("  computed %g [%g,%g] for ", 
		bestPt, info -> lower_ [index], info -> upper_ [index]); 
	reference_ -> print ();
	printf ("\n");
      }
    }
  }

  if (pseudoMultType_ == PROJECTDIST) {

    if (chosen) {
      downEstimate_ = brdistDn;
      upEstimate_   = brdistUp;
    } 
    else downEstimate_ = upEstimate_ = 1.;
  }

  if (brPts)  free (brPts);
  if (brDist) free (brDist);

  return bestPt;
}


#define TOL 0.

/// fix nonlinear coordinates of current integer-nonlinear feasible solution
double CouenneVarObject::feasibleRegion (OsiSolverInterface *solver, 
					 const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  assert (index >= 0);

  double val = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val-TOL);
  solver -> setColUpper (index, val+TOL);

  return 0.;
}


/// are we on the bad or good side of the expression?
bool CouenneVarObject::isCuttable () const {

  const std::set <int>                &deplist = problem_ -> Dependence () [reference_ -> Index ()];
  const std::vector <CouenneObject *> &objects = problem_ -> Objects ();

  for (std::set <int>::const_iterator depvar = deplist. begin ();
       depvar != deplist. end (); ++depvar)
    if (!(objects [*depvar] -> isCuttable ()))
      return false;

  return (!(reference_ -> isInteger ()));
}

