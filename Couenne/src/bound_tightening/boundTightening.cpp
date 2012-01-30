/* $Id$
 *
 * Name:    boundTightening.cpp
 * Author:  Pietro Belotti
 * Purpose: tighten bounds prior to convexification cuts
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneProblemElem.hpp"
#include "BonBabInfos.hpp"
#include "BonCbc.hpp"

// max # bound tightening iterations
#define THRES_IMPROVED 0

using namespace Couenne;

// core of the bound tightening procedure

bool CouenneProblem::btCore (t_chg_bounds *chg_bds) const {

  if (!chg_bds) {

    chg_bds = new t_chg_bounds [nVars ()];

    for (int i=0; i < nVars (); i++) 

    if (Var (i) -> Multiplicity () > 0) {

      chg_bds [i].setLower (t_chg_bounds::CHANGED);
      chg_bds [i].setUpper (t_chg_bounds::CHANGED);
    }
  }

  // Bound propagation and implied bounds ////////////////////

  int   ntightened = 0,
      nbwtightened = 0,
      niter = 0;

  bool first = true;

  installCutOff ();

  // check if bt cuts the optimal solution -- now and after bound tightening
  bool contains_optimum = false;

  if (optimum_ != NULL) {
    contains_optimum = true;
    for (int i=0; i < nVars (); i++)
      if ((optimum_ [i] < Lb (i) * (1 - COUENNE_EPS) - COUENNE_EPS) ||
	  (optimum_ [i] > Ub (i) * (1 + COUENNE_EPS) + COUENNE_EPS)) {
	/*printf ("won't check BT: %d [%g,%g] (%g) -- %g\n", 
		i, Lb (i), Ub (i), optimum_ [i],
		CoinMax (- optimum_ [i] + (Lb (i) * (1 - COUENNE_EPS) - COUENNE_EPS),
		optimum_ [i] - (Ub (i) * (1 + COUENNE_EPS) + COUENNE_EPS)));*/
	contains_optimum = false;
	break;
      }
  }

  if (max_fbbt_iter_)  do {

    if (CoinCpuTime () > maxCpuTime_)
      break;

    // propagate bounds to auxiliary variables

    //    if ((nbwtightened > 0) || (ntightened > 0))
    ntightened = ((nbwtightened > 0) || first) ? 
      tightenBounds (chg_bds) : 0;

    // implied bounds. Call also at the beginning, as some common
    // expression may have non-propagated bounds

    // if last call didn't signal infeasibility
    nbwtightened = ((ntightened > 0) || ((ntightened==0) && first)) ? impliedBounds (chg_bds) : 0;

    if (first)
      first = false;

    if ((ntightened < 0) || (nbwtightened < 0)) {
      Jnlst () -> Printf (Ipopt::J_ITERSUMMARY, J_BOUNDTIGHTENING, "infeasible BT\n");
      return false;
    }

    // continue if EITHER procedures gave (positive) results, as
    // expression structure is not a tree.

    if (contains_optimum) {
      for (int i=0; i<nVars (); i++)
	if ((optimum_[i] < Lb(i) - COUENNE_EPS * (1. + CoinMin (fabs(optimum_ [i]), fabs (Lb(i))))) ||
	    (optimum_[i] > Ub(i) + COUENNE_EPS * (1. + CoinMin (fabs(optimum_ [i]), fabs (Ub(i)))))) {
	  printf ("bound tightening CUTS optimum: x%d [%e,%e] val = %e, violation = %e\n", 
		  i, Lb (i), Ub (i), optimum_ [i],
		  CoinMax (- optimum_ [i] + Lb (i),
		  	     optimum_ [i] - Ub (i)));
	  contains_optimum = false;
	}
    }

    // double width = 0.;

    // int 
    //   ninf  = 0,
    //   ndinf = 0;

    // for (int i=nOrigVars (); i--;)

    //   if ((Lb (i) < -COUENNE_INFINITY/1e3) && 
    // 	  (Ub (i) >  COUENNE_INFINITY/1e3)) ndinf++;
    //   else if ((Lb (i) < -COUENNE_INFINITY/1e3) ||
    // 	       (Ub (i) >  COUENNE_INFINITY/1e3)) ninf++;
    //   else width += Ub (i) - Lb (i);

    // printf ("pass %5d: %5d fwd, %5d bwd, %5d inf, %5d double inf, width: %g\n", 
    // 	    niter, ntightened, nbwtightened, ninf, ndinf, width);

  } while (((ntightened > 0) || (nbwtightened > 0)) && 
	   (ntightened + nbwtightened > THRES_IMPROVED) &&
	   ((max_fbbt_iter_ < 0) || (niter++ < max_fbbt_iter_)));

  fbbtReachedIterLimit_ = ((max_fbbt_iter_ > 0) && (niter >= max_fbbt_iter_));

  // TODO: limit should depend on number of constraints, that is,
  // bound transmission should be documented and the cycle should stop
  // as soon as no new constraint subgraph has benefited from bound
  // transmission. 
  //
  // BT should be more of a graph spanning procedure that moves from
  // one vertex to another when either tightening procedure has given
  // some result. This should save some time especially
  // w.r.t. applying implied bounds to ALL expressions just because
  // one single propagation was found.

  for (int i = 0, j = nVars (); j--; i++) 
    if (Var (i) -> Multiplicity () > 0) {

      // final test 
      if ((Lb (i) > Ub (i) + COUENNE_BOUND_PREC * (1 + CoinMin (fabs (Lb (i)), fabs (Ub (i))))) || 
	  (Ub (i) < - MAX_BOUND) ||
	  (Lb (i) >   MAX_BOUND)) {

	Jnlst () -> Printf (Ipopt::J_ITERSUMMARY, J_BOUNDTIGHTENING, "final test: infeasible BT\n");
	return false;
      }

      // sanity check. Ipopt gives an exception when Lb (i) is above Ub (i)
      if (Lb (i) > Ub (i)) {
	CouNumber swap = Lb (i);
	Lb (i) = Ub (i);
	Ub (i) = swap;
      }
    }

  return true;
}


/// procedure to strengthen variable bounds. Return false if problem
/// turns out to be infeasible with given bounds, true otherwise.

bool CouenneProblem::boundTightening (t_chg_bounds *chg_bds, 
				      Bonmin::BabInfo * babInfo) const {

  //  double startTime = CoinCpuTime ();
  //
  // #define SMALL_BOUND 1e4
  //   for (int i=nOrigVars (); i--;) {
  //     if (Lb (i) < -SMALL_BOUND) Lb (i) = -SMALL_BOUND;
  //     if (Ub (i) >  SMALL_BOUND) Ub (i) =  SMALL_BOUND;
  //   } 

  Jnlst () -> Printf (Ipopt::J_ITERSUMMARY, J_BOUNDTIGHTENING,
		      "Feasibility-based Bound Tightening\n");

  int objInd = Obj (0) -> Body () -> Index ();

  /////////////////////// MIP bound management /////////////////////////////////

  if ((objInd >= 0) && babInfo && (babInfo -> babPtr ())) {

    CouNumber UB      = babInfo  -> babPtr () -> model (). getObjValue(),
              LB      = babInfo  -> babPtr () -> model (). getBestPossibleObjValue (),
              primal0 = Ub (objInd), 
              dual0   = Lb (objInd);

    // Bonmin assumes minimization. Hence, primal (dual) is an UPPER
    // (LOWER) bound.
    
    if ((UB < COUENNE_INFINITY) && 
	(UB < primal0 - COUENNE_EPS)) { // update primal bound (MIP)

      Ub (objInd) = UB;
      chg_bds [objInd].setUpper(t_chg_bounds::CHANGED);
    }

    if ((LB > - COUENNE_INFINITY) && 
	(LB > dual0 + COUENNE_EPS)) { // update dual bound
      Lb (objInd) = LB;
      chg_bds [objInd].setLower(t_chg_bounds::CHANGED);
    }
  }

  return btCore (chg_bds);

  //printf ("Total cpu time = %e\n", CoinCpuTime () - startTime);
  //exit (-1);
  //return retval;
}


/// reduced cost bound tightening
int CouenneProblem::redCostBT (const OsiSolverInterface *psi,
			       t_chg_bounds *chg_bds) const {
  
  int 
    nchanges = 0,
    objind   = Obj (0) -> Body () -> Index ();

  if (objind < 0) 
    return 0;

  CouNumber
    UB = getCutOff (), //babInfo -> babPtr () -> model (). getObjValue(),
    LB = Lb (objind);  //babInfo -> babPtr () -> model (). getBestPossibleObjValue ();

  //////////////////////// Reduced cost bound tightening //////////////////////

  if ((LB > -COUENNE_INFINITY) && 
      (UB <  COUENNE_INFINITY)) {

    const double 
      *X  = psi -> getColSolution (),
      *L  = psi -> getColLower    (),
      *U  = psi -> getColUpper    (),
      *RC = psi -> getReducedCost ();

    if (jnlst_ -> ProduceOutput (Ipopt::J_MATRIX, J_BOUNDTIGHTENING)) {
      printf ("REDUCED COST BT (LB=%g, UB=%g):\n", LB, UB);
      for (int i=0, j=0; i < nVars (); i++) {
	if (Var (i) -> Multiplicity () <= 0)
	  continue;
	printf ("%3d %7e [%7e %7e] c %7e ", i, X [i], L [i], U [i], RC [i]);
	if (!(++j % 3))
	  printf ("\n");
      }
      printf ("-----------\n");
    }

    int ncols = psi -> getNumCols ();

    for (int i=0; i<ncols; i++)
      if ((i != objind) && 
	  (Var (i) -> Multiplicity () > 0)) {

	CouNumber
	  x  = X  [i],
	  l  = L  [i],
	  u  = U  [i],
	  rc = RC [i];

#define CLOSE_TO_BOUNDS 1.e-15

	if ((fabs (rc)  < CLOSE_TO_BOUNDS) || 
	    (fabs (l-u) < CLOSE_TO_BOUNDS)) // no need to check
	  continue;

	bool isInt = Var (i) -> isInteger ();

	if ((rc >= 0.) && 
	    (fabs (x-l) <= CLOSE_TO_BOUNDS)) {

	  if (LB + (u-l)*rc > UB) {

	    CouNumber newUb = l + (UB-LB) / rc; // which is surely < u
	    newUb = !isInt ? newUb : floor (newUb + COUENNE_EPS);

	    if (newUb < Ub (i)) {

	      Ub (i) = newUb;

	      nchanges++;
	      chg_bds [i].setLower (t_chg_bounds::CHANGED);
	    }
	  }

	} else if ((rc <= 0.) && 
		   (fabs (x-u) <= CLOSE_TO_BOUNDS)) {

	  if (LB - (u-l) * rc > UB) {

	    CouNumber newLb = u + (UB-LB) / rc; // recall rc < 0 here

	    newLb = !isInt ? newLb : ceil (newLb - COUENNE_EPS);

	    if (newLb > Lb (i)) {

	      Lb (i) = newLb;

	      nchanges++;
	      chg_bds [i].setUpper (t_chg_bounds::CHANGED);
	    }
	  }
	}
      }

    if (jnlst_ -> ProduceOutput (Ipopt::J_MATRIX, J_BOUNDTIGHTENING)) {
      printf ("AFTER reduced cost bt:\n");
      for (int i=0, j=0; i < nVars (); ++i) {
	if (Var (i) -> Multiplicity () <= 0)
	  continue;
	printf ("%3d [%7e %7e] ", i, Lb (i), Ub (i));
	if (!(++j % 4))
	  printf ("\n");
      }
      printf ("-----------\n");
    }
  }

  return nchanges;
}
