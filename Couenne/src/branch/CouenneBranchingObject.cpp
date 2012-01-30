/* $Id$
 *
 * Name:    CouenneBranchingObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

#include "CouenneCutGenerator.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

using namespace Ipopt;
using namespace Couenne;

// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);

/** \brief Constructor. 
 *
 * Get a variable as an argument and set value_ through a call to
 * operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (OsiSolverInterface *solver,
						const OsiObject * originalObject,
						JnlstPtr jnlst, 
						CouenneCutGenerator *cutGen,
						CouenneProblem *problem,
						expression *var, 
						int way, 
						CouNumber brpoint, 
						bool doFBBT, bool doConvCuts):

  OsiTwoWayBranchingObject (solver, originalObject, way, brpoint),

  cutGen_       (cutGen),
  problem_      (problem),
  variable_     (var),
  jnlst_        (jnlst),
  doFBBT_       (doFBBT),
  doConvCuts_   (doConvCuts),
  downEstimate_ (0.),
  upEstimate_   (0.),
  simulate_     (false) {

  // This two-way branching rule is only applied when both lower and
  // upper bound are finite. Otherwise, a CouenneThreeWayBranchObj is
  // used (see CouenneThreeWayBranchObj.hpp).
  //
  // The rule is as follows:
  //
  // - if x is well inside the interval (both bounds are infinite or
  // there is a difference of at least COU_NEAR_BOUND), set
  // value_ to x;
  //
  // - otherwise, try to get far from bounds by setting value_ to a
  // convex combination of current and midpoint
  //
  // TODO: consider branching value that maximizes distance from
  // current point (how?)

  CouNumber lb, ub;

  variable_ -> getBounds (lb, ub);

  // bounds may have tightened and may exclude value_ now, update it

  value_ = (fabs (brpoint) < large_bound) ? brpoint : (*variable_) ();

  if   (lb < -large_bound)
    if (ub >  large_bound) value_ = 0.;                                                   // ]-inf,+inf[
    else                   value_ = ((value_ < -COUENNE_EPS) ? (AGGR_MUL * (-1+value_)) : // ]-inf,u]
				     (value_ >  COUENNE_EPS) ? 0. : -AGGR_MUL);
  else
    if (ub >  large_bound) value_ = ((value_ >  COUENNE_EPS) ? (AGGR_MUL *  (1+value_)) : // [l,+inf[
				     (value_ < -COUENNE_EPS) ? 0. :  AGGR_MUL);
    else {                                                                                // [l,u]
      double margin = fabs (ub-lb) * closeToBounds;
      if      (value_ < lb + margin) value_ = lb + margin;
      else if (value_ > ub - margin) value_ = ub - margin;
    }

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		    "Branch: x%-3d will branch on %g (cur. %g) [%g,%g]; firstBranch_ = %d\n", 
		    variable_ -> Index (), value_, (*variable_) (), lb, ub, firstBranch_);
}


/** \brief Execute the actions required to branch, as specified by the
 *	   current state of the branching object, and advance the
 *         object's state.
 *
 *         Returns change in guessed objective on next branch
 */

double CouenneBranchingObject::branch (OsiSolverInterface * solver) {

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		    "::branch: at %.20e\n", value_);

  // way = 0 if "<=" node, 
  //       1 if ">=" node

  int 
    way   = (!branchIndex_) ? firstBranch_ : !firstBranch_,
    index = variable_ -> Index ();

  bool 
    integer    = variable_ -> isInteger (),
    infeasible = false;

  CouNumber
    l      = solver -> getColLower () [index],
    u      = solver -> getColUpper () [index],
    brpt   = value_;

  if      (brpt < l) brpt = l;
  else if (brpt > u) brpt = u;

  // anticipate test on bounds w.r.t. test on integrality.

  if   (l <  -large_bound) {
    if (u <=  large_bound) // ]-inf,u]
      brpt = ((u < -COUENNE_EPS) ? CoinMax (CoinMax (brpt, .5 * (l+u)), AGGR_MUL * (-1. + brpt)) : 
	      (u >  COUENNE_EPS) ? 0.                                                            : -AGGR_MUL);
    else brpt = 0.;
  } else
    if (u >  large_bound) // [l,+inf[
      brpt = ((l >  COUENNE_EPS) ? CoinMin (CoinMin (brpt, .5 * (l+u)), AGGR_MUL * ( 1. + brpt)) : 
	      (l < -COUENNE_EPS) ? 0.                                                            :  AGGR_MUL);
    else {                // [l,u] (finite)

      CouNumber point = default_alpha * brpt + (1. - default_alpha) * (l + u) / 2.;

      if      ((point-l) / (u-l) < closeToBounds) brpt = l + (u-l) * closeToBounds;
      else if ((u-point) / (u-l) < closeToBounds) brpt = u + (l-u) * closeToBounds;
      else                                        brpt = point;
    }

  // If brpt is integer and the variable is constrained to be integer,
  // there will be a valid but weak branching. Modify brpt depending
  // on way and on the bounds on the variable, so that the usual
  // integer branching will be performed.

  if (integer && 
      ::isInteger (brpt)) {

    // Look at all possible cases (l,u are bounds, b is the branching
    // point. l,u,b all integer):
    //
    // 1) l <  b <  u: first branch on b +/- 1 depending on branch
    // direction, right branch on b;
    //
    // 2) l <= b <  u: LEFT branch on b, RIGHT branch on b+1
    // 
    // 3) l <  b <= u: LEFT branch on b-1, RIGHT branch on b

    // assert ((brpt - l > .5) || 
    // 	    (u - brpt > .5));

    if ((brpt - l > .5) &&
	(u - brpt > .5)) {// brpt is integer interior point of [l,u]

      if (!branchIndex_) { // if this is the first branch operation
	if (!way) brpt -= (1. - COUENNE_EPS);
	else      brpt += (1. - COUENNE_EPS);
      }
      else { // adjust brpt so that interval covers integer value
	     // necessary for nvs13.nl 
	if (!way) brpt += COUENNE_EPS;
	else      brpt -= COUENNE_EPS;
      }
    } 
    else {
      if (u - brpt > .5) {
	if (way) {
	  brpt += (1. - COUENNE_EPS);
	}
	else { // adjust brpt so that interval covers integer value
	  brpt += COUENNE_EPS;
	}
      }
      else { 
	if (brpt - l > .5) {
	  if (!way) {
	    brpt -= (1. - COUENNE_EPS);
	  }
	  else { // adjust brpt so that interval covers integer value
	    brpt -= COUENNE_EPS;
	  }
	}
	else { // u == l == brpt; must still branch to fix variables in object
	  // but one branch is infeasible
	  if(way) {
	    solver->setColLower(index, u+1); // infeasible
	  }
	}
      }
    }
  }

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "Branching: x%-3d %c= %.20e\n", 
		    index, way ? '>' : '<', integer ? (way ? ceil (brpt) : floor (brpt)) : brpt);

  if (way) {
    if      (brpt < l)             jnlst_->Printf (J_STRONGWARNING, J_BRANCHING, "Nonsense up-br: [ %.8f ,(%.8f)] -> %.8f\n",l,u,value_);
    else if (brpt < l+COUENNE_EPS) jnlst_->Printf (J_STRONGWARNING, J_BRANCHING, "## WEAK  up-br: [ %.8f ,(%.8f)] -> %.8f\n",l,u,value_);
  } else {
    if      (brpt > u)             jnlst_->Printf (J_STRONGWARNING, J_BRANCHING, "Nonsense dn-br: [(%.8f), %.8f ] -> %.8f\n",l,u,value_);
    else if (brpt > u-COUENNE_EPS) jnlst_->Printf (J_STRONGWARNING, J_BRANCHING, "## WEAK  dn-br: [(%.8f), %.8f ] -> %.8f\n",l,u,value_);
  }

  /*
  double time = CoinCpuTime ();
  jnlst_ -> Printf (J_VECTOR, J_BRANCHING,"[vbctool] %02d:%02d:%02d.%02d_I x%d%c=%g_[%g,%g]\n",
		    (int) (floor(time) / 3600), 
		    (int) (floor(time) / 60) % 60, 
		    (int) floor(time) % 60, 
		    (int) ((time - floor (time)) * 100),
		    index, way ? '>' : '<', integer ? ((way ? ceil (brpt): floor (brpt))) : brpt,
		    solver -> getColLower () [index],
		    solver -> getColUpper () [index]);
  */

  t_chg_bounds *chg_bds = NULL;

  // Apply core branching decision (usually new variable bound)

  branchCore (solver, index, way, integer, brpt, chg_bds);

  //CouenneSolverInterface *couenneSolver = dynamic_cast <CouenneSolverInterface *> (solver);
  //CouenneProblem *p = cutGen_ -> Problem ();
  //couenneSolver -> CutGen () -> Problem ();

  int 
    nvars  = problem_ -> nVars (),
    objind = problem_ -> Obj (0) -> Body () -> Index ();

  //CouNumber &estimate = way ? upEstimate_ : downEstimate_;
  CouNumber estimate = 0.;//way ? upEstimate_ : downEstimate_;

  if ((doFBBT_ && problem_ -> doFBBT ()) ||
      (doConvCuts_ && simulate_ && cutGen_)) {

    problem_ -> domain () -> push (solver); // have to alloc+copy

    if (            doFBBT_ &&           // this branching object should do FBBT, and
	problem_ -> doFBBT ()) {         // problem allowed to do FBBT

      problem_ -> installCutOff ();

      if (!problem_ -> btCore (chg_bds)) // done FBBT and this branch is infeasible
	infeasible = true;               // ==> report it

      else {

	const double
	  *lb = solver -> getColLower (),
	  *ub = solver -> getColUpper ();

	//CouNumber newEst = problem_ -> Lb (objind) - lb [objind];
	estimate = CoinMax (0., objind >= 0 ? (problem_ -> Lb (objind) - lb [objind]) : 0.);

	for (int i=0; i<nvars; i++) {
	  if (problem_ -> Lb (i) > lb [i]) solver -> setColLower (i, problem_ -> Lb (i));
	  if (problem_ -> Ub (i) < ub [i]) solver -> setColUpper (i, problem_ -> Ub (i));
	}
      }
    }

    if (!infeasible && doConvCuts_ && simulate_ && cutGen_) { 

      // generate convexification cuts before solving new node's LP
      int nchanged, *changed = NULL;
      OsiCuts cs;

      // sparsify structure with info on changed bounds and get convexification cuts
      sparse2dense (nvars, chg_bds, changed, nchanged);
      cutGen_ -> genRowCuts (*solver, cs, nchanged, changed, chg_bds);
      free (changed);

      solver -> applyCuts (cs);
    }

    //delete [] chg_bds;

    problem_ -> domain () -> pop ();
  }

  if (chg_bds) delete [] chg_bds;

  // next time do other branching
  branchIndex_++;

  return (infeasible ? COIN_DBL_MAX : estimate); // estimated change in objective function
}
