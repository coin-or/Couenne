/* $Id$
 *
 * Name:    CouenneBranchingObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

//#include "CouenneSolverInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneCutGenerator.hpp"

class CouenneCutGenerator;

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

  value_ = (*variable_) ();

  if (fabs (brpoint) < COUENNE_INFINITY) 
    value_ = brpoint;

  CouNumber lb, ub;
  var -> getBounds (lb, ub);

  // do not branch too close to bounds
  if ((lb > -COUENNE_INFINITY) && (ub < COUENNE_INFINITY)) {
    if      ((value_ - lb) / (ub-lb) < closeToBounds) value_ = lb + (ub-lb) * closeToBounds;
    else if ((ub - value_) / (ub-lb) < closeToBounds) value_ = ub + (lb-ub) * closeToBounds;
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

  if (way) {
    if      (value_ < l)             
      jnlst_->Printf(J_STRONGWARNING, J_BRANCHING, "Nonsense up-br: [ %.8f ,(%.8f)] -> %.8f\n",l,u,value_);
    else if (value_ < l+COUENNE_EPS) 
      jnlst_->Printf(J_STRONGWARNING, J_BRANCHING, "## WEAK  up-br: [ %.8f ,(%.8f)] -> %.8f\n",l,u,value_);
  } else {
    if      (value_ > u)             
      jnlst_->Printf(J_STRONGWARNING, J_BRANCHING, "Nonsense dn-br: [(%.8f), %.8f ] -> %.8f\n",l,u,value_);
    else if (value_ > u+COUENNE_EPS) 
      jnlst_->Printf(J_STRONGWARNING, J_BRANCHING, "## WEAK  dn-br: [(%.8f), %.8f ] -> %.8f\n",l,u,value_);
  }

  if ((brpt < l) || (brpt > u))
    brpt = 0.5 * (l+u);

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "Branching: x%-3d %c= %g\n", 
		    //printf ("Branching: x%-3d %c= %g\n", 
		    index, way ? '>' : '<', integer ? (way ? ceil (brpt) : floor (brpt)) : brpt);

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

  if (!way) solver -> setColUpper (index, integer ? floor (brpt) : brpt); // down branch
  else      solver -> setColLower (index, integer ? ceil  (brpt) : brpt); // up   branch

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

    problem_ -> domain () -> push (nvars,
				   solver -> getColSolution (), 
				   solver -> getColLower    (), 
				   solver -> getColUpper    ()); // have to alloc+copy

    t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

    if (!way) chg_bds [index].setUpper (t_chg_bounds::CHANGED);
    else      chg_bds [index].setLower (t_chg_bounds::CHANGED);

    if (            doFBBT_ &&           // this branching object should do FBBT, and
	problem_ -> doFBBT ()) {         // problem allowed to do FBBT

      problem_ -> installCutOff ();

      if (!problem_ -> btCore (chg_bds)) // done FBBT and this branch is infeasible
	infeasible = true;        // ==> report it
      else {

	const double
	  *lb = solver -> getColLower (),
	  *ub = solver -> getColUpper ();

	//CouNumber newEst = problem_ -> Lb (objind) - lb [objind];
	estimate = CoinMax (0., problem_ -> Lb (objind) - lb [objind]);

	//if (newEst > estimate) 
	//estimate = newEst;

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

    delete [] chg_bds;

    problem_ -> domain () -> pop ();
  }

  // next time do other branching
  branchIndex_++;

  return (infeasible ? COIN_DBL_MAX : estimate); // estimated change in objective function
}
