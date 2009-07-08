/* $Id$
 *
 * Name:    CouenneComplBranchingObject.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneComplBranchingObject.hpp"
#include "CouenneCutGenerator.hpp"

// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);


/** \brief Constructor. 
 *
 * Get a variable as an argument and set value_ through a call to
 * operator () of that exprAux.
*/

CouenneComplBranchingObject::CouenneComplBranchingObject (OsiSolverInterface *solver,
							  const OsiObject * originalObject,
							  JnlstPtr jnlst, 
							  CouenneCutGenerator *c,
							  CouenneProblem *p,
							  expression *var, 
							  expression *var2, 
							  int way, 
							  CouNumber brpoint, 
							  bool doFBBT, bool doConvCuts):

  CouenneBranchingObject (solver, originalObject, jnlst, c, p, var, way, brpoint, doFBBT, doConvCuts),
  variable2_    (var2) {

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		    "Complem. Branch: x%-3d = 0 or x%-3d = 0\n", 
		    way ? variable_ -> Index () : variable2_ -> Index ());
}


/** \brief Execute the actions required to branch, as specified by the
 *	   current state of the branching object, and advance the
 *         object's state.
 *
 *         Returns change in guessed objective on next branch
 */

double CouenneComplBranchingObject::branch (OsiSolverInterface * solver) {

  // way = 0 if "x1=0" node, 
  //       1 if "x2=0" node

  int 
    way   = (!branchIndex_) ? firstBranch_ : !firstBranch_,
    index = way ? variable2_ -> Index () : variable_ -> Index ();

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "Branching: x%-3d = 0\n", 
		    //printf ("complementarity Branching: x%-3d = 0\n", 
		    way ? variable2_ -> Index () : variable_ -> Index ());

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

  ////////////////////////////////////////////////////////////////////

  solver -> setColLower (index, 0.);
  solver -> setColUpper (index, 0.);

  ////////////////////////////////////////////////////////////////////

  //CouenneSolverInterface *couenneSolver = dynamic_cast <CouenneSolverInterface *> (solver);
  //CouenneProblem *p = couenneSolver -> CutGen () -> Problem ();

  int 
    nvars  = problem_ -> nVars (),
    objind = problem_ -> Obj (0) -> Body () -> Index ();

  problem_ -> domain () -> push (nvars,
			  solver -> getColSolution (), 
			  solver -> getColLower    (), 
			  solver -> getColUpper    ()); // have to alloc+copy

  //CouNumber &estimate = way ? upEstimate_ : downEstimate_;
  CouNumber estimate = 0.;//way ? upEstimate_ : downEstimate_;

  t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

  chg_bds [index].setUpper (t_chg_bounds::CHANGED);
  chg_bds [index].setLower (t_chg_bounds::CHANGED);

  bool infeasible = false;

  if (     doFBBT_ &&           // this branching object should do FBBT, and
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

  // next time do other branching
  branchIndex_++;

  return (infeasible ? COIN_DBL_MAX : estimate); // estimated change in objective function
}
