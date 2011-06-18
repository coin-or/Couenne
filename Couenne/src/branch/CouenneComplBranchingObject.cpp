/* $Id$
 *
 * Name:    CouenneComplBranchingObject.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneObject.hpp"
#include "CouenneComplBranchingObject.hpp"

using namespace Ipopt;
using namespace Couenne;

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
							  bool doFBBT, bool doConvCuts, int sign):

  CouenneBranchingObject (solver, originalObject, jnlst, c, p, var, way, brpoint, doFBBT, doConvCuts),
  variable2_    (var2),
  sign_         (sign) {

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

  int 
    ind0 = variable_  -> Index (),
    ind1 = variable2_ -> Index ();

  if (!sign_) { // constraint x1 x2 = 0
    solver -> setColLower (index, 0.);
    solver -> setColUpper (index, 0.);
  } else {

    if (sign_ < 0) // it is a x1 x2 <= 0
      if (!way) {  // branch on second orthant first
	solver -> setColUpper (ind0, 0.);
	solver -> setColLower (ind1, 0.);
      } else {     // branch on fourth orthant first
	solver -> setColLower (ind0, 0.);
	solver -> setColUpper (ind1, 0.);
      }
    else           // it is a x1 x2 >= 0
      if (!way) {  // branch on first orthant first
	solver -> setColLower (ind0, 0.);
	solver -> setColLower (ind1, 0.);
      } else {     // branch on third orthant first
	solver -> setColUpper (ind0, 0.);
	solver -> setColUpper (ind1, 0.);
      }
  }

  ////////////////////////////////////////////////////////////////////

  //CouenneSolverInterface *couenneSolver = dynamic_cast <CouenneSolverInterface *> (solver);
  //CouenneProblem *p = couenneSolver -> CutGen () -> Problem ();

  int 
    nvars  = problem_ -> nVars (),
    objind = problem_ -> Obj (0) -> Body () -> Index ();

  //CouNumber &estimate = way ? upEstimate_ : downEstimate_;
  CouNumber estimate = 0.;//way ? upEstimate_ : downEstimate_;

  bool infeasible = false;

  // only bother doing all of the below if the allocated and pushed
  // stuff will be really used
  if ((doFBBT_ && problem_ -> doFBBT ()) ||
      (doConvCuts_ && simulate_ && cutGen_)) { 

    problem_ -> domain () -> push (nvars,
				   solver -> getColSolution (), 
				   solver -> getColLower    (), 
				   solver -> getColUpper    ()); // have to alloc+copy

    t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

    if (!sign_) {
      chg_bds [index].setLower (t_chg_bounds::CHANGED);
      chg_bds [index].setUpper (t_chg_bounds::CHANGED);
    } else {
      chg_bds [ind0].setLower (t_chg_bounds::CHANGED);
      chg_bds [ind0].setUpper (t_chg_bounds::CHANGED);

      chg_bds [ind1].setLower (t_chg_bounds::CHANGED);
      chg_bds [ind1].setUpper (t_chg_bounds::CHANGED);
    }

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
