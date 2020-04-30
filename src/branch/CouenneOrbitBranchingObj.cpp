/* $Id$
 *
 * Name:    CouenneOrbitBranchingObj.cpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 *          
 * Purpose: Branching object for auxiliary variables
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneOrbitBranchingObj.hpp"

using namespace Couenne;

namespace Couenne {
class CouenneCutGenerator;
}

// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);


/** \brief Constructor. 
 *
 * Get a variable as an argument and set value_ through a call to
 * operator () of that exprAux.
*/

CouenneOrbitBranchingObj::CouenneOrbitBranchingObj (OsiSolverInterface *solver,
						    const OsiObject * originalObject,
						    JnlstPtr jnlst, 
						    CouenneCutGenerator *cutGen,
						    CouenneProblem *problem,
						    expression *var, 
						    int way, 
						    CouNumber brpoint, 
						    bool doFBBT, bool doConvCuts):

  CouenneBranchingObject (solver, originalObject, jnlst, cutGen, problem, 
			  var, way, brpoint, doFBBT, doConvCuts) {


}



/** \brief Execute the actions required to branch, as specified by the
 *	   current state of the branching object, and advance the
 *         object's state.
 *
 *         Returns change in guessed objective on next branch
 */

double CouenneOrbitBranchingObj::branch (OsiSolverInterface * solver) {

  return 0;
}
