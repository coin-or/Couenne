/* $Id:$
 *
 * Name:    BranchCore.cpp
 * Authors: Jim Ostrowski
 * Purpose: Branching step with symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Common Public License (CPL)
 */

//#include "CoinHelperFunctions.hpp"

//#include "OsiRowCut.hpp"

//#include "config_couenne.h"

//#include "CouenneSolverInterface.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
//#include "CouenneCutGenerator.hpp"

#include "CouenneProblem.hpp"

using namespace Couenne;

/** \brief Execute the core of the branch --- need to separate because
    of include conflicts with other packages' config_*.h
 */

void CouenneBranchingObject::branchCore (OsiSolverInterface *solver, int indVar, int way, bool integer, double brpt) {

#ifdef COIN_HAS_NTY
  problem_ -> ChangeBounds( solver -> getColLower (),  solver -> getColUpper (), solver -> getNumCols ());
  problem_ -> Compute_Symmetry();
  
  std::vector< int > branch_orbit;
  //  problem_ -> Print_Orbits();
  //printf("branching on var %i \n", indVar);
  branch_orbit = problem_ -> Find_Orbit (indVar);
  /*
  if(branch_orbit.size() >= 2){
    printf("branching on orbit of size %i \n", branch_orbit.size());
    problem_ -> Print_Orbits();
  }
  */
#endif

  if (!way) solver -> setColUpper (indVar, integer ? floor (brpt) : brpt); // down branch
  else {

#ifdef COIN_HAS_NTY
    for (std::vector<int>::iterator it = branch_orbit.begin(); it!=branch_orbit.end(); ++it) 
      solver -> setColLower (*it, integer ? ceil  (brpt) : brpt); // up   branch
#else
    solver -> setColLower (indVar, integer ? ceil (brpt) : brpt); // up branch
#endif

  }
}
