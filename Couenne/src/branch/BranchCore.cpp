/* $Id:$
 *
 * Name:    BranchCore.cpp
 * Authors: Jim Ostrowski
 * Purpose: Branching step with symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

/** \brief Execute the core of the branch --- need to separate code
    because of include conflicts with other packages' config_*.h
 */

void CouenneBranchingObject::branchCore (OsiSolverInterface *solver, int indVar, int way, bool integer, double brpt) {

  /// only perform orbital branching if
  ///
  /// 1) Nauty has been made available through configure
  /// 2) The orbital_branching option has been set to yes

#ifdef COIN_HAS_NTY
  if (problem_ -> orbitalBranching ()) {

    problem_ -> ChangeBounds (solver -> getColLower (),  
			      solver -> getColUpper (), 
			      solver -> getNumCols  ());

    problem_ -> Compute_Symmetry();
  
    std::vector< int > *branch_orbit;

    //  problem_ -> Print_Orbits();
    //printf("branching on var %i \n", indVar);

    branch_orbit = problem_ -> Find_Orbit (indVar);

    /*
      if(branch_orbit.size() >= 2){
      printf("branching on orbit of size %i \n", branch_orbit.size());
      problem_ -> Print_Orbits();
      }
    */

    if (!way) {

      jnlst_ -> Printf (J_ERROR, J_BRANCHING, 
			"Branch: x%d <= %g [%g,%g]\n", 
			indVar, 
			integer ? floor (brpt) : brpt,
			solver -> getColLower () [indVar], 
			solver -> getColUpper () [indVar]);

      solver -> setColUpper (indVar, integer ? floor (brpt) : brpt); // down branch

    } else {

      jnlst_ -> Printf (J_ERROR, J_BRANCHING, 
			"Branch Symm:");

      for (std::vector<int>::iterator it = branch_orbit -> begin (); it != branch_orbit -> end (); ++it)  {

	jnlst_ -> Printf (J_ERROR, J_BRANCHING, 
			  " x%d >= %g; ", 
			  *it, 
			  integer ? ceil (brpt) : brpt,
			  solver -> getColLower () [*it], 
			  solver -> getColUpper () [*it]);

	solver -> setColLower (*it, integer ? ceil  (brpt) : brpt); // up   branch
      }

      jnlst_ -> Printf (J_ERROR, J_BRANCHING, "\n");
    }

    return;
  }

#endif

  /// plain (non-orbital) branching

  if (!way) {

    jnlst_ -> Printf (J_ERROR, J_BRANCHING, 
		      "Branch: x%d <= %g [%g,%g]\n", 
		      indVar, 
		      integer ? floor (brpt) : brpt,
		      solver -> getColLower () [indVar], 
		      solver -> getColUpper () [indVar]);

    solver -> setColUpper (indVar, integer ? floor (brpt) : brpt); // down branch
  } else {

    jnlst_ -> Printf (J_ERROR, J_BRANCHING, 
		      "Branch: x%d >= %g [%g,%g]\n", 
		      indVar, 
		      integer ? ceil (brpt) : brpt,
		      solver -> getColLower () [indVar], 
		      solver -> getColUpper () [indVar]);
    
    solver -> setColLower (indVar, integer ? ceil (brpt) : brpt); // up branch
  }
}
