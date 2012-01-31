/* $Id$
 *
 * Name:    BranchCore.cpp
 * Authors: Jim Ostrowski
 * Purpose: Branching step with symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneProblem.hpp"

extern int nOrbBr;

using namespace Ipopt;
using namespace Couenne;

/** \brief Execute the core of the branch --- need to separate code
    because of include conflicts with other packages' config_*.h
 */

void CouenneBranchingObject::branchCore (OsiSolverInterface *solver, int indVar, int way, bool integer, double brpt,
					 t_chg_bounds *&chg_bds) {

  /// only perform orbital branching if
  ///
  /// 1) Nauty has been made available through configure
  /// 2) The orbital_branching option has been set to yes

  if ((doFBBT_ && problem_ -> doFBBT ()) ||
      (doConvCuts_ && simulate_ && cutGen_))
    chg_bds = new t_chg_bounds [problem_ -> nVars ()];

#ifdef COIN_HAS_NTY

  if (problem_ -> orbitalBranching ()) {

    std::vector< int > *branch_orbit = problem_ -> Find_Orbit (indVar);

    // if(branch_orbit -> size() >= 2){
    //   printf("branching on orbit of size %i \n", branch_orbit -> size());
    //   problem_ -> Print_Orbits();
    // }

    if (!way) {

      // DOWN BRANCH: xi <= brpt

      if (jnlst_ -> ProduceOutput (J_ERROR, J_BRANCHING)) {

	printf ("Branch: x%d <= %g [%g,%g]\n", 
		indVar, 
		integer ? floor (brpt) : brpt,
		solver -> getColLower () [indVar], 
		solver -> getColUpper () [indVar]);

	if (problem_ -> bestSol ()) {

	  if ((solver  -> getColUpper () [indVar] > problem_ -> bestSol () [indVar]) &&
	      (brpt                               < problem_ -> bestSol () [indVar]))

	    printf ("Branching rule EXCLUDES optimal solution\n");
	  else
	    for (int i=0; i<problem_ -> nVars (); i++)

	      if ((solver -> getColLower () [indVar] > problem_ -> bestSol () [indVar] + COUENNE_EPS) ||
		  (solver -> getColUpper () [indVar] < problem_ -> bestSol () [indVar] - COUENNE_EPS))

		{printf ("This node does not include optimal solution\n"); break;}
	}
      }

      // BRANCHING RULE -------------------------------------------------------------------------

      solver -> setColUpper (indVar, integer ? floor (brpt) : brpt); // down branch, x [indVar] <= brpt
      if (chg_bds) chg_bds [indVar].setUpper (t_chg_bounds::CHANGED);

    } else {

      // UP BRANCH: xi >= brpt for all i in symmetry group

      jnlst_ -> Printf (J_ERROR, J_BRANCHING, "Branch Symm (%d vars):", branch_orbit -> size ());

      if (branch_orbit -> size () > 1)
	nOrbBr ++;

      bool 
	brExclude   = false, 
	nodeExclude = false;

      for (std::vector<int>::iterator it = branch_orbit -> begin (); it != branch_orbit -> end (); ++it)  {

	assert (*it < problem_ -> nVars ());

	//if (*it >= problem_ -> nVars ()) 
	//continue;

	if (jnlst_ -> ProduceOutput (J_ERROR, J_BRANCHING)) {
	  printf (" x%d>%g [%g,%g]", 
		  *it, 
		  integer ? ceil (brpt) : brpt,
		  solver -> getColLower () [*it], 
		  solver -> getColUpper () [*it]);

	  if (problem_ -> bestSol () &&
	      (solver  -> getColLower () [*it] < problem_ -> bestSol () [*it]) &&
	      (brpt                            > problem_ -> bestSol () [*it]) && !brExclude)

	    brExclude = true;

	  if (problem_ -> bestSol ()) {

	    for (int i=0; i<problem_ -> nVars (); i++)

	      if (((solver -> getColLower () [indVar] > problem_ -> bestSol () [indVar] + COUENNE_EPS) ||
		   (solver -> getColUpper () [indVar] < problem_ -> bestSol () [indVar] - COUENNE_EPS))) {

		nodeExclude = true;
		break;
	      }
	  }
	}

	// BRANCHING RULE -------------------------------------------------------------------------
	if ((integer ? ceil  (brpt) : brpt) > solver -> getColLower () [*it]) {

	  solver -> setColLower (*it, integer ? ceil  (brpt) : brpt); // up branch, x [indVar] >= brpt
	  if (chg_bds) chg_bds [*it].setLower (t_chg_bounds::CHANGED);
	}
      }

      if (jnlst_ -> ProduceOutput (J_ERROR, J_BRANCHING)) {
	if (brExclude)   printf (" (Branching EXCLUDES optimal solution)");
	if (nodeExclude) printf (" (This node does not contain optimal solution)");
	printf ("\n");
      }
    }

    return;
  }

#endif

  /// plain (non-orbital) branching

  if (!way) {

    if (jnlst_ -> ProduceOutput (J_ERROR, J_BRANCHING)) {

      printf ("Branch: x%d <= %g [%g,%g] (opt %g)\n", 
	      indVar, 
	      integer ? floor (brpt) : brpt,
	      solver -> getColLower () [indVar], 
	      solver -> getColUpper () [indVar],
	      problem_ -> bestSol () ? problem_ -> bestSol () [indVar] : 0.);

      if (problem_ -> bestSol ()) {

	if ((solver  -> getColUpper () [indVar] >= problem_ -> bestSol () [indVar]) &&
	    (brpt                               <  problem_ -> bestSol () [indVar]))

	  printf ("Branching EXCLUDES optimal solution\n");
	else
	  for (int i=0; i<problem_ -> nVars (); i++)

	    if ((solver -> getColLower () [i] > problem_ -> bestSol () [i] + COUENNE_EPS) ||
		(solver -> getColUpper () [i] < problem_ -> bestSol () [i] - COUENNE_EPS))

	      {printf ("This node does not contain optimal solution: x%d in [%g,%g] (%g)\n", 
		       i, solver -> getColLower () [i], solver -> getColUpper () [i], problem_ -> bestSol () [i]); break;}
      }
    }

    // BRANCHING RULE -------------------------------------------------------------------------
    solver -> setColUpper (indVar, integer ? floor (brpt + COUENNE_EPS) : brpt); // down branch, x [indVar] <= brpt
    if (chg_bds) chg_bds [indVar].setUpper (t_chg_bounds::CHANGED);

  } else {

    if (jnlst_ -> ProduceOutput (J_ERROR, J_BRANCHING)) {

      printf ("Branch: x%d >= %g [%g,%g] (opt %g)\n", 
	      indVar, 
	      integer ? ceil (brpt) : brpt,
	      solver -> getColLower () [indVar], 
	      solver -> getColUpper () [indVar],
	      problem_ -> bestSol () ? problem_ -> bestSol () [indVar] : 0.);

      if (problem_ -> bestSol ()) {

	if ((solver  -> getColLower () [indVar] <= problem_ -> bestSol () [indVar]) &&
	    (brpt                               >  problem_ -> bestSol () [indVar]))

	  printf ("Branching EXCLUDES optimal solution\n");

	else
	  for (int i=0; i<problem_ -> nVars (); i++)

	    if ((solver -> getColLower () [indVar] > problem_ -> bestSol () [indVar] + COUENNE_EPS) ||
		(solver -> getColUpper () [indVar] < problem_ -> bestSol () [indVar] - COUENNE_EPS))

	      {printf ("This node does not contain optimal solution: x%d in [%g,%g] (%g)\n", 
		       i, solver -> getColLower () [i], solver -> getColUpper () [i], problem_ -> bestSol () [i]); break;}
      }
    }

    // BRANCHING RULE -------------------------------------------------------------------------
    solver -> setColLower (indVar, integer ? ceil (brpt - COUENNE_EPS) : brpt); // up branch, x [indVar] >= brpt
    if (chg_bds) chg_bds [indVar].setLower (t_chg_bounds::CHANGED);
  }
}
