/* $Id$
 *
 * Name:    updateBranchInfo.cpp
 * Author:  Pietro Belotti
 * Purpose: get new bounds from parents' bounds + branching rules
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneTypes.hpp"
#include "CouenneProblem.hpp"
#include "BonAuxInfos.hpp"

namespace Couenne {

  //#define DEBUG

/// Get changed bounds due to branching
void updateBranchInfo (const OsiSolverInterface &si, CouenneProblem *p, 
		       t_chg_bounds *chg_bds, const CglTreeInfo &info) {

  int ncols = p -> nVars ();

  if ((info.inTree) && (info.pass==0)) {

    // We are anywhere in the B&B tree but at the root node. Check,
    // through the auxiliary information, which bounds have changed
    // from the parent node.

    OsiBabSolver *auxinfo = dynamic_cast <OsiBabSolver *> (si.getAuxiliaryInfo ());

    bool 
      have_parent_lower = false,
      have_parent_upper = false;

    if (auxinfo && (auxinfo -> extraCharacteristics () & 2)) {

      // get previous bounds
      const double 
	*beforeLower = auxinfo -> beforeLower (),
	*beforeUpper = auxinfo -> beforeUpper ();

      if (beforeLower || beforeUpper) {

	// get current bounds
	const double 
	  *nowLower = p -> Lb (),
	  *nowUpper = p -> Ub (); //si.getColUpper ();

	if (beforeLower) {

	  have_parent_lower = true;

	  for (int i=0; i < ncols; i++)
	    if (*nowLower++ >= *beforeLower++ + COUENNE_EPS)
	      chg_bds [i].setLower (t_chg_bounds::CHANGED);
	}

	if (beforeUpper) {

	  have_parent_upper = true;

	  for (int i=0; i < ncols; i++)
	    if (*nowUpper++ <= *beforeUpper++ - COUENNE_EPS)
	      chg_bds [i].setUpper (t_chg_bounds::CHANGED);
	}
      } 
    }

    // not all bounds are available

    if (!have_parent_lower || 
	!have_parent_upper) {

#ifdef DEBUG
      printf ("### Warning: could not access parent node's %sbounds in generateCuts()\n",
	      have_parent_lower ? "upper " : have_parent_upper ? "lower " : "");
#endif

      // have to assume ALL bounds have changed
      if (!have_parent_lower) 
	for (int i=0; i < ncols; i++) 
	  chg_bds [i].setLower (t_chg_bounds::CHANGED);

      if (!have_parent_upper) 
	for (int i=0; i < ncols; i++) 
	  chg_bds [i].setUpper (t_chg_bounds::CHANGED);
    }
  }
}
}
