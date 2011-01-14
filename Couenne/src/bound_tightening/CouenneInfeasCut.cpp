/* $Id$
 *
 * Name:    CouenneInfeasCut.cpp
 * Author:  Pietro Belotti
 * Purpose: An infeasible cut to tell the node solver this node is infeasible -- implementation
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiCuts.hpp"

///
/// Add a fictitious cut 1<= x_0 <= -1 as a signal to the node solver
/// that this node is deemed infeasible by this cut generator (most
/// likely a bound tightener).
///

void WipeMakeInfeas (OsiCuts &cs) {

  //for (int i=cs.sizeRowCuts(); i--;) cs. eraseRowCut (i);
  //for (int i=cs.sizeColCuts(); i--;) cs. eraseColCut (i);

  OsiColCut *infeascut = new OsiColCut;

  if (infeascut) {
    int i=0;
    double upper = -1., lower = +1.;
    infeascut -> setLbs (1, &i, &lower);
    infeascut -> setUbs (1, &i, &upper);
    cs.insert (infeascut);
    delete infeascut;
  }
}

///
/// Check whether the previous cut generators have added an infeasible
/// cut.
///

bool isWiped (OsiCuts &cs) {

  if (cs.sizeColCuts () == 0)
  //(cs.sizeColCuts () != 1))
    return false;

  CoinPackedVector 
    lbs = cs.colCutPtr (cs.sizeColCuts () - 1) -> lbs (),
    ubs = cs.colCutPtr (cs.sizeColCuts () - 1) -> ubs ();

  return ((lbs.getNumElements ()  == 1)  &&
	  (ubs.getNumElements ()  == 1)  &&
	  (*(lbs.getIndices   ()) == 0)  &&
	  (*(lbs.getElements  ()) == 1.) &&
	  (*(ubs.getIndices   ()) == 0)  &&
	  (*(ubs.getElements  ()) == -1.));
}
