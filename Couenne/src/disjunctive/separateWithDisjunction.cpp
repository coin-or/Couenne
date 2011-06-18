/* $Id$
 *
 * Name:    separateWithDisjunction.cpp
 * Author:  Pietro Belotti
 * Purpose: generate cuts of disjunction
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneDisjCuts.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

using namespace Ipopt;
using namespace Couenne;

/// generate row cuts for one side of the  disjunction
int CouenneDisjCuts::separateWithDisjunction (OsiCuts *cuts, 
					      OsiSolverInterface &si, 
					      OsiCuts &cs, 
					      const CglTreeInfo &info) const {

  if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS) && 
      ((cuts -> sizeRowCuts ()) ||
       (cuts -> sizeColCuts ()))) {

    printf ("applying unilateral cuts:\n");

    if (cuts -> sizeRowCuts ()) {
      printf ("Row\n");
      for (int i=0; i < cuts -> sizeRowCuts (); i++) cuts -> rowCutPtr (i) -> print ();
    }

    if (cuts -> sizeColCuts ()) {
      printf (" Col\n");
      for (int i=0; i < cuts -> sizeColCuts (); i++) cuts -> colCutPtr (i) -> print ();
    }
  }

  int ncols = si.getNumCols ();
  t_chg_bounds *chg = new t_chg_bounds [ncols]; // all init'd automatically to UNCHANGED
  CouenneProblem *p = couenneCG_ -> Problem ();

  p -> domain () -> push (&si);

  // apply cuts
  for (int i = cuts -> sizeColCuts (); i--;) {

    const CoinPackedVector
      &lb = cuts -> colCutPtr (i) -> lbs (),
      &ub = cuts -> colCutPtr (i) -> ubs ();

    const int 
      *lind = lb.getIndices (),
      *uind = ub.getIndices ();
 
    const double
      *lval = lb.getElements (),    *oLB = si.getColLower (),
      *uval = ub.getElements (),    *oUB = si.getColUpper ();

    for (int j=lb.getNumElements (); j--; lind++, lval++)
      if (*lval > oLB [*lind] + COUENNE_EPS) {
	p -> Lb (*lind) = *lval;
	chg [*lind].setLower (t_chg_bounds::CHANGED);
      }

    for (int j=ub.getNumElements (); j--; uind++, uval++)
      if (*uval < oUB [*uind] - COUENNE_EPS) {
	p -> Ub (*uind) = *uval;
	chg [*uind].setUpper (t_chg_bounds::CHANGED);
      }
  }

  int *changed = NULL, // will be created within sparse2dense
    nchanged = 0;

  sparse2dense (ncols, chg, changed, nchanged);

  couenneCG_ -> genRowCuts (si, *cuts, nchanged, changed, chg);

  p -> domain () -> pop ();

  delete [] chg;

  if (changed)
    free (changed);

  return COUENNE_FEASIBLE;
}
