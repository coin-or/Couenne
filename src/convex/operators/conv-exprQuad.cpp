/* $Id$
 *
 * Name:    conv-exprQuad.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: implementation of convexification methods for exprQuad
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"

#include "CouenneCutGenerator.hpp"

#include "CouenneExprAux.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprBQuad.hpp"

using namespace Couenne;

/// Get lower and upper bound of an expression (if any)
void exprQuad::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLBQuad (this);
  ub = new exprUBQuad (this);

  /*printf ("generated quad bounds:\n  ");
  lb -> print (); printf (" [%g]\n  ", (*lb) ());
  ub -> print (); printf (" [%g]\n",   (*ub) ());*/
}


/// Get lower and upper bound of an expression (if any)
void exprQuad::getBounds (CouNumber &lb, CouNumber &ub) {
  expression::getBounds (lb, ub);
}



// generate equality between *this and *w
void exprQuad::generateCuts (expression *w, //const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg, 
			     int wind, CouNumber lb, CouNumber ub) {

  if (((!(cg -> isFirst ())) &&                     // unless a convexification was never created,
       (fabs ((*this) () - (*w) ()) < COUENNE_EPS)) // do we really need a convexification cut?
      || !alphaConvexify (cg -> Problem ()))        // ... or a new alpha-convexification?
    return;

  /*int 
    nrc = cs.sizeRowCuts (),
    ncc = cs.sizeColCuts ();*/

  // generate linear cuts for convex quadratic [upper|lower]-envelope
  // of this expression
  quadCuts (w, cs, cg);

  /*if (cs.sizeRowCuts () > nrc) {
    printf ("------------------ constraint row cuts\n");
    for (int i=nrc; i<cs.sizeRowCuts (); i++) 
      cs.rowCutPtr (i) -> print ();
  }
  if (cs.sizeColCuts () > nrc) {
    printf ("================== constraint col cuts\n");
    for (int i=ncc; i<cs.sizeColCuts (); i++) 
      cs.colCutPtr (i) -> print ();
      }*/
}
