/* $Id$
 *
 * Name:    conv-exprMul-genCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: method to convexify multiplications
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

/// generate convexification cut for constraint w = x*y

void exprMul::generateCuts (expression *w, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {

  expression *xe = arglist_ [0];
  expression *ye = arglist_ [1];

  int wi = w  -> Index (), 
      xi = xe -> Index (), 
      yi = ye -> Index ();

  // if expression is x*c or c*y, with c constant from the problem
  // definition or from the branching rules, the expression is
  // linear. Add one convexification equality constraint.

  // check if either operand is constant

  bool is0const = (xe -> Type () == CONST),
       is1const = (ye -> Type () == CONST);

  CouNumber c0=0, c1=0;

  // is one of the two constant?

  if (is0const) c0 = xe -> Value ();  
  if (is1const) c1 = ye -> Value ();  

  // compute bounds

  CouNumber xl, xu, yl, yu, wl, wu;

  xe -> getBounds (xl, xu);
  ye -> getBounds (yl, yu);
  w  -> getBounds (wl, wu);

  if (lbw > wl) wl = lbw;
  if (ubw < wu) wu = ubw;

  // check if either operator got constant because of the branching
  // rules: 

  bool i0s, i1s = i0s = false;

  // TODO: Fix this!

  // x...

  if (!is0const && ((xu-xl) < COUENNE_EPS)) {

    if (is1const) i0s = (fabs (c1)                 * (xu-xl) < COUENNE_EPS);
    else          i0s = ((fabs (yu) + fabs (yl))   * (xu-xl) < COUENNE_EPS);

    if (i0s) 
      c0 = 0.5 * (xl+xu);
  }

  // ...and y

  if (!is1const && ((yu-yl) < COUENNE_EPS)) {

    if (is0const) i1s = (fabs (c0)                 * (yu-yl) < COUENNE_EPS);
    else          i1s = ((fabs (xu) + fabs (xl))   * (yu-yl) < COUENNE_EPS);

    if (i1s) 
      c1 = 0.5 * (yl+yu);
  }

  if (i0s) is0const = true;
  if (i1s) is1const = true;

  // right now c0 and c1 only have a value if the corresponding
  // expression is constant

  if (is0const || is1const) {

    if (cg -> isFirst () ||            // if first call or
	((xe -> Type () != CONST) &&   // neither term is a defined constant
	 (ye -> Type () != CONST))) {  // (=> implied by branching rule)

      if (is0const && is1const)
	// w = c0*c1, which is either because the intervals got very
	// narrow, or because these are indeed constant (which should
	// have been dealt with in simplify(), but who knows...)
	cg -> createCut (cs, c0 * c1, cg -> Problem () -> Var (wi) -> sign (), wi, 1.);

      else {

	CouNumber coe;
	int ind;

	if (is0const) {coe = c0; ind = yi;} // c*y
	else          {coe = c1; ind = xi;} // x*c

	cg -> createCut (cs, 0., cg -> Problem () -> Var (wi) -> sign (), wi, 1., ind, -coe);
      }
    }

    return;
  }

  enum auxSign sign = cg -> Problem () -> Var (wi) -> sign ();

  // add different cuts, to cut out current point in bounding box but
  // out of the hyperbola's belly

  CouNumber x0 = (*(arglist_ [0])) (),
            y0 = (*(arglist_ [1])) ();

  unifiedProdCuts (cg, cs, 
		   xi, x0,      xl, xu, 
		   yi, y0,      yl, yu,
		   wi, (*w) (), wl, wu,
		   //sign == expression::AUX_LEQ ? -COIN_DBL_MAX : wl, 
		   //sign == expression::AUX_GEQ ?  COIN_DBL_MAX : wu,
		   chg, sign);
}
