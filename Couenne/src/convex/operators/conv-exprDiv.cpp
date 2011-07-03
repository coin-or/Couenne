/* $Id$
 *
 * Name:    conv-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization and convexification methods for divisions
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprOp.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

// Create standard formulation of this expression
exprAux *exprDiv::standardize (CouenneProblem *p, bool addAux) {

  exprOp::standardize (p);
  return (addAux ? (p -> addAuxiliary (this)) : new exprAux (this, p -> domain ()));
}


// generate convexification cut for constraint w = x/y
void exprDiv::generateCuts (expression *w, //const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  // compute y bounds

  CouNumber yl, yu;
  arglist_ [1] -> getBounds (yl, yu);

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index (),
      wi = w            -> Index ();

  bool cLW,  cRW,  cLY,  cRY = 
       cLW = cRW = cLY = true;

  if (!(cg -> isFirst ()) && chg) {
    cLW = (chg [wi].lower() != t_chg_bounds::UNCHANGED);
    cRW = (chg [wi].upper() != t_chg_bounds::UNCHANGED);
    cLY = (chg [yi].lower() != t_chg_bounds::UNCHANGED);
    cRY = (chg [yi].upper() != t_chg_bounds::UNCHANGED);
  }

  // no convexification for terms x/y where y=0 is internal to the
  // bounding box

  if ((yl < -0.) && 
      (yu >  0.)) return;   

  CouNumber k;

  enum auxSign sign = cg -> Problem () -> Var (wi) -> sign ();

  // special case #1: y is almost constant (nonzero) --> y = k. We
  // only need a single plane w >/</= x/k.
  if ((fabs (yl-yu) < COUENNE_EPS) && 
      ((fabs (k = ((yl+yu) / 2)) > COUENNE_EPS))) {
    if (cLY || cRY)
      cg -> createCut (cs, 0., sign, wi, 1., xi, -1./k);
    return;
  }

  CouNumber wl, wu;
  w -> getBounds (wl, wu);

  if (lbw > wl) wl = lbw;
  if (ubw < wu) wu = ubw;

  // special case #2: w is almost constant (nonzero) --> w = x/y = k. We
  // only need a single plane x >/</= y*k.

  if ((fabs (wl-wu) < COUENNE_EPS) &&
      ((k = fabs (wl+wu) / 2) > COUENNE_EPS) &&
      // extra condition: either y's bounds are both pos or both neg,
      // or this is an equality
      ((sign==expression::AUX_EQ) || (yl > 0.) || (yu < 0.))) { 

    if (cLW || cRW) {
      if (sign==expression::AUX_EQ || (yl > 0.)) cg -> createCut (cs, 0., sign, yi,  k, xi, -1.);
      else                                       cg -> createCut (cs, 0., sign, yi, -k, xi,  1.);
    }

    return;
  }

  CouNumber xl, xu;
  arglist_ [0] -> getBounds (xl, xu);

  if ((fabs (xl-xu) < COUENNE_EPS) &&
      (fabs (yl-yu) < COUENNE_EPS) &&
      (fabs (wl-wu) < COUENNE_EPS))
    return; // not much to do here...

  // same as product, just a change in coordinates

  //CouNumber *x = w -> domain () -> x ();
  CouNumber *x = cg -> Problem () -> X ();

  bool 
    ineqFullOrthantF = (((sign == expression::AUX_LEQ) && (yl >  0.)) || ((sign == expression::AUX_GEQ) && (yu < -0.))),
    ineqFullOrthantB = (((sign == expression::AUX_LEQ) && (yu < -0.)) || ((sign == expression::AUX_GEQ) && (yl >  0.)));

  unifiedProdCuts (cg, cs,
		   wi, x [wi], wl, wu,
		   yi, x [yi], yl, yu,
		   xi, x [xi], 
		   ineqFullOrthantF ? -COIN_DBL_MAX : xl, 
		   ineqFullOrthantB ?  COIN_DBL_MAX : xu, 
		   chg, 
		   ineqFullOrthantF ? expression::AUX_GEQ :
		   ineqFullOrthantB ? expression::AUX_LEQ : 
 	 	                      expression::AUX_EQ);
}
