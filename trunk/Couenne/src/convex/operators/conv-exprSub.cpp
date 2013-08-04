/* $Id$
 *
 * Name:    exprSub.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for the Subtraction class
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

// generate equality between *this and *w
void exprSub::generateCuts (expression *w, //const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, 
			    int wind, CouNumber lb, CouNumber ub) {

  if (!(cg -> isFirst ()))
    return;

  // only add one cut at the beginning

  expression *x = arglist_ [0];
  expression *y = arglist_ [1];

  int wi = w -> Index ();
  int xi = x -> Index ();
  int yi = y -> Index ();

  // If this aux is fixed, don't write 
  //
  // "- w + ax = -b" but just
  // 
  // "ax = -b+ w0" 
  //
  // with w0 its constant value

  CouNumber vlb, vub;
  w -> getBounds (vlb, vub);
  bool uselessAux = (vub < vlb + COUENNE_EPS); 

  // TODO: generalize to sign!= ::EQ

  if ((wind >= 0) || uselessAux) {
    wi = -1; // do not insert w's index if specified as input
    lb = ub = vlb;
  }
  else lb = ub = 0;

  if (xi < 0) {
    CouNumber x0 = x -> Value ();
    lb -= x0;
    ub -= x0;
  }

  if (yi < 0) {
    CouNumber y0 = y -> Value ();
    lb += y0;
    ub += y0;
  }

  enum auxSign sign = cg -> Problem () -> Var (w -> Index ()) -> sign ();

  if      (sign == expression::AUX_GEQ) lb = -COIN_DBL_MAX;
  else if (sign == expression::AUX_LEQ) ub =  COIN_DBL_MAX;

  cg -> createCut (cs, lb, ub, wi, -1., xi, 1., yi, -1., true);
}
