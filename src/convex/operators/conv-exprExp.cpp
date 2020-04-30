/* $Id$
 *
 * Name:    conv-exprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the exponential operator
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExprExp.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprPow.hpp"

#include "CouenneProblem.hpp"

using namespace Couenne;

// generate convexification cut for constraint w = this

void exprExp::generateCuts (expression *aux, //const OsiSolverInterface &si, 
			    OsiCuts &cs,  const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  CouNumber l, u;
  argument_ -> getBounds (l, u);

  int w_ind = aux       -> Index (),
      x_ind = argument_ -> Index ();

  bool cL = !chg || (chg [x_ind].lower() != t_chg_bounds::UNCHANGED) || cg -> isFirst ();
  bool cR = !chg || (chg [x_ind].upper() != t_chg_bounds::UNCHANGED) || cg -> isFirst ();

  enum auxSign sign = cg -> Problem () -> Var (w_ind) -> sign ();

  if (fabs (u-l) < COUENNE_EPS) {  // bounds very close, convexify with a single line

    CouNumber x0 = 0.5 * (u+l), ex0 = exp (x0);
    if (cL || cR)
      cg -> createCut (cs, ex0 * (1 - x0), sign, w_ind, 1., x_ind, - ex0);

    return;
  }

  CouNumber x = (cg -> isFirst ()) ? 
                 0 : powNewton ((*argument_) (), (*aux) (), exp, exp, exp);

  // upper segment

  if ((sign != expression::AUX_GEQ)
      && (cL || cR) 
      && (u < log (COUENNE_INFINITY) ) 
      && (l > -    COUENNE_INFINITY / 1e4)) { // tame lower bound

    CouNumber expl     = exp (l),
              oppslope = (expl - exp (u)) / (u - l);

    cg -> createCut (cs, expl + oppslope*l, -1, 
		     w_ind, 1., 
		     x_ind, oppslope);
  }

  // no need to continue if this is an expression of the form y<=e^x
  // (the upper segment is needed only)
  if (sign == expression::AUX_LEQ)
    return;

  // add tangent points: first choose sampling points

  const CouNumber logMC = log (COU_MAX_COEFF);

  // add tangents with finite coefficients
  if (l < - logMC) l = - logMC;
  if (u >   logMC) u =   logMC;

  // approximate the exponential function from below
  cg -> addEnvelope (cs, +1, exp, exp, w_ind, x_ind, x, l, u, chg, true);
}
