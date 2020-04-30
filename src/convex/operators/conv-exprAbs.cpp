/* $Id$
 *
 * Name:    conv-exprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for |f(x)|
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"

#include "OsiSolverInterface.hpp"
#include "CouenneTypes.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprAbs.hpp"
#include "CouenneExprAux.hpp"

using namespace Couenne;

// generate convexification cut for constraint w = |x|

void exprAbs::generateCuts (expression *w, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {

  int w_ind = w         -> Index (),
      x_ind = argument_ -> Index ();

  CouNumber l, u;
  argument_ -> getBounds (l, u);

  enum auxSign sign = cg -> Problem () -> Var (w_ind) -> sign ();

  bool
    cbase  = !chg || cg -> isFirst (),
    cLeft  = cbase || (chg [x_ind].lower() != t_chg_bounds::UNCHANGED),
    cRight = cbase || (chg [x_ind].upper() != t_chg_bounds::UNCHANGED);

  // if l, u have the same sign, then w = x (l >= 0) or w = -x (u <= 0)

  if      (l >= -0) {if (cLeft)  cg -> createCut (cs, 0., sign, w_ind, 1., x_ind, -1.);}
  else if (u <=  0) {if (cRight) cg -> createCut (cs, 0., sign, w_ind, 1., x_ind, +1.);}
  else {

    // add two global cuts: w >= x and w >= -x
    if (cg -> isFirst () && sign != expression::AUX_LEQ) {
      cg -> createCut (cs, 0., +1, w_ind, 1., x_ind, -1.);
      cg -> createCut (cs, 0., +1, w_ind, 1., x_ind,  1.);
    }

    // otherwise check if at most one of the bounds is infinite: even
    // so, we can still add a plane, whose slope is 1 (x unbounded
    // from above) or -1 (from below)

    if (sign != expression::AUX_GEQ) {

      if (l > - COUENNE_INFINITY) {
	if (u < COUENNE_INFINITY) { // the upper approximation has slope other than -1, 1

	  CouNumber slope = (u+l) / (u-l); // should be stable, l < 0 < u

	  // add an upper segment, which depends on the lower/upper bounds
	  if (cLeft || cRight) 
	    cg -> createCut (cs, -l*(slope+1.), -1, w_ind, 1., x_ind, -slope);
	}
	else // slope = 1
	  if (cLeft) cg -> createCut (cs, -2*l, -1, w_ind, 1., x_ind, -1.);
      }
      else if (u < COUENNE_INFINITY) // slope = -1
	if (cRight) cg -> createCut (cs, 2*u, -1, w_ind, 1., x_ind, 1.);
    }
  }
}
