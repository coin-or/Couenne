/* $Id$
 *
 * File:    addEnvelope.cpp
 * Author:  Pietro Belotti, Carnegie Mellon University
 * Purpose: Add generic envelope to convex function based on function and its derivative
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiRowCut.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneFunTriplets.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

void CouenneCutGenerator::addEnvelope (OsiCuts &cs, int sign,
				       unary_function f,      // function to be linearized
				       unary_function fprime, // derivative of f
				       int w_ind, int x_ind, 
				       CouNumber x, CouNumber l, CouNumber u,
				       t_chg_bounds *chg,
				       bool is_global) const {

  simpletriplet st (f, fprime, fprime); // don't really care about third parameter...
  addEnvelope (cs, sign, &st, w_ind, x_ind, x, l, u, chg, is_global);
}


void CouenneCutGenerator::addEnvelope (OsiCuts &cs, int sign,
				       funtriplet *ft,
				       int w_ind, int x_ind, 
				       CouNumber x, CouNumber l, CouNumber u,
				       t_chg_bounds *chg,
				       bool is_global) const {
  CouNumber opp_slope = - ft -> Fp (x);

  // TODO: remove check of !firstcall_ if point is available already

  // if bounds are very close, convexify with a single line

  bool cLeft  = !chg || (chg [x_ind].lower() != t_chg_bounds::UNCHANGED) || firstcall_;
  bool cRight = !chg || (chg [x_ind].upper() != t_chg_bounds::UNCHANGED) || firstcall_;

  if (fabs (u - l) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (u+l), fp0 = ft -> Fp (x0);
    if (cLeft || cRight) 
      createCut (cs, ft -> F(x0) - fp0 * x0, 0, w_ind, 1., x_ind, - fp0);
    return;
  }

  // Add tangent in any case

  if (((!firstcall_) || ((x >= l) && (x <= u)))
      && !CoinIsnan (opp_slope) 
      && (fabs (opp_slope) < COUENNE_INFINITY)) {

    if (!(problem_ -> Var (x_ind) -> isInteger ()))

      // normal cut
      createCut (cs, ft -> F (x) + opp_slope * x, sign, w_ind, 1., 
		 x_ind, opp_slope, -1, 0., is_global);

    else {

      // If the independent variable is integer, the envelope cut can
      // be made tighter: it is the line through
      //
      // (x1,y1) = (floor (x), ft -> F (floor (x))) and
      // (x2,y2) = (ceil  (x), ft -> F (ceil  (x)))
      //
      // therefore the inequality will be
      //
      // y + ax >=<  b
      //
      // where the sign is determined by the variable sign, and since
      // the line has equation
      //
      // x  - x1     y  - y1
      // -------  =  -------
      // x2 - x1     y2 - y1
      //
      //                  y1 - y2              y1 - y2
      // we have that a = ------- and b = y1 + -------  x1
      //                  x2 - x1              x2 - x1
      //
      // Thanks to Sergey for gently encouraging me to do this :-)

      CouNumber x1, x2, y1, y2;

      if ((x1 = floor (x)) < l)
	x1 = ceil (l - COUENNE_EPS);

      y1 = ft -> F (x1);

      x2 = ceil (x);

      if (fabs (x2-x1) < COUENNE_EPS)
	x2 += 1.;

      y2 = ft -> F (x2);

      if ((x2 > u)   ||
	  CoinIsnan (y1) || CoinIsnan (y2) || 
	  !CoinFinite (y1) || !CoinFinite (y2)) // fall back to non-integer cut

	createCut (cs, ft -> F (x) + opp_slope * x, sign, w_ind, 1., 
		   x_ind, opp_slope, -1, 0., is_global);

      else {

	CouNumber 
	  slope = (y1-y2) / (x2-x1),
	  rhs   = y1 + slope * x1;

	createCut (cs, rhs, sign, w_ind, 1., 
		   x_ind, slope, -1, 0., is_global);
      }

      // TODO: if the DEPENDENT variable is integer also, in principle
      // the cut can be made tighter, but not as simply.

    }
  }

  // If this is the first call, add a set of cuts by dividing the
  // interval in a grid and add a cut at each point of the grid

  if ((convtype_ == UNIFORM_GRID) || firstcall_) {

    if (cLeft || cRight) {
      // now add tangent at each sampling point

      CouNumber 
	sample = l, 
	step   = (u-l) / nSamples_;

      //    printf ("[%.4f %.4f], step = %.4f, %d samples\n", 
      //	    l, u, step, nSamples_);

      for (int i = 0; i <= nSamples_; i++) {

	opp_slope = - ft -> Fp (sample);

	if ((fabs (opp_slope) < COUENNE_INFINITY) && 
	    (fabs (sample-x) > COUENNE_EPS)) { // do not add twice cut at current point


	  // same as above, check integrality of argument

	  if (!(problem_ -> Var (x_ind) -> isInteger ()))

	    createCut (cs, ft -> F (sample) + opp_slope * sample, sign, 
		       w_ind, 1.,
		       x_ind, opp_slope, 
		       -1, 0., is_global);
	  else {

	    CouNumber x1, x2, y1, y2;

	    if ((x1 = floor (sample)) < l)
	      x1 = ceil (l - COUENNE_EPS);

	    y1 = ft -> F (x1);

	    x2 = ceil (sample);

	    if (fabs (x2-x1) < COUENNE_EPS)
	      x2 += 1.;

	    y2 = ft -> F (x2);

	    if ((x2 > u)   ||
		CoinIsnan (y1) || CoinIsnan (y2) || 
		!CoinFinite (y1) || !CoinFinite (y2)) // fall back to non-integer cut

	      createCut (cs, ft -> F (sample) + opp_slope * sample, sign, w_ind, 1., 
			 x_ind, opp_slope, -1, 0., is_global);
	  }
      	}

	//	printf ("  Uniform %d: ", i); cut -> print ();

	sample += step;
      }
    }
  }
  else if (convtype_ != CURRENT_ONLY) {

    CouNumber sample = x;

    if (fabs (opp_slope) < COUENNE_INFINITY)
      createCut (cs, ft -> F (x) + opp_slope * x, sign, 
		 w_ind, 1.,
		 x_ind, opp_slope, 
		 -1, 0.,
		 is_global);
      //      printf ("  Current tangent: "); cut -> print ();

    for (int i = 0; i <= nSamples_/2; i++) {

      sample += (x-l) / nSamples_;
      opp_slope = - ft -> Fp (sample);

      if (fabs (opp_slope) < COUENNE_INFINITY)
	createCut (cs, ft -> F (sample) + opp_slope * sample, sign, 
		   w_ind, 1.,
		   x_ind, opp_slope, 
		   -1, 0.,
		   is_global);
	//	printf ("  neighbour -%d: ", i); cut -> print ();
    }

    sample = x;

    for (int i = 0; i <= nSamples_/2; i++) {

      sample += (u-x) / nSamples_;
      opp_slope = - ft -> Fp (sample);

      createCut (cs, ft -> F(sample) + opp_slope * sample, sign, 
		 w_ind, 1.,
		 x_ind, opp_slope, 
		 -1, 0.,
		 is_global);
	//	printf ("  neighbour  %d: ", i); cut -> print ();
    }
  }
}
