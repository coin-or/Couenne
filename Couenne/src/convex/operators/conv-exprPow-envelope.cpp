/* $Id$
 *
 * Name:    conv-exprPow-envelope.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "CouenneTypes.hpp"
#include "rootQ.hpp"
#include "exprPow.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"
#include "funtriplets.hpp"


// adds convex (upper/lower) envelope to a power function

void addPowEnvelope (const CouenneCutGenerator *cg, OsiCuts &cs,
		     int wi, int xi,
		     CouNumber x, CouNumber y,
		     CouNumber k, 
		     CouNumber l, CouNumber u,
		     int sign) {

  // set x to get a deeper cut (so that we get a tangent which is
  // orthogonal with line through current- and tangent point)

  powertriplet pt (k);

  if (!(cg -> isFirst ()))
    x = powNewton (x, y, &pt);

  if      (x<l) x=l;
  else if (x>u) x=u;

  // limit the bounds for the envelope

  CouNumber powThres = (k<=1) ? COU_MAX_COEFF: pow (COU_MAX_COEFF, 1./k),
            step     = (1 + log (1. + (double) (cg -> nSamples ()))) * powThres / COU_MAX_COEFF;

  // If the bounds are too large, the linearization cuts might have
  // large coefficients. To prevent that, Couenne used to set very
  // small fictitious bounds, resulting in
  //
  // 1) still valid cuts, but
  // 2) a very abrupt change in their coefficients; 
  // 3) cuts that may result in an LP solution far away from x (this
  // behavior recalls that of bundle methods for NDO);
  // 4) a non-exact linearization at the bounds (in theory, necessary
  // for convergence...).
  //
  // New values for l and u, if necessary, are therefore set to the
  // maximum bounds if l and/or u are beyond them.
  //
  // Thanks to Sergey for pointing this out.

  if (l < - powThres + 1) {
    l = - powThres + 1; // keeps bounds reasonably large 
    //l = x - step;
    if (u > powThres - 1) {
      u = powThres - 1;
    //u = x + step;
    }
  } else 
    if (u > powThres - 1) {
      u = powThres - 1;
      //u = x + step;
    }

  // convex envelope
  cg -> addEnvelope (cs, sign, &pt, 
		     wi, xi, x, l, u);
}
