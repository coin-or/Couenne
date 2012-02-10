/* $Id$
 *
 * Name:    BTPerfIndicator.cpp
 * Author:  Pietro Belotti
 * Purpose: Measures performance of BT in terms of shrunken bounds -- implementation
 *
 * (C) Pietro Belotti, 2011.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneTypes.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "math.h"
#include "CouenneBTPerfIndicator.hpp"

using namespace Couenne;

/// 
void CouenneBTPerfIndicator::update (CouNumber *lb, CouNumber *ub, int depth) const {

  assert (oldLB_ != NULL && 
	  oldUB_ != NULL);

  double
    curWei = 1./ CoinMax (1., (double) depth),
    newWS  = weightSum_ + curWei;

  // Compute improvement in bounds:

  int 
    nFixed  = 0, 
    nShr    = 0, 
    nShrDbl = 0, 
    nPrInf  = 0;

  double ratio = 0.;

  CouNumber *optimum = problem_ -> bestSol ();

  for (int i=problem_ -> nVars (); i--;) {

    CouNumber 
      olb = oldLB_ [i], oub = oldUB_ [i],
      nlb = lb     [i], nub = ub     [i];

    // Check if optimal solution violated

    if (optimum && 
	((optimum [i] < nlb && optimum [i] > olb) ||
	 (optimum [i] > nub && optimum [i] < oub)))
	
      printf (" %30s cuts optimum at x_%d=%e: [%e,%e] --> [%e,%e], diff:%e\n",
	      name_.c_str (),
	      i, optimum [i], 
	      olb, oub,
	      nlb, nub,
	      CoinMax (nlb - optimum [i],
		       optimum [i] - nub));

    // Check if bound has worsened rather than improved

    if ((((olb > -COUENNE_INFINITY / 1e4) && (nlb < olb - COUENNE_EPS)) ||  
	 ((oub <  COUENNE_INFINITY / 1e4) && (nub > oub + COUENNE_EPS))) && 
	((nlb <= nub + 2 - 1e-5) || i != 0)) // check this is not a wiping cut
	
      printf (" %30s makes bound worse (x%d): [%e,%e] --> [%e,%e], diff:%e\n",
	      name_.c_str (),
	      i, 
	      olb, oub,
	      nlb, nub,
	      CoinMax (olb - nlb,
		       nub - oub));

    if (fabs (nlb - nub) <  COUENNE_EPS) {

      if (fabs (olb - oub) >= COUENNE_EPS)     /// case 1: variable got fixed
	++nFixed;

    } else if (nlb >= nub + COUENNE_EPS)

      ++ nPrInf;

    else if (fabs (nlb - nub) <  COUENNE_INFINITY) { // finite bounds

      if (fabs (olb - oub) >= COUENNE_INFINITY) { /// case 2: interval got finite from [-inf,u] or [l,+inf]

	if ((olb < -COUENNE_INFINITY) &&
	    (oub >  COUENNE_INFINITY))

	  nShr += 2;

	else ++nShr;

      } else  {                                    /// case 3: both old and new were finite

	// printf ("{%g}\t", ratio);

	ratio += (log (oub - olb) - 
		  log (nub - nlb));

	// printf ("%g\tlog (%g-%g) - log (%g-%g) = log %g - log %g = %g - %g = %g ===> %g, %g [%g,%g,%g]\n", 
	// 	boundRatio_      * weightSum_,
	// 	oub, olb, nub, nlb, oub - olb, nub - nlb, 
	// 	log (oub - olb),  log (nub - nlb), 
	// 	log (oub - olb) - log (nub - nlb), 

	// 	ratio, boundRatio_, weightSum_, curWei, newWS);

      }
    } else if ((fabs (olb) > COUENNE_INFINITY) &&
	       (fabs (oub) > COUENNE_INFINITY))  /// case 4: interval became [-inf,u] or [l,+inf] from [-inf,+inf]
      ++nShrDbl;
  }
    
  // Update 

  ++nRuns_;

  ratio /= log (2); // so that it is readable in terms of halvings of the bound interval

  boundRatio_      = (boundRatio_      * weightSum_ + curWei * ratio)   / newWS;
  shrunkInf_       = (shrunkInf_       * weightSum_ + curWei * nShr)    / newWS;
  shrunkDoubleInf_ = (shrunkDoubleInf_ * weightSum_ + curWei * nShrDbl) / newWS;

  nProvedInfeas_   += nPrInf;
  nFixed_          += nFixed;

  //nProvedInfeas_   = (nProvedInfeas_   * weightSum_ + curWei * nPrInf)  / newWS;
  //nFixed_          = (nFixed_          * weightSum_ + curWei * nFixed)  / newWS;

  weightSum_ = newWS;

  delete [] oldLB_;
  delete [] oldUB_;

  oldLB_ = NULL;
  oldUB_ = NULL;
}
