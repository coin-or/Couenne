/* $Id$
 *
 * Name:    CouenneFPpool.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Pool of MILP- (and why not? NLP-) feasible solutions for
 *          restart use in the Feasibility Pump
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneFPpool.hpp"

using namespace Couenne;

/// CouenneProblem-aware constructor
CouenneFPsolution::CouenneFPsolution (CouenneProblem *p, CouNumber *x):
  x_          (CoinCopyOfArray (x, p -> nVars ())),
  n_          (p -> nVars ()),
  nNLPinf_    (0),
  nIinf_      (0),
  objVal_     (0.),
  maxNLinf_   (0.),
  maxIinf_    (0.) {


}


/// independent constructor --- must provide other data as no
/// CouenneProblem to compute them
CouenneFPsolution::CouenneFPsolution (CouNumber *x, 
				      int n, 
				      int nNLPinf,
				      int nIinf,
				      CouNumber objVal,
				      CouNumber maxNLinf,
				      CouNumber maxIinf):
  x_          (CoinCopyOfArray (x, n)),
  n_          (n),
  nNLPinf_    (nNLPinf),
  nIinf_      (nIinf),
  objVal_     (objVal),
  maxNLinf_   (maxNLinf),
  maxIinf_    (maxIinf) {}


/// basic comparison procedure -- what to compare depends on user's choice
int CouenneFPsolution::compare (const CouenneFPsolution &other, enum what_to_compare comparedTerm) const {

  return 0;
}
