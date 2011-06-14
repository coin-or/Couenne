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
#include "CouenneExprVar.hpp"
#include "CouenneExprAux.hpp"

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

  for (std::vector <exprVar *>::iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end ();
       ++i) {

    CouNumber 
      vval = (**i) ();

    if ((*i) -> Multiplicity () <= 0)
      continue;

    if ((*i) -> isInteger ()) {

      double inf = CoinMax (vval - floor (vval + COUENNE_EPS),
			    ceil (vval - COUENNE_EPS) - vval);

      if (inf > COUENNE_EPS) {

	++nIinf_;

	if (inf > maxIinf_) 
	  maxIinf_ = inf;
      }
    }

    if (((*i) -> Type () == AUX) &&
	((*i) -> Image () -> Linearity () > LINEAR)) {

      double 
	diff = 0.,
	fval = (*((*i) -> Image ())) ();

      if      ((*i) -> sign () != expression::AUX_GEQ) diff = CoinMax (diff, vval - fval);
      else if ((*i) -> sign () != expression::AUX_LEQ) diff = CoinMax (diff, fval - vval);

      if (diff > COUENNE_EPS) {

	++nNLPinf_;

	if (diff > maxNLinf_)
	  maxNLinf_ = diff;
      }
    }
  }
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
