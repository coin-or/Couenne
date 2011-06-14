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
  nNLinf_     (0),
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

	++nNLinf_;

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
				      int nNLinf,
				      int nIinf,
				      CouNumber objVal,
				      CouNumber maxNLinf,
				      CouNumber maxIinf):

  x_          (CoinCopyOfArray (x, n)),
  n_          (n),
  nNLinf_     (nNLinf),
  nIinf_      (nIinf),
  objVal_     (objVal),
  maxNLinf_   (maxNLinf),
  maxIinf_    (maxIinf) {}

/// copy constructor
CouenneFPsolution::CouenneFPsolution (const CouenneFPsolution &src):
  x_          (src.x_ ? CoinCopyOfArray (src.x_, src.n_) : NULL),
  n_          (src.n_),
  nNLinf_     (src.nNLinf_),
  nIinf_      (src.nIinf_),
  objVal_     (src.objVal_),
  maxNLinf_   (src.maxNLinf_),
  maxIinf_    (src.maxIinf_) {}


/// assignment
CouenneFPsolution &CouenneFPsolution::operator= (const CouenneFPsolution &src) {

  x_         = src.x_ ? CoinCopyOfArray (src.x_, src.n_) : NULL;
  n_         = src.n_;
  nNLinf_    = src.nNLinf_;
  nIinf_     = src.nIinf_;
  objVal_    = src.objVal_;
  maxNLinf_  = src.maxNLinf_;
  maxIinf_   = src.maxIinf_;

  return *this;
}

/// destructor
CouenneFPsolution::~CouenneFPsolution () {

  if (x_)
    delete [] x_;
}

/// basic comparison procedure -- what to compare depends on user's choice
bool CouenneFPsolution::compare (const CouenneFPsolution &other, enum what_to_compare comparedTerm) const {

  switch (comparedTerm) {

  case SUM_NINF: return (nNLinf_   + nIinf_   < other.nNLinf_   + other.nIinf_);
  case SUM_INF:  return (maxNLinf_ + maxIinf_ < other.maxNLinf_ + other.maxIinf_);
  case OBJVAL:   return (objVal_              < other.objVal_);
  }

  printf ("CouenneFPsolution::compare: bad compared term\n");
  abort ();
}


/// copy constructor
CouenneFPpool::CouenneFPpool (const CouenneFPpool &src):
  queue_ (src.queue_) {}


/// assignment
CouenneFPpool &CouenneFPpool::operator= (const CouenneFPpool &src) {

  queue_ = src.queue_;

  return *this;
}

