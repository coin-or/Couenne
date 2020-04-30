/* $Id$
 *
 * Name:    CouenneGlobalCutOff.cpp
 * Author:  Pietro Belotti, Lehigh University
 *          Andreas Waechter, IBM
 * Purpose: a cutoff that replicates itself, implementation
 *
 * (C) Carnegie-Mellon University, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneGlobalCutOff.hpp"

#include "CoinFinite.hpp"
#include "CoinHelperFunctions.hpp"

using namespace Couenne;

GlobalCutOff::GlobalCutOff (): 
  cutoff_ (COIN_DBL_MAX), 
  sol_    (NULL), 
  size_   (0), 
  valid_  (false) {}

GlobalCutOff::GlobalCutOff (double c, const double *s, int n): 
  cutoff_ (c),
  sol_    (NULL),
  size_   (n),
  valid_  (false) {
  if (s) {
    sol_   = CoinCopyOfArray (s, n);
    size_  = n;
    valid_ = true;
  }
}


GlobalCutOff::~GlobalCutOff () 
{if (sol_) delete [] sol_;}


void GlobalCutOff::setCutOff (const CouenneProblem *p, double cutoff, const double *s) {

  cutoff_ = cutoff;

  valid_ = (s != NULL);

  if (s) {

    if (!sol_) 
      sol_ = new CouNumber [size_ = p -> nVars ()];

    CoinCopyN (s, p -> nOrigVars (), sol_); // fill first variables with values from NLP
    CoinFillN (sol_ + p -> nOrigVars (), p -> nVars () - p -> nOrigVars (), 0.); // pad with zeros

    p -> getAuxs (sol_);
  }
}
