/* $Id$
 *
 * Name:    CouenneGlobalCutOff.hpp
 * Author:  Pietro Belotti, Lehigh University
 *          Andreas Waechter, IBM
 * Purpose: a cutoff that replicates itself
 *
 * (C) Carnegie-Mellon University, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_GLOBAL_CUTOFF_HPP
#define COUENNE_GLOBAL_CUTOFF_HPP

#include "CouenneProblem.hpp"

namespace Couenne {

  class GlobalCutOff {

  private:

    GlobalCutOff (const GlobalCutOff&);

    double  cutoff_; ///< Value of the best solution
    double *sol_;    ///< Best solution
    int     size_;   ///< Size of the vector stored in sol (should be #var of reformulation)
    bool    valid_;  ///< Stored solution corresponds to cutoff

  public:

    GlobalCutOff ();
    GlobalCutOff (double c, const double *s=NULL, int n=0);
    ~GlobalCutOff ();

    void setCutOff (const CouenneProblem *p, double cutoff, const double *s=NULL);

    inline double  getCutOff    () const {return cutoff_;}
    inline double *getCutOffSol () const {return sol_;}
  };
}

#endif
