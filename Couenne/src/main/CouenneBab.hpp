/* $Id$
 *
 * Name:    CouenneBab.hpp
 * Author:  Pietro Belotti
 * Purpose: B&B object  
 * Created: 2012-01-25
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEBAB_H
#define COUENNEBAB_H

#include "BonCbc.hpp"
#include "BonBabSetupBase.hpp"

namespace Couenne {

  class CouenneProblem;

  class CouenneBab: public Bonmin::Bab {

  public:

    CouenneBab (); ///< Constructor

    virtual ~CouenneBab(); ///< Destructor

    void setProblem (CouenneProblem *p);

    virtual void branchAndBound (Bonmin::BabSetupBase &s); ///< Carry out branch and bound

    /// Get the best solution known to the problem (is optimal if
    /// MipStatus is FeasibleOptimal).  If no solution is known
    /// returns NULL.
    const double * bestSolution() const;

    /// Return objective value of the bestSolution
    double bestObj() const;

  protected:

    CouenneProblem *problem_;

  };
}

#endif
