/*
 *
 * Name:    CouennePSDcon.hpp
 * Author:  Pietro Belotti
 * Purpose: define the class of positive semidefinite constraints
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouennePSDcon_hpp
#define CouennePSDcon_hpp

#include "CouenneConfig.h"
#include "CouenneProblemElem.hpp"
#include <iostream>

namespace Couenne {

  class CouenneExprMatrix;
  class CouenneProblem;
  class Domain;

  /// Class to represent positive semidefinite constraints //////////////////

  class COUENNELIB_EXPORT CouennePSDcon: public CouenneConstraint {

  protected:

    CouenneExprMatrix *X_; ///< contains indices of matrix X \succeq 0

  public:

    /// Constructor
    CouennePSDcon  (CouenneExprMatrix *X): 
      CouenneConstraint (),
      X_                (X) {}

    /// Destructor
    ~CouennePSDcon ();

    /// Copy constructor
    CouennePSDcon (const CouennePSDcon &c, Domain *d = NULL);

    /// Assignment operator
    CouennePSDcon &operator= (const CouennePSDcon &c);

    /// Cloning method
    inline CouenneConstraint *clone (Domain *d = NULL) const
    {return new CouennePSDcon (*this, d);}

    /// return X
    CouenneExprMatrix *getX () const {return X_;}

    /// Decompose body of constraint through auxiliary variables
    exprAux *standardize (CouenneProblem *);

    /// Print constraint
    void print (std::ostream & = std::cout);
  };
}

#endif
