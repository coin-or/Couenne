/* $Id$
 *
 * Name:    CouennePSDcon.hpp
 * Author:  Pietro Belotti
 * Purpose: define the class of positive semidefinite constraints
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouennePSDcon_hpp
#define CouennePSDcon_hpp

#include "CouenneProblemElem.hpp"
#include <iostream>

namespace Couenne {

  class CouenneSparseMatrix;
  class CouenneProblem;
  class Domain;

  /// Class to represent positive semidefinite constraints //////////////////

  class CouennePSDcon: public CouenneConstraint {

  protected:

    CouenneSparseMatrix *X_; ///< contains indices of matrix X \succeq 0

  public:

    /// Constructor
    CouennePSDcon  (CouenneSparseMatrix *X): 
      CouenneConstraint (),
      X_                (X) {}

    /// Destructor
    ~CouennePSDcon ();

    /// Copy constructor
    CouennePSDcon (const CouennePSDcon &c, Domain *d = NULL);

    /// Cloning method
    inline CouenneConstraint *clone (Domain *d = NULL) const
    {return new CouennePSDcon (*this, d);}

    /// return X
    CouenneSparseMatrix *getX () const {return X_;}

    /// Decompose body of constraint through auxiliary variables
    exprAux *standardize (CouenneProblem *);

    /// Print constraint
    void print (std::ostream & = std::cout);
  };
}

#endif
