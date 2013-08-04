/* $Id$
 *
 * Name:    exprStore.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of a storage class for expressions
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRSTORE_HPP
#define COUENNE_EXPRSTORE_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "CouenneExprCopy.hpp"

namespace Couenne {

/// storage class for previously evaluated expressions

class exprStore: public exprCopy {

 protected:

  /// Value of the (previously evaluated) expression
  CouNumber value_;

 public:

  /// Constructor
  exprStore (expression *copy):
    exprCopy (copy) {}

  /// Store constructor -- Must go
  exprStore (const exprStore &e, Domain *d = NULL):
    exprCopy (e, d) {
    //copy_  = e.Original () -> clone ();
  }

  /// Destructor
  virtual ~exprStore () 
  {copy_ = NULL;}

  /// Printing
  virtual void print (std::ostream &out = std::cout, 
		      bool descend      = false) const;

  /// Cloning method
  virtual inline expression *clone (Domain *d = NULL) const
  {return new exprStore (*this, d);}

  /// function for evaluating the expression -- returns value of
  /// exprCopy pointed to, which returns a value stored from a
  /// previous evaluation
  virtual inline CouNumber operator () () 
  {return (copy_ -> Value ());}
};

}

#endif
