/* $Id$
 *
 * Name:    exprClone.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the clone class (different from exprCopy in
 *          that evaluation is propagated)
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRCLONE_HPP
#define COUENNE_EXPRCLONE_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "CouenneExprCopy.hpp"

namespace Couenne {

/// expression clone (points to another expression) 

class exprClone: public exprCopy {

 public:

  /// Constructor
  exprClone  (expression *copy): 
    exprCopy (getOriginal (copy)) {}

  /// copy constructor
  exprClone (const exprClone &e, Domain *d = NULL):
  //{copy_ = e.copy_;}// d = e.domain_;}
  exprCopy (e, d) {}

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprClone (*this, d);}
  /*{
    if (copy_ != Original ()) return copy_ -> clone (d);
    else                      return new exprClone (*this, d);
    }*/
  //{return copy_ -> Original () -> clone (d);}

  /// Destructor
  virtual ~exprClone () {
    //if (copy_ != Original ())
    //delete copy_;
    copy_ = NULL; // unlink pointer so that ~exprCopy does nothing
  }

  /// Printing
  virtual void print (std::ostream &out = std::cout, 
		      bool descend      = false) const;

  /// value
  inline CouNumber Value () const 
  {return copy_ -> Value ();}

  /// null function for evaluating the expression
  inline CouNumber operator () () 
  {return ((*copy_) ());}
};

}

#endif
