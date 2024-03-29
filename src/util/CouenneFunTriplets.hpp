/*
 *
 * Name:    funtriplets.hpp
 * Author:  Pietro Belotti
 * Purpose: class for representing a function and its first- and second-order derivative
 *
 * (C) Carnegie-Mellon University, 2007-10
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef Funtriplets_hpp
#define Funtriplets_hpp

#include <math.h>

#include "CouenneExprPow.hpp"

namespace Couenne
{

///
class COUENNELIB_EXPORT funtriplet {

public:

  /// Basic constructor
  funtriplet () {}

  /// Destructor
  virtual ~funtriplet () {}

  virtual CouNumber F     (CouNumber x) = 0; //< main funtion
  virtual CouNumber Fp    (CouNumber x) = 0; //< first-order derivative of main funtion
  virtual CouNumber Fpp   (CouNumber x) = 0; //< second-order derivative of main funtion
  virtual CouNumber FpInv (CouNumber x) = 0; //< inverse of the first-order derivative
};


///
class COUENNELIB_EXPORT simpletriplet: public funtriplet {

protected:

  unary_function f_;   //< the function
  unary_function fp_;  //< the first-order derivative
  unary_function fpp_; //< the second-order derivative
  unary_function fpI_; //< the inverse of the first-order derivative

public:

  /// Basic constructor
  simpletriplet (unary_function f   = NULL,
		 unary_function fp  = NULL,
		 unary_function fpp = NULL,
		 unary_function fpI = NULL):
    f_   (f),
    fp_  (fp),
    fpp_ (fpp),
    fpI_ (fpI) {}

  /// Destructor
  virtual ~simpletriplet () {}

  virtual CouNumber F     (CouNumber x) {return f_   (x);} //< main funtion
  virtual CouNumber Fp    (CouNumber x) {return fp_  (x);} //< first-order derivative
  virtual CouNumber Fpp   (CouNumber x) {return fpp_ (x);} //< second-order derivative
  virtual CouNumber FpInv (CouNumber x) {return fpI_ (x);} //< inverse of first-order derivative
};


///
class COUENNELIB_EXPORT powertriplet: public funtriplet {

protected:

  CouNumber exponent_;    //< defines the power function triplet
  bool      issignpower_; //< determines if signed power

public:

  /// Basic constructor
  powertriplet (CouNumber exponent, bool signpower = false):
    exponent_ (exponent), issignpower_ (signpower) {}

  /// Destructor
  virtual ~powertriplet () {}

  virtual CouNumber F   (CouNumber x)
  {return safe_pow (x, exponent_, issignpower_);}                                   //< main funtion

  virtual CouNumber Fp  (CouNumber x)
  {return exponent_ * safe_pow (issignpower_ ? fabs(x) : x, exponent_ - 1);}  //< first-order derivative

  virtual CouNumber Fpp (CouNumber x)
  {return exponent_ * (exponent_ - 1) * safe_pow (x, exponent_ - 2, issignpower_);} //< second-order derivative

  virtual CouNumber FpInv (CouNumber x)
  {return safe_pow (x / exponent_, 1 / (exponent_ - 1), issignpower_);} //< inverse of first derivative
};


///
class COUENNELIB_EXPORT kpowertriplet: public powertriplet {

protected:

  CouNumber mult_; //< pre-multiplier

public:

  /// Basic constructor
  kpowertriplet (CouNumber exponent, CouNumber k):
    powertriplet (exponent),
    mult_ (k) {}

  /// Destructor
  virtual ~kpowertriplet () {}

  virtual CouNumber F   (CouNumber x)  //< main funtion
  {return mult_ * safe_pow (x, exponent_);}

  virtual CouNumber Fp  (CouNumber x)  //< first-order derivative
  {return mult_ * exponent_ * safe_pow (x, exponent_ - 1);}

  virtual CouNumber Fpp (CouNumber x)  //< second-order derivative
  {return mult_ * exponent_ * (exponent_ - 1) * safe_pow (x, exponent_ - 2);}

  virtual CouNumber FpInv (CouNumber x)
  {return safe_pow (x / (mult_ * exponent_), 1 / (exponent_ - 1));} //< inverse of first derivative
};

}
#endif
