/*
 *
 * Name:    exprBMul.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of multiplications
 *
 * (C) Carnegie-Mellon University, 2006.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRBMUL_H
#define COUENNE_EXPRBMUL_H

#include "CouenneExprOp.hpp"
#include "CouenneConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

namespace Couenne {

#define MUL_ZERO 1e-20
#define MUL_INF  sqrt (COIN_DBL_MAX)

/// product that avoids NaN's
inline CouNumber safeProd (CouNumber a, CouNumber b) {

  if (a >  MUL_INF) return (b < -MUL_ZERO) ? -COIN_DBL_MAX : (b > MUL_ZERO) ?  COIN_DBL_MAX : 0.;
  if (a < -MUL_INF) return (b < -MUL_ZERO) ?  COIN_DBL_MAX : (b > MUL_ZERO) ? -COIN_DBL_MAX : 0.;

  if (b >  MUL_INF) return (a < -MUL_ZERO) ? -COIN_DBL_MAX : (a > MUL_ZERO) ?  COIN_DBL_MAX : 0.;
  if (b < -MUL_INF) return (a < -MUL_ZERO) ?  COIN_DBL_MAX : (a > MUL_ZERO) ? -COIN_DBL_MAX : 0.;

  return a*b;
}


/// class to compute lower bound of a product based on the bounds of
/// both factors

class COUENNELIB_EXPORT exprLBMul: public exprOp {

 public:

  /// Constructors, destructor
  exprLBMul  (expression **al, int n):
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprLBMul (clonearglist (d), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "LB_Mul";}
};


/// compute sum

inline CouNumber exprLBMul::operator () () {

  CouNumber n = (*(arglist_ [0])) ();
  CouNumber N = (*(arglist_ [1])) ();
  CouNumber d = (*(arglist_ [2])) ();
  CouNumber D = (*(arglist_ [3])) ();

  if (d>=0)
    if   (n>=0) return safeProd (n,d);
    else        return safeProd (n,D);
  else // d <= 0
    if (N>0) {
      CouNumber Nd = safeProd (N,d), nD;
      if (n<0 && D>0 &&
	  (Nd > (nD = safeProd (n,D)))) return nD;
      else                              return Nd;
    }
    else
      if (D>0) return safeProd (n,D);
      else     return safeProd (N,D);
}


/// class to compute upper bound of a product based on the bounds of
/// both factors

class COUENNELIB_EXPORT exprUBMul: public exprOp {

 public:

  /// Constructors, destructor
  exprUBMul  (expression **al, int n):
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprUBMul (clonearglist (d), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "UB_Mul";}
};


/// compute sum

inline CouNumber exprUBMul::operator () () {

  //  exprOp:: operator () ();

  CouNumber n = (*(arglist_ [0])) ();
  CouNumber N = (*(arglist_ [1])) ();
  CouNumber d = (*(arglist_ [2])) ();
  CouNumber D = (*(arglist_ [3])) ();

  if (d>0)
    if (N<0) return safeProd (N,d);
    else     return safeProd (N,D);
  else // d <= 0
    if (n<0) {
      CouNumber nd = safeProd (n,d), ND;
      if (N>0 && D>0 &&
	  ((ND = safeProd (N,D)) > nd)) return ND;
      else                              return nd;
    }
    else
      if (D>0) return safeProd (N,D);
      else     return safeProd (n,D);
}

}

#endif
