/*
 *
 * Name:    exprConst.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprConst
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRCONST_HPP
#define COUENNE_EXPRCONST_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"

namespace Couenne {

/// constant-type operator

class COUENNELIB_EXPORT exprConst: public expression {

private:

  /// the value of this constant
  CouNumber value_;

public:

  /// node type
  inline enum nodeType Type () const
  {return CONST;}

  /// value of expression
  inline CouNumber Value () const
  {return value_;}

  /// Constructor
  exprConst (CouNumber value):
    value_ (value) {}

  /// Copy constructor
  exprConst (const exprConst &e, Domain *d = NULL)
  {value_ = e.value_;}

  /// Cloning method
  virtual inline expression *clone (Domain *d = NULL) const
  {return new exprConst (value_);}

  /// I/O
  void print (std::ostream &out = std::cout,
	      bool = false) const
  {out << value_;}

  /// return constant's value
  inline CouNumber operator() ()
  {return value_;}

  /// differentiation
  inline expression *differentiate (int)
  {return new exprConst (0.);}

  /// dependence on variable set
  inline int dependsOn (int *ind, int n, enum dig_type type = STOP_AT_AUX)
  {return 0;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
  {return ((fabs (value_) < COUENNE_EPS) ? ZERO: CONSTANT);}

  /// Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) {
    lower = new exprConst (value_);
    upper = new exprConst (value_);
  }

  /// Get value of lower and upper bound of an expression (if any)
  inline void getBounds (CouNumber &lower, CouNumber &upper)
  {lower = upper = value_;}

  /// generate convexification cut for constraint w = this
  void generateCuts (expression *, //const OsiSolverInterface &,
		     OsiCuts &, const CouenneCutGenerator *,
		     t_chg_bounds * = NULL, int = -1,
		     CouNumber = -COUENNE_INFINITY,
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual inline enum expr_type code ()
  {return COU_EXPRCONST;}

  /// is this expression integer?
  virtual inline bool isInteger ()
  {return Couenne::isInteger (value_);}

  /// used in rank-based branching variable choice
  virtual inline int rank ()
  {return 0;}
};

}

#endif
