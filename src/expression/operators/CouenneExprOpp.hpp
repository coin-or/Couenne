/*
 *
 * Name:    exprOpp.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPROPP_HPP
#define COUENNE_EXPROPP_HPP

#include "CouennePrecisions.hpp"
#include "CouenneExprUnary.hpp"

namespace Couenne {

/// operator opp: returns the opposite of a number

inline CouNumber opp (CouNumber arg)
{return - arg;}


/// class opposite, \f$ -f(x) \f$

class COUENNELIB_EXPORT exprOpp: public exprUnary {

 public:

  /// Constructors, destructor
  exprOpp (expression *al):
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprOpp (argument_ -> clone (d));}

  /// the operator's function
  inline unary_function F ()
    {return opp;}

  /// Output
  void print (std::ostream &out,
	      bool descend) const;

  /// return l-2 norm of gradient at given point
  inline CouNumber gradientNorm (const double *x)
  {return (argument_ -> Index () < 0) ? 0. : 1.;}

  /// differentiation
  expression *differentiate (int index);

  /// simplification
  virtual expression *simplify ();

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return argument_ -> Linearity ();}

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Get value of lower and upper bound of an expression (if any)
  void getBounds (CouNumber &, CouNumber&);

  /// special version for linear constraints
  virtual void generateCuts (expression *, //const OsiSolverInterface &,
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY,
			     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code ()
    {return COU_EXPROPP;}

  /// is this expression integer?
  bool isInteger ()
    {return argument_ -> isInteger ();}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

  /// standardization (to deal with complex arguments)
  exprAux *standardize (CouenneProblem *, bool addAux = true);
};

}

#endif
