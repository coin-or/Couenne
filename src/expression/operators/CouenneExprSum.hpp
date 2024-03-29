/*
 *
 * Name:    exprSum.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRSUM_H
#define COUENNE_EXPRSUM_H

#include <vector>

#include "CouenneExprOp.hpp"

namespace Couenne {

/// class sum, \f$ \sum_{i=1}^n f_i(x) \f$

class COUENNELIB_EXPORT exprSum: public exprOp {

 public:

  /// Constructors, destructor
  exprSum  (expression ** = NULL, int = 0);

  /// Constructor with two elements
  exprSum (expression *, expression *);

  /// Empty destructor
  virtual ~exprSum () {}

  /// Cloning method
  virtual expression *clone (Domain *d = NULL) const
    {return new exprSum (clonearglist (d), nargs_);}

  /// Print operator
  std::string printOp () const
    {return "+";}

  /// Function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// Differentiation
  virtual expression *differentiate (int index);

  /// Simplification
  virtual expression *simplify ();

  /// Get a measure of "how linear" the expression is:
  virtual int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (CouNumber &, CouNumber &);

  /// Reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// Special version for linear constraints
  virtual void generateCuts (expression *, //const OsiSolverInterface &,
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY,
			     CouNumber =  COUENNE_INFINITY);

  /// Code for comparison
  virtual enum expr_type code ()
    {return COU_EXPRSUM;}

  /** Implied bound.
   *  An expression
   *
   *  \f$w = a0 + \sum_{i\in I1} a_i x_i + \sum_{i\in I2} a_i x_i\f$
   *
   *  is given such that all \f$a_i\f$ are positive for \f$i \in I1\f$ and
   *  negative for \f$i \in I2\f$. If the bounds on \f$w \in [l,u]\f$, implied
   *  bounds on all \f$x_i, i\in I1 \cup I2\f$ are as follows:
   *
   *  \f$\forall i\in I1\f$
   *    \f$x_i \ge (l - a0 - \sum_{j\in I1 | j\neq i} a_j u_j - \sum_{j\in I2}        a_j l_j) / a_i\f$
   *    \f$x_i \le (u - a0 - \sum_{j\in I1 | j\neq i} a_j l_j - \sum_{j\in I2}        a_j u_j) / a_i\f$
   *
   *  \f$\forall i\in I2\f$
   *    \f$x_i \ge (u - a0 - \sum_{j\in I1} a_j u_j       - \sum_{j\in I2 | j\neq i} a_j l_j) / a_i\f$
   *    \f$x_i \le (l - a0 - \sum_{j\in I1} a_j l_j       - \sum_{j\in I2 | j\neq i} a_j u_j) / a_i\f$,
   *
   *  where \f$l_i\f$ and \f$u_i\f$ are lower and upper bound,
   *  respectively, of \f$x_i\f$. We also have to check if some of
   *  these bounds are infinite.
   */
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

  /// Checks for quadratic terms in the expression and returns an
  /// exprQuad if there are enough to create something that can be
  /// convexified
  exprAux *createQuadratic (CouenneProblem *);

protected:

  /// inferring bounds on factors of a product
  int impliedBoundSum (CouNumber wl,
		       CouNumber wu,
		       std::vector <CouNumber> &xl,
		       std::vector <CouNumber> &xu,
		       std::vector <std::pair <int, CouNumber> > &nl,
		       std::vector <std::pair <int, CouNumber> > &nu);
};


/// compute sum

inline CouNumber exprSum::operator () () {

  CouNumber ret = 0;

  expression **al = arglist_;

  for (int n = nargs_; n--;)
    ret += (**al++) ();

  return ret;
}

}

#endif
