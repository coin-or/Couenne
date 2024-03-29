/*
 *
 * Name:    exprMin.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of $\f{\textrm{argmin}_{i\in I} y_i}$
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRMIN_H
#define COUENNE_EXPRMIN_H

#include "CouenneExprOp.hpp"
#include "CouenneExprCopy.hpp"
#include "CouenneExprStore.hpp"

// TODO: CouenneObject might take expression of the form
// min(x1,x2,...,xn) and branch it by creating two subsets N1 and N2
// of N={1,2,...,n} and then impose the constraints
//
// (xi <= xj for i in N1, j in N2) OR
// (xi >= xj for i in N1, j in N2)

namespace Couenne {

/// class for minima

class COUENNELIB_EXPORT exprMin: public exprOp {

 public:

  /// Constructor
  exprMin  (expression **al, int n):
    exprOp (al, n) {}

  /// Constructor with only two arguments
  exprMin  (expression *el0, expression *el1):
    exprOp (new expression * [4], 4) {
    arglist_ [0] = new exprCopy (el0); arglist_ [1] = new exprStore (arglist_ [0]);
    arglist_ [2] = new exprCopy (el1); arglist_ [3] = new exprStore (arglist_ [2]);
  }

  /// Cloning method
  exprMin *clone (Domain *d = NULL) const
    {return new exprMin (clonearglist (d), nargs_);}

  /// Print operator
  std::string printOp () const
    {return "min";}

  /// Print operator
  enum pos printPos () const
    {return PRE;}

  /// Function for the evaluation of the expression
  CouNumber operator () ();

  /// Differentiation
  inline expression *differentiate (int)
    {return NULL;}

  /// Simplification
  inline expression *simplify ()
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return NONLINEAR;}

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual inline exprAux *standardize (CouenneProblem *, bool addAux = true)
    {return NULL;}

  /// Generate equality between *this and *w
  void generateCuts (expression *w, //const OsiSolverInterface &si,
		     OsiCuts &cs, const CouenneCutGenerator *cg,
		     t_chg_bounds * = NULL, int = -1,
		     CouNumber = -COUENNE_INFINITY,
		     CouNumber =  COUENNE_INFINITY);

  /// Code for comparisons
  virtual enum expr_type code ()
  {return COU_EXPRMIN;}
};


/// Compute minimum

inline CouNumber exprMin::operator () () {

  CouNumber best_val = (*(arglist_ [0])) ();
  int best_ind = 0;

  for (int ind = 2; ind < nargs_; ind += 2) {

    CouNumber val = (*(arglist_ [ind])) ();

    if (val < best_val) {
      best_ind = ind;
      best_val = val;
    }
  }

  return (*(arglist_ [best_ind + 1])) ();
}

}

#endif
