/*
 *
 * Name:    exprDiv.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of divisions
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRDIV_HPP
#define COUENNE_EXPRDIV_HPP

#include "CouenneExprOp.hpp"
#include "CouennePrecisions.hpp"

namespace Couenne {

#define BR_NEXT_ZERO 1e-3
#define BR_MULT      1e-3

/// class for divisions, \f$ \frac{f(x)}{g(x)} \f$

class COUENNELIB_EXPORT exprDiv: public exprOp {

 public:

  /// Constructor
  exprDiv (expression **al, int n = 2):
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// Constructor with two arguments given explicitly
  exprDiv (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  /// Cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprDiv (clonearglist (d), nargs_);}

  /// Print operator
  std::string printOp () const
    {return "/";}

  /// Function for the evaluation of the expression
  inline CouNumber operator () ();

  /// return l-2 norm of gradient at given point
  CouNumber gradientNorm (const double *x);

  /// Differentiation
  expression *differentiate (int index);

  /// Simplification
  expression *simplify ();

  /// Get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity () {

    if (arglist_ [1] -> Type () == CONST)
      return arglist_ [0] -> Linearity ();
    else return NONLINEAR;
  }

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&lb, expression *&ub);

  /// Get value of lower and upper bound of an expression (if any)
  void getBounds (CouNumber &lb, CouNumber &ub);

  /// Reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// Generate equality between *this and *w
  void generateCuts (expression *w, //const OsiSolverInterface &si,
		     OsiCuts &cs, const CouenneCutGenerator *cg,
		     t_chg_bounds * = NULL, int = -1,
		     CouNumber = -COUENNE_INFINITY,
		     CouNumber =  COUENNE_INFINITY);

  /// Code for comparisons
  virtual enum expr_type code () {return COU_EXPRDIV;}

  /// is this expression integer?
  bool isInteger ();

  /// Implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (const CouenneObject *obj,
				  const OsiBranchingInformation *info,
				  expression * &var,
				  double * &brpts,
 				  double * &brDist, // distance of current LP
					  	    // point to new convexifications
				  int &way);

  /// compute $y^{lv}$ and  $y^{uv}$ for Violation Transfer algorithm
  virtual void closestFeasible (expression *varind,
				expression *vardep,
				CouNumber &left,
				CouNumber &right) const;

  /// can this expression be further linearized or are we on its
  /// concave ("bad") side
  virtual bool isCuttable (CouenneProblem *problem, int index) const
  {return false;} // concave on both sides, as for products
};


/// Compute division

inline CouNumber exprDiv::operator () ()
  {return ((*(*arglist_)) () / (*(arglist_ [1])) ());}


#define SAFE_COEFFICIENT 1e9

/// check if bounding box is suitable for a multiplication/division
/// convexification constraint

inline bool is_boundbox_regular (CouNumber b1, CouNumber b2) {

  // Why SAFE_COEFFICIENT and not COUENNE_INFINITY? Because
  // OsiRowCut::set[LU]b do not work for values more than
  // SAFE_COEFFICIENT and apparently makes the convexification
  // infeasible.
  return
    (fabs (b1)    < SAFE_COEFFICIENT) &&
    (fabs (b2)    < SAFE_COEFFICIENT) &&
    (fabs (b1*b2) < SAFE_COEFFICIENT);
    //    && ((fabs (b1) > COUENNE_EPS) || (fabs (b2) > COUENNE_EPS));
}

}

#endif
