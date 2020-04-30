/* $Id$
 *
 * Name:    CouenneExprEvenPow.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers with even exponent
 *
 * (C) Pietro Belotti 2011
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPREVENPOW_HPP
#define COUENNE_EXPREVENPOW_HPP

#include <math.h>

#include "CouenneExprOp.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprConst.hpp"

namespace Couenne {

  /// Power of an expression (binary operator) with even exponent, \f$
  /// f(x)^k\f$ with \f$ k\in \mathbb Z\f$ constant even

  class exprEvenPow: public exprPow {

  public:

    /// Constructor
    exprEvenPow (expression **al, int n = 2): 
      exprPow (al, n) {} //< non-leaf expression, with argument list

    /// Constructor with only two arguments
    exprEvenPow (expression *arg0, expression *arg1):
      exprPow (arg0, arg1) {}

    /// cloning method
    expression *clone (Domain *d = NULL) const
    {return new exprEvenPow (clonearglist (d), nargs_);}

    /// function for the evaluation of the expression
    CouNumber operator () ();

    /// Get lower and upper bound of an expression (if any)
    void getBounds (expression *&, expression *&);

    /// Get value of lower and upper bound of an expression (if any)
    void getBounds (CouNumber &lb, CouNumber &ub);

    /// reduce expression in standard form, creating additional aux
    /// variables (and constraints)
    exprAux *standardize (CouenneProblem *p, bool addAux = true);

    /// generate equality between *this and *w
    void generateCuts (expression *w, //const OsiSolverInterface &si, 
		       OsiCuts &cs, const CouenneCutGenerator *cg, 
		       t_chg_bounds * = NULL, int = -1, 
		       CouNumber = -COUENNE_INFINITY, 
		       CouNumber =  COUENNE_INFINITY);

    /// return an index to the variable's argument that is better fixed
    /// in a branching rule for solving a nonconvexity gap
    expression *getFixVar () 
    {return arglist_ [0];}

    /// code for comparison
    virtual enum expr_type code () 
    {return COU_EXPREVENPOW;}

    /// implied bound processing
    bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

    /// set up branching object by evaluating many branching points for
    /// each expression's arguments
    virtual CouNumber selectBranch (const CouenneObject *obj, 
				    const OsiBranchingInformation *info,
				    expression * &var, 
				    double * &brpts, 
				    double * &brDist, // distance of current LP
				    // point to new convexifications
				    int &way);

    /// can this expression be further linearized or are we on its
    /// concave ("bad") side
    virtual bool isCuttable (CouenneProblem *problem, int index) const;
  };

  /// compute power
  inline CouNumber exprEvenPow::operator () () {
    //  return (currValue_ = safe_pow (base, exponent));
    return (safe_pow ((**arglist_) (), (*(arglist_ [1])) ()));
  }
}

#endif
