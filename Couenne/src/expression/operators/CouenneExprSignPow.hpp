/* $Id$
 *
 * Name:    CouenneExprSignPow.hpp
 * Author:  Pietro Belotti
 * Purpose: Defines signed powers, i.e., functions of the form x|x|^k with k in R
 *
 * (C) Pietro Belotti 2011
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRSIGNPOW_HPP
#define COUENNE_EXPRSIGNPOW_HPP

#include <math.h>

#include "CouenneExprPow.hpp"

namespace Couenne {

  class funtriplet;

  /// Power of an expression (binary operator), \f$
  /// f(x)|f(x)|^{k-1}\f$ with \f$ k \in \mathbb R\f$ constant

  class exprSignPow: public exprSignPow {

  public:

    /// Constructor
    exprSignPow (expression **al, int n = 2): 
      exprPow (al, n) {} //< non-leaf expression, with argument list

    /// Constructor with only two arguments
    exprSignPow (expression *arg0, expression *arg1):
      exprPow (arg0, arg1) {}

    /// cloning method
    expression *clone (Domain *d = NULL) const
    {return new exprSignPow (clonearglist (d), nargs_);}

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

    /// code for comparison
    virtual enum expr_type code () 
    {return COU_EXPRSIGNPOW;}

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
  inline CouNumber exprSignPow::operator () () {
    //  return (currValue_ = safe_pow (base, exponent));
    return (safe_pow ((**arglist_) (), (*(arglist_ [1])) ()));
  }
}

#endif
