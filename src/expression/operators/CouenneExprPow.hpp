/*
 *
 * Name:    CouenneExprPow.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers, specifically with fractional
 *          exponent (odd/even/signed: see appropriate files)
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRPOW_HPP
#define COUENNE_EXPRPOW_HPP

#include <math.h>
#include <stdio.h>

#include "CouenneExprOp.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprConst.hpp"

namespace Couenne {

class funtriplet;


/// Power of an expression (binary operator), \f$ f(x)^k\f$ with \f$ k\f$ constant

class COUENNELIB_EXPORT exprPow: public exprOp {

 private:

  /// do we mean a signed power function: sign(arg0) * |arg0|^arg1 (assumes that arg1 is constant) 
  bool issignpower_;

 public:

  /// Constructor
  exprPow (expression **al, int n = 2, bool signpower = false): 
    exprOp (al, n), issignpower_ (signpower) {} //< non-leaf expression, with argument list

  /// Constructor with only two arguments
  exprPow (expression *arg0, expression *arg1, bool signpower = false):
    exprOp (arg0, arg1), issignpower_ (signpower) {}

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprPow (clonearglist (d), nargs_, issignpower_);}

  /// print operator positioning
  virtual enum pos printPos () const
  {return issignpower_ ? PRE : INSIDE ;}

  /// print operator
  virtual std::string printOp () const
  {return issignpower_ ? "signpower" : "^";}

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// return l-2 norm of gradient at given point
  virtual CouNumber gradientNorm (const double *x);

  /// differentiation
  virtual expression *differentiate (int index); 

  /// simplification
  virtual expression *simplify ();

  /// get a measure of "how linear" the expression is
  virtual int Linearity ();

  /// is this expression integer?
  virtual bool isInteger ();

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Get value of lower and upper bound of an expression (if any)
  virtual void getBounds (CouNumber &lb, CouNumber &ub);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// generate equality between *this and *w
  virtual void generateCuts (expression *w, //const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar () 
  {return arglist_ [0];}

  /// code for comparison
  virtual enum expr_type code () 
  {return (issignpower_ ? COU_EXPRSIGNPOW : COU_EXPRPOW);}

  /// implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (const CouenneObject *obj, 
				  const OsiBranchingInformation *info,
				  expression * &var, 
				  double * &brpts, 
 				  double * &brDist, // distance of current LP
					  	    // point to new convexifications
				  int &way);

  /// compute $y^{lv}$ and $y^{uv}$ for Violation Transfer algorithm
  virtual void closestFeasible (expression *varind,
				expression *vardep, 
				CouNumber &left,
				CouNumber &right) const;

  /// can this expression be further linearized or are we on its
  /// concave ("bad") side
  virtual bool isCuttable (CouenneProblem *problem, int index) const;

  /// return whether this expression corresponds to a signed integer power 
  virtual bool isSignpower () const { return issignpower_; } 
};


/// compute power and check for integer-and-odd inverse exponent

inline CouNumber safe_pow (CouNumber base, 
			   CouNumber exponent, bool signpower = false) {

  double 
    lbase     = base,
    lexponent = exponent,
    retval    = 0.;

  if (lbase < 0.) {

    int rndexp = COUENNE_round (lexponent);

    if (((fabs (lexponent - rndexp) < COUENNE_EPS) ||
	 ((fabs (lexponent) > COUENNE_EPS) && 
	  (fabs (1. / lexponent - (rndexp = COUENNE_round (1. / lexponent))) < COUENNE_EPS)))) {
      if ((rndexp % 2) || signpower)
	retval = (- pow (- lbase, lexponent)); // x^k, x negative, k odd  or signed power
      else retval = pow (- lbase, lexponent);  // x^k, x negative, k even or signed power
    }
    else retval =  0.; // this is incorrect but avoids nan's
  }
  else if (fabs (lbase) >= COUENNE_INFINITY) {

    if (lbase <= -COUENNE_INFINITY) {

      int intk = COUENNE_round (lexponent);

      if ((fabs (lexponent - intk) < COUENNE_EPS) && (intk % 2 || signpower))
	retval = (lexponent < 0.) ? 0. : -COUENNE_INFINITY;
    }
    else retval = (lexponent < 0.) ? 0. : COUENNE_INFINITY;
  }
  else
    retval = (pow (lbase, lexponent));

  return (CouNumber) (retval);
}


/// compute power
inline CouNumber exprPow::operator () () {
  //  return (currValue_ = safe_pow (base, exponent));
  return safe_pow ((**arglist_) (), (*(arglist_ [1])) (), issignpower_);
}


/// add upper/lower envelope to power in convex/concave areas
COUENNELIB_EXPORT
void addPowEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int,
		     CouNumber, CouNumber, CouNumber, CouNumber, CouNumber, int, bool = false);

/// find proper tangent point to add deepest tangent cut
COUENNELIB_EXPORT
CouNumber powNewton (CouNumber, CouNumber, unary_function, unary_function, unary_function);

/// find proper tangent point to add deepest tangent cut
COUENNELIB_EXPORT
CouNumber powNewton (CouNumber, CouNumber, funtriplet *);

}

#endif
