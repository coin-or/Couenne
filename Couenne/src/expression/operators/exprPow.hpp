/* $Id: exprPow.hpp 154 2009-06-16 18:52:53Z pbelotti $ */
/*
 * Name:    exprPow.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRPOW_HPP
#define COUENNE_EXPRPOW_HPP

#include <math.h>

#include "exprOp.hpp"
#include "exprMul.hpp"
#include "exprClone.hpp"
#include "exprConst.hpp"

class funtriplet;


/// Power of an expression (binary operator)

class exprPow: public exprOp {

 public:

  /// Constructor
  exprPow (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// Constructor with only two arguments
  exprPow (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprPow (clonearglist (d), nargs_);}

  /// print operator
  std::string printOp () const
    {return "^";}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// return l-2 norm of gradient at given point
  CouNumber gradientNorm (const double *x);

  /// differentiation
  expression *differentiate (int index); 

  /// simplification
  expression *simplify ();

  /// get a measure of "how linear" the expression is
  int Linearity ();

  /// is this expression integer?
  bool isInteger ();

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Get value of lower and upper bound of an expression (if any)
  void getBounds (CouNumber &lb, CouNumber &ub);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// generate equality between *this and *w
  void generateCuts (expression *w, const OsiSolverInterface &si, 
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
    {return COU_EXPRPOW;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

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
};


/// compute power and check for integer-and-odd inverse exponent

inline CouNumber safe_pow (CouNumber base, 
			   CouNumber exponent) {

  long double 
    lbase     = base,
    lexponent = exponent,
    retval    = 0.;

  if (lbase < 0.) {

    register int rndexp;

    if (((fabs (lexponent - (rndexp = COUENNE_round (lexponent))) < COUENNE_EPS) ||
	 ((fabs (lexponent) > COUENNE_EPS) && 
	  (fabs (1. / lexponent - (rndexp = COUENNE_round (1. / lexponent))) < COUENNE_EPS)))) {
      if (rndexp % 2)
	retval = (- pow (- lbase, lexponent)); // x^k, x negative, k odd
      else retval = pow (-lbase, lexponent);   // x^k, x negative, k even
    }
    else retval =  0.; // this is incorrect but avoids nan's
  }
  else if (fabs (lbase) >= COUENNE_INFINITY) {

    if (lbase <= -COUENNE_INFINITY) {

      register int intk = COUENNE_round (lexponent);

      if ((fabs (lexponent - intk) < COUENNE_EPS) && (intk % 2))
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
  return (safe_pow ((**arglist_) (), (*(arglist_ [1])) ()));
}


/// add upper/lower envelope to power in convex/concave areas
void addPowEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int,
		     CouNumber, CouNumber, CouNumber, CouNumber, CouNumber, int);


/// find proper tangent point to add deepest tangent cut
CouNumber powNewton (CouNumber, CouNumber, unary_function, unary_function, unary_function);

/// find proper tangent point to add deepest tangent cut
CouNumber powNewton (CouNumber, CouNumber, funtriplet *);

#endif
