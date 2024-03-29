/*
 *
 * Name:    exprSin.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRSIN_HPP
#define COUENNE_EXPRSIN_HPP

#include <math.h>
#include <assert.h>

#include "CouenneExprUnary.hpp"
#include "CouenneExprConst.hpp"

namespace Couenne {

/// specify which trigonometric function is dealt with in trigEnvelope
enum cou_trig {COU_SINE, COU_COSINE};


/// normalize angle within [0,b] (typically, pi or 2pi)
inline CouNumber modulo (CouNumber a, CouNumber b)
  {return a - b * floor (a/b);}


/// generalized procedure for both sine and cosine
COUENNELIB_EXPORT
CouNumber trigSelBranch (const CouenneObject *obj,
			 const OsiBranchingInformation *info,
			 expression * &var,
			 double * &brpts,
			 double * &brDist, // distance of current LP
				           // point to new convexifications
			 int &way,
			 enum cou_trig type);


/// generalized implied bound procedure for sine/cosine
COUENNELIB_EXPORT
bool trigImpliedBound (enum cou_trig, int, int, CouNumber *, CouNumber *, t_chg_bounds *);


/// class for \f$ \sin f(x)\f$
class COUENNELIB_EXPORT exprSin: public exprUnary {

 public:

  /// Constructors, destructor
  exprSin (expression *al):
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprSin (argument_ -> clone (d));}

  //// the operator's function
  inline unary_function F ()
  {return sin;}

  /// print operator
  std::string printOp () const
  {return "sin";}

  /// return l-2 norm of gradient at given point
  inline CouNumber gradientNorm (const double *x) {
    return (argument_ -> Index () < 0) ?
      0. : fabs (cos (x [argument_ -> Index ()]));
  }

  /// differentiation
  expression *differentiate (int index);

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Get value of lower and upper bound of an expression (if any)
  void getBounds (CouNumber &lb, CouNumber &ub);

  /// generate equality between *this and *w
  void generateCuts (expression *w, //const OsiSolverInterface &si,
		     OsiCuts &cs, const CouenneCutGenerator *cg,
		     t_chg_bounds * = NULL, int = -1,
		     CouNumber = -COUENNE_INFINITY,
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code ()
  {return COU_EXPRSIN;}

  /// implied bound processing
  bool impliedBound (int index, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign = expression::AUX_EQ) {

    bool impl = trigImpliedBound (COU_SINE, index, argument_ -> Index (), l, u, chg);

    if (impl && argument_ -> isInteger ()) {

      int ind = argument_ -> Index ();
      assert (ind >= 0);
      l [ind] = ceil  (l [ind] - COUENNE_EPS);
      u [ind] = floor (u [ind] + COUENNE_EPS);
    }

    return impl;
  }

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (const CouenneObject *obj,
				  const OsiBranchingInformation *info,
				  expression * &var,
				  double * &brpts,
 				  double * &brDist, // distance of current LP
					  	    // point to new convexifications
				  int &way)

  {return trigSelBranch (obj, info, var, brpts, brDist, way, COU_SINE);}

  /// closest feasible points in function in both directions
  virtual void closestFeasible (expression *varind, expression *vardep,
				CouNumber& left, CouNumber& right) const;

  /// can this expression be further linearized or are we on its
  /// concave ("bad") side
  virtual bool isCuttable (CouenneProblem *problem, int index) const
  {return false;}
};

}

#endif
