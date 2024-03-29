/*
 *
 * Name:    CouenneProblemElem.hpp
 * Author:  Pietro Belotti
 * Purpose: define the classes used by class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_PROBLEM_ELEM_HPP
#define COUENNE_PROBLEM_ELEM_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"

namespace Couenne {

/** Class to represent nonlinear constraints
 *
 *  It consists of an expression as the body and two range expressions
 *  as lower- and upper bounds.
 *
 *  A general constraint is defined as lb_ <= body_ <= ub_, where all
 *  three components are expressions, depending on variables,
 *  auxiliaries and bounds. If the constraint is 2 <= exp (x1+x2) <=
 *  4, then:
 *
 *  body_ = exp (x1+x2), that is,
 *
 *  new exprExp (new exprSum (new exprVar (1), new exprVar (2))
 *
 *  while lb_ = new exprConst (2.) and ub_ = new exprConst (4.).
 */

class COUENNELIB_EXPORT CouenneConstraint {

 protected:

  expression *body_; ///< Body of constraint
  expression *lb_;   ///< Lower bound (expression)
  expression *ub_;   ///< Upper bound (expression)

 public:

  /// Constructor
  CouenneConstraint  (expression *body = NULL,
  	              expression *lb   = NULL,
		      expression *ub   = NULL):
    body_     (body),
    lb_       (lb),
    ub_       (ub) {

    if (!lb_)
      if (!ub_) {
	lb_ = new exprConst (0.);
	ub_ = new exprConst (0.);
      }
      else         lb_ = new exprConst (- COUENNE_INFINITY);
    else if (!ub_) ub_ = new exprConst   (COUENNE_INFINITY);
  }

  /// Destructor
  virtual ~CouenneConstraint () {
    delete body_;
    delete lb_;
    delete ub_;
  }

  /// Copy constructor
  CouenneConstraint  (const CouenneConstraint &c, Domain *d = NULL):
    body_  (c.Body () -> clone (d)),
    lb_    (c.Lb   () -> clone (d)),
    ub_    (c.Ub   () -> clone (d)) {}

  /// Cloning method
  virtual inline CouenneConstraint *clone (Domain *d = NULL) const
  {return new CouenneConstraint (*this, d);}

  // Get constraint's elements
  virtual inline expression *Lb   () const {return lb_;}   ///< Expression of lower bound
  virtual inline expression *Ub   () const {return ub_;}   ///< Expression of upper bound
  virtual inline expression *Body () const {return body_;} ///< Expression of body of constraint

  /// Set body of constraint
  virtual inline expression *Body (expression *newBody)
  {body_ = newBody; return body_;}

  /// decompose body of constraint through auxiliary variables
  virtual exprAux *standardize (CouenneProblem *);

  /// print constraint
  virtual void print (std::ostream & = std::cout);
};



/**
 * Objective function
 *
 * It consists of an expression only. We only assume minimization
 * problems (proper sign changes are applied upon reading)
 *
 */

class COUENNELIB_EXPORT CouenneObjective {

 protected:

  /// expression to optimize
  expression *body_;

 public:

  /// constructor
  CouenneObjective (expression *body):
    body_ (body) {}

  /// destructor
  ~CouenneObjective ()
  {delete body_;}

  /// copy constructor
  CouenneObjective  (const CouenneObjective &o, Domain *d = NULL):
    body_  (o.body_ -> clone (d)) {}

  /// cloning method
  inline CouenneObjective *clone (Domain *d = NULL) const
  {return new CouenneObjective (*this, d);}

  /// get body
  inline expression *Body () const
  {return body_;}

  /// Set body
  expression *Body (expression *newBody)
  {body_ = newBody; return body_;}

  /// Get standard form of this objective function
  inline exprAux *standardize (CouenneProblem *p)
  {return body_ -> standardize (p);}

  /// Print to iostream
  void print (std::ostream &out = std::cout) {
    out << "min ";
    body_ -> print (out);
    out << std::endl;
  }
};

}

#endif
