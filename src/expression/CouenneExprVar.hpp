/*
 *
 * Name:    exprVar.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprVar for variables
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRVAR_HPP
#define COUENNE_EXPRVAR_HPP

#include <iostream>
#include <set>

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneDomain.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
  class Journalist;
}

namespace Bonmin {
  class BabSetupBase;
}

namespace Couenne {

  typedef Ipopt::SmartPtr<Ipopt::Journalist> JnlstPtr;
  typedef Ipopt::SmartPtr<const Ipopt::Journalist> ConstJnlstPtr;


class CouenneCutGenerator;

/// variable-type operator
///
/// All variables of the expression must be objects of this class or
/// of the derived exprAux class

class COUENNELIB_EXPORT exprVar: public expression {

 protected:

  int varIndex_;   ///< The index of the variable.
  Domain *domain_; ///< Pointer to a descriptor of the current point/bounds.

 public:

  /// Node type.
  virtual inline enum nodeType Type () const
  {return VAR;}

  /// Constructor.
  exprVar (int varIndex, Domain *d = NULL):
    varIndex_ (varIndex),
    domain_   (d) {}

  /// destructor
  virtual ~exprVar () {}

  /// Copy constructor
  exprVar (const exprVar &e, Domain *d = NULL):
    varIndex_ (e.varIndex_),
    domain_   (d) {}

  /// Cloning method
  virtual inline exprVar *clone (Domain *d = NULL) const
  {return new exprVar (*this, d);}

  /// Get variable index in problem
  inline int Index () const
  {return varIndex_;}

  // Bounds
  virtual expression *Lb (); ///< Get lower bound expression
  virtual expression *Ub (); ///< Get upper bound expression

  // Bounds
  virtual inline CouNumber &lb () {return domain_ -> lb (varIndex_);} ///< Get/set lower bound value
  virtual inline CouNumber &ub () {return domain_ -> ub (varIndex_);} ///< Get/set upper bound value

  /// print
  virtual void print (std::ostream &out = std::cout,
		      bool = false) const
  {out << "x_" << varIndex_;}

  /// return the value of the variable
  virtual inline CouNumber operator () ()
  {return domain_ -> x (varIndex_);}

  /// return l-2 norm of gradient at given point
  virtual inline CouNumber gradientNorm (const double *x)
  {return 1.;}

  /// differentiation
  virtual inline expression *differentiate (int index)
  {return new exprConst ((index == varIndex_) ? 1. : 0.);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual inline int DepList (std::set <int> &deplist,
			      enum dig_type type = ORIG_ONLY) {

    if (deplist.find (varIndex_) == deplist.end ()) {
      deplist.insert (varIndex_);
      return 1;
    }
    return 0;
  }

  /// set bounds depending on both branching rules and propagated
  /// bounds. To be used after standardization
  virtual inline void crossBounds () {}

  /// simplify
  virtual inline expression *simplify ()
  {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.hpp)
  virtual inline int Linearity ()
  {return LINEAR;}

  /// is this expression defined as an integer?
  virtual inline bool isDefinedInteger ()
  {return false;}

  /// is this variable integer?
  virtual inline bool isInteger () {

    // we definitely want to be very conservative here. A continuous
    // variable is integer if and only if its current value is really
    // integer (no round-off) -- see globallib/bearing and term x^(-3.55)

    CouNumber lb = domain_ -> lb (varIndex_);
    return (lb == domain_ -> ub (varIndex_)) && (COUENNE_round (lb) == lb);
  }

  /// Get expressions of lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Get value of lower and upper bound of an expression (if any)
  virtual void getBounds (CouNumber &lb, CouNumber &ub);

  /// Get values of lower and upper bound of an expression (if any)
  //virtual inline void getBounds (CouNumber &lb, CouNumber &ub) {
  //lb = domain_ -> lb (varIndex_);
  //ub = domain_ -> ub (varIndex_);
  //}

  /// generate cuts for expression associated with this auxiliary
  virtual void generateCuts (//const OsiSolverInterface &,
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY,
			     CouNumber =  COUENNE_INFINITY) {}

  /// generate convexification cut for constraint w = this
  virtual void generateCuts (expression *w,
			     //const OsiSolverInterface &si,
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY,
			     CouNumber =  COUENNE_INFINITY);

  /// code for comparison
  virtual inline enum expr_type code ()
  {return COU_EXPRVAR;}

  /// implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *, enum auxSign = expression::AUX_EQ);

  /// rank of an original variable is always one
  virtual inline int rank ()
  {return 1;}

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);

  /// is this variable fixed?
  virtual inline bool isFixed ()
  {return (fabs (lb () - ub ()) < COUENNE_EPS);}

  /// link this variable to a domain
  virtual inline void linkDomain (Domain *d)
  {domain_ = d;}

  /// return pointer to variable domain
  virtual inline Domain *domain ()
  {return domain_;}

  // empty for compatibility
  virtual inline void decreaseMult () {}

  /// Disable variable (empty for compatibility with exprAux)
  virtual inline void zeroMult () {}

  /// Set this variable as integer (empty for compatibility with exprAux)
  virtual inline void setInteger (bool value) {}

  /// either CONVEX, CONCAVE, AFFINE, or NONCONVEX
  virtual enum convexity convexity () const
  {return AFFINE;}

  /// return proper object to handle expression associated with this
  /// variable (NULL if this is not an auxiliary)
  virtual CouenneObject *properObject (CouenneCutGenerator *c,
				       CouenneProblem *p,
				       Bonmin::BabSetupBase *base,
				       JnlstPtr jnlst_);

  /// return its sign in the definition constraint
  virtual inline enum auxSign sign () const
  {return expression::AUX_UNDEF;}
};

}

#endif
