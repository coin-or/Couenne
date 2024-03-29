/*
 *
 * Name:    exprUnary.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class for univariate functions
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_EXPRUNARY_HPP
#define COUENNE_EXPRUNARY_HPP

#include <iostream>

#include "CouenneExpression.hpp"
#include "CouenneTypes.hpp"

namespace Couenne {

/// zero function (used by default by exprUnary)
inline CouNumber zero_fun (CouNumber x)
{return 0.;}


/// expression class for unary functions (sin, log, etc.)
///
/// univariate operator-type expression: requires single argument. All
/// unary functions are derived from this base class, which has a lot
/// of common methods that need not be re-implemented by any
/// univariate class.

class COUENNELIB_EXPORT exprUnary: public expression {

 protected:

  /// single argument taken by this expression
  expression *argument_;

 public:

  /// node type
  virtual inline enum nodeType Type () const
  {return UNARY;}

  /// Constructor
  exprUnary  (expression *argument):
    argument_ (argument)        //< non-leaf expression, with argument list
  {}

  /// the operator itself (e.g. sin, log...)
  virtual inline unary_function F ()
  {return zero_fun;}

  /// Destructor
  virtual ~exprUnary ()
  {if (argument_) delete argument_;}

  /// return number of arguments
  inline int nArgs () const
  {return 1;}

  /// return argument
  virtual inline expression *Argument () const
  {return argument_;}

  /// return pointer to argument
  virtual inline expression **ArgPtr ()
  {return &argument_;}

  /// print this expression to iostream
  virtual void print (std::ostream &out = std::cout, bool = false) const;

  /// print position (PRE, INSIDE, POST)
  virtual enum pos printPos () const
  {return PRE;}

  /// print operator
  virtual std::string printOp () const
  {return "?";}

  /// compute value of unary operator
  virtual inline CouNumber operator () ()
  {return (F ()) ((*argument_) ());}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual inline int DepList (std::set <int> &deplist, enum dig_type type = ORIG_ONLY)
  {return argument_ -> DepList (deplist, type);}

  /// simplification
  expression *simplify ();

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  /// for general univariate functions, return nonlinear.
  virtual inline int Linearity ()
  {return NONLINEAR;}

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *, bool addAux = true);

  /// type of operator
  virtual inline enum expr_type code ()
  {return COU_EXPRUNARY;}

  /// is this expression integer?
  virtual bool isInteger ();

  /// compare two unary functions
  virtual int compare (exprUnary &);

  /// used in rank-based branching variable choice
  virtual inline int rank ()
  {return (argument_ -> rank ());}

  /// fill in dependence structure
  virtual inline void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g)
  {argument_ -> fillDepSet (dep, g);}

  /// replace variable with other
  virtual void replace (exprVar *, exprVar *);

  /// empty function to redirect variables to proper variable vector
  virtual inline void realign (const CouenneProblem *p)
  {argument_ -> realign (p);}
};

}

#endif
