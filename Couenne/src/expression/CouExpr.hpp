/* $Id$
 *
 * Name:    CouExpr.hpp
 * Author:  Pietro Belotti
 * Purpose: Container class for expressions
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouExpr_hpp
#define CouExpr_hpp

#include "CouenneExpression.hpp"

namespace Couenne {

class CouExpr {

private:

  expression *expr_;

public:

  CouExpr (expression *e):
    expr_ (e) {}

  CouExpr (const CouExpr &e) {
    expr_ = e.expr_ -> clone ();
  }

  CouExpr &operator=(CouExpr &e) {
    expr_ = e.expr_ -> clone ();
    return *this;
  }

  expression *Expression () const 
  {return expr_;}
};


CouExpr operator+(CouExpr &e1, CouExpr &e2);
CouExpr &operator/(CouExpr &e1, CouExpr &e2);
CouExpr &operator%(CouExpr &e1, CouExpr &e2);
CouExpr &operator-(CouExpr &e1, CouExpr &e2);
CouExpr &operator*(CouExpr &e1, CouExpr &e2);
CouExpr &operator^(CouExpr &e1, CouExpr &e2);

CouExpr &sin(CouExpr &e);
CouExpr &cos(CouExpr &e);
CouExpr &log(CouExpr &e);
CouExpr &exp(CouExpr &e);

CouExpr &operator+(CouNumber &e1, CouExpr &e2);
CouExpr &operator/(CouNumber &e1, CouExpr &e2);
CouExpr &operator%(CouNumber &e1, CouExpr &e2);
CouExpr &operator-(CouNumber &e1, CouExpr &e2);
CouExpr &operator*(CouNumber &e1, CouExpr &e2);
CouExpr &operator^(CouNumber &e1, CouExpr &e2);

CouExpr &sin(CouNumber &e);
CouExpr &cos(CouNumber &e);
CouExpr &log(CouNumber &e);
CouExpr &exp(CouNumber &e);

CouExpr &operator+(CouExpr &e1, CouNumber &e2);
CouExpr &operator/(CouExpr &e1, CouNumber &e2);
CouExpr &operator%(CouExpr &e1, CouNumber &e2);
CouExpr &operator-(CouExpr &e1, CouNumber &e2);
CouExpr &operator*(CouExpr &e1, CouNumber &e2);
CouExpr &operator^(CouExpr &e1, CouNumber &e2);

}

#endif
