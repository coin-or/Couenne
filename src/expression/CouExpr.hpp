/*
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

class COUENNELIB_EXPORT CouExpr {

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


COUENNELIB_EXPORT CouExpr operator+(CouExpr &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator/(CouExpr &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator%(CouExpr &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator-(CouExpr &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator*(CouExpr &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator^(CouExpr &e1, CouExpr &e2);

COUENNELIB_EXPORT CouExpr &sin(CouExpr &e);
COUENNELIB_EXPORT CouExpr &cos(CouExpr &e);
COUENNELIB_EXPORT CouExpr &log(CouExpr &e);
COUENNELIB_EXPORT CouExpr &exp(CouExpr &e);

COUENNELIB_EXPORT CouExpr &operator+(CouNumber &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator/(CouNumber &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator%(CouNumber &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator-(CouNumber &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator*(CouNumber &e1, CouExpr &e2);
COUENNELIB_EXPORT CouExpr &operator^(CouNumber &e1, CouExpr &e2);

COUENNELIB_EXPORT CouExpr &sin(CouNumber &e);
COUENNELIB_EXPORT CouExpr &cos(CouNumber &e);
COUENNELIB_EXPORT CouExpr &log(CouNumber &e);
COUENNELIB_EXPORT CouExpr &exp(CouNumber &e);

COUENNELIB_EXPORT CouExpr &operator+(CouExpr &e1, CouNumber &e2);
COUENNELIB_EXPORT CouExpr &operator/(CouExpr &e1, CouNumber &e2);
COUENNELIB_EXPORT CouExpr &operator%(CouExpr &e1, CouNumber &e2);
COUENNELIB_EXPORT CouExpr &operator-(CouExpr &e1, CouNumber &e2);
COUENNELIB_EXPORT CouExpr &operator*(CouExpr &e1, CouNumber &e2);
COUENNELIB_EXPORT CouExpr &operator^(CouExpr &e1, CouNumber &e2);

}

#endif
