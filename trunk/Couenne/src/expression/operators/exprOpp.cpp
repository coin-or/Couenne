/* $Id$
 *
 * Name:    exprOpp.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprOpp.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

// find bounds of -x given bounds on x
void exprOpp::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  lb = new exprOpp (uba);
  ub = new exprOpp (lba);
}


// find bounds of -x given bounds on x
void exprOpp::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lba, uba;
  argument_ -> getBounds (lba, uba);

  lb = -uba;
  ub = -lba;
}


// differentiation
inline expression *exprOpp::differentiate (int index) 
{return new exprOpp (argument_ -> differentiate (index));}


/// implied bound processing for expression w = -x, upon change in
/// lower- and/or upper bound of w, whose index is wind
bool exprOpp::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  int ind = argument_ -> Index ();

  bool 
    res    = false, 
    argInt = argument_ -> isInteger ();

  CouNumber 
    wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
    wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  if (updateBound (-1, l + ind, argInt ? ceil  (- wu - COUENNE_EPS) : - wu)) {
    res = true; 
    chg [ind].setLower(t_chg_bounds::CHANGED);
  }

  if (updateBound ( 1, u + ind, argInt ? floor (- wl + COUENNE_EPS) : - wl)) {
    res = true; 
    chg [ind].setUpper(t_chg_bounds::CHANGED);
  }

  return res;
}


/// simplification

expression *exprOpp::simplify () {

  expression *subst = exprUnary::simplify (); // simplify what's inside first

  if (subst)
    return subst;

  // check if this is a -(-f(x))
  if (argument_ -> code () == COU_EXPROPP) {
    // leak. don't clone, use exprClone
    expression *ret = argument_ -> Argument () -> clone ();
    delete argument_;
    argument_ = NULL;
    return ret;
  }

  // check if this is a -(const)
  if (argument_ -> Type () == CONST) {
    expression *ret = new exprConst (- argument_ -> Value ());
    delete argument_;
    argument_ = NULL;
    return ret;
  }

  return NULL;
}

// print 
void exprOpp::print (std::ostream &out, 
		       bool descend) const {

  //if (printPos () == PRE)  out << printOp ();
  out << "(-"; 
  argument_ -> print (out, descend); 
  out << ")";
  //if (printPos () == POST) out << printOp ();
}
