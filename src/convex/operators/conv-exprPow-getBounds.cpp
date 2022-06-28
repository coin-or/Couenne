/*
 *
 * Name:    conv-exprPow-getBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: method to get lower and upper bounds of a power x^y
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <math.h>

#include "CouenneTypes.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMax.hpp"
#include "CouenneExprMin.hpp"
#include "CouenneExprOpp.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

// compute expressions for lower and upper bounds of a power x^y,
// based on the lower/upper bounds of x and y

void exprPow::getBounds (expression *&lb, expression *&ub) {

  // We have a standardized expression of the form w = x^y, where x or
  // y could be constant. Let us study each case separately.

  assert (arglist_ [0] -> Type () != CONST);

  // x is not constant, so it has (possibly different) lower and
  // upper bounds. The expression is x^b, b constant (the case x^y
  // has been decomposed by simplify() into exp(y log x).

  expression *lbbase, *ubbase;
  arglist_ [0] -> getBounds (lbbase, ubbase);

  //    printf ("ubbase = "); ubbase -> print (std::cout); printf ("\n");

  if (arglist_ [1] -> Type () == CONST) {

    // expression = x^b, b!=0. There are four cases:
    //
    // 1) b   is integer and odd or signpower (cube, x^5, etc)
    // 2) b   is integer and even (square, x^8, etc)
    // 3) 1/b is integer and odd  (cube root, x^(1/7), etc)
    // 4) 1/b is integer and even (square root, x^(1/4), etc)
    // 5) none of the above
    //
    // For all of these, need to check if the exponent is negative...

    CouNumber expon = arglist_ [1] -> Value ();
    int rndexp;

    bool isInt =  fabs (expon - (rndexp = COUENNE_round (expon))) < COUENNE_EPS,
      isInvInt = !isInt &&
      ((fabs (expon) > COUENNE_EPS) &&
       (fabs (1/expon - (rndexp = COUENNE_round (1/expon))) < COUENNE_EPS));

    if (issignpower_ || ((isInt || isInvInt) && (rndexp % 2) && (rndexp > 0))) {

      // the exponent is integer (or inverse integer), odd or signpower and
      // positive, hence the function is monotone non decreasing

      lb = new exprPow (lbbase, new exprConst (expon), issignpower_);
      ub = new exprPow (ubbase, new exprConst (expon), issignpower_);
    }
    else {

      // the exponent is either negative, integer even and not signpower, or fractional
      assert(!issignpower_);

      expression **all = new expression * [6];

      all [0] = new exprOpp   (lbbase);
      all [2] = new exprConst (0.);
      all [4] = ubbase;

      if (expon > 0)
	all    [1] = new exprPow (new exprClone (lbbase), new exprConst (expon));
      else all [1] = new exprPow (new exprClone (ubbase), new exprConst (expon));

      // all [3] is lower bound when lbbase <= 0 <= ubbase

      if (expon > COUENNE_EPS) all [3] = new exprConst (0.);
      else if (isInt || isInvInt) {
	if (rndexp % 2)
	  all [3] = new exprConst (-COUENNE_INFINITY);
	else all [3] = new exprMin (new exprClone (all [1]),
				    new exprPow (new exprClone (lbbase),
						 new exprConst (expon)));
      }
      else all [3] = new exprClone (all [1]);

      // all [5] is the lower bound value when lbbase <= ubbase <= 0

      if (expon > COUENNE_EPS) {
	if (isInt && !(rndexp % 2))
	  all [5] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	else all [5] = new exprConst (0.);
      }
      else {
	if (isInt || isInvInt) {
	  if (rndexp % 2)
	    all    [5] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	  else all [5] = new exprPow (new exprClone (lbbase), new exprConst (expon));
	}
	else all [5] = new exprConst (0.);
      }

      lb = new exprMin (all, 6);

      // And now the upper bound ///////////////////////////////////

      if (expon > 0) {

	// special case: upper bound depends to variable bounds only:
	// $max {lb^k, ub^k}$

	ub = new exprMax (new exprPow (new exprClone (lbbase), new exprConst (expon)),
			  new exprPow (new exprClone (ubbase), new exprConst (expon)));

      } else { // from this point on, expon < 0

	expression **alu = new expression * [6];

	alu [0] = new exprClone (all [0]);
	alu [2] = new exprConst (0.);
	alu [4] = new exprClone (ubbase);

	//if ((isInt || isInvInt) && !(rndexp % 2))
	//alu    [1] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	//else

	// if negative exponent and base has nonnegative lower bound,
	// the upper bound can only be lb^k
	alu [1] = new exprPow (new exprClone (lbbase), new exprConst (expon));

	// alu [3] is upper bound when lbbase <= 0 <= ubbase

	//if (expon < - COUENNE_EPS)
	alu [3] = new exprConst (COUENNE_INFINITY);
	//else if (isInt && !(rndexp % 2))
	//alu [3] = new exprPow (new exprMax (new exprClone (lbbase), new exprClone (ubbase)),
	//new exprConst (expon));
	//else alu [3] = new exprPow (new exprClone (ubbase), new exprConst (expon));

	// alu [5] is the upper bound value when lbbase <= ubbase <= 0

	/*if (expon > COUENNE_EPS) {

	  if (isInt && !(rndexp % 2))
	    alu [5] = new exprPow (new exprClone(ubbase), new exprConst(expon));
	  else alu [5] = new exprConst (0.);
	}
	else {*/
	if (isInt || isInvInt)
	  alu [5] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	else alu [5] = new exprConst (COUENNE_INFINITY);
	  //}

	ub = new exprMin (alu, 6);
      }
    }
  }
  else // should NOT happen, exponent is not constant...
    printf ("exprPow::getBounds(): Warning, exponent not constant\n");

  /*CouNumber l, u;
  arglist_ [0] -> getBounds (l,u);

  printf ("pow::bound: [");
  lb -> print (); printf ("=%g, ", (*lb) ());
  ub -> print (); printf ("=%g [%g,%g]\n", (*ub) (), l, u);*/
}


// get value of lower and upper bound for the expression
void exprPow::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lba, uba, k = (*(arglist_ [1])) ();
  arglist_ [0] -> getBounds (lba, uba);
  int intk;

  bool
    isInt    =           fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
    isInvInt = !isInt && fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS;

  if (!isInt && (!isInvInt || !(intk % 2 || issignpower_))) {

    // if exponent is fractional or 1/even, the base better be nonnegative

    if (lba < 0.) lba = 0.;
    if (uba < 0.) uba = 0.;
  }

  if (isInt && !(intk % 2 || issignpower_) && (k > 0)) { // x^{2h}

    if (uba < 0) {
      lb = safe_pow (-uba, k);
      ub = safe_pow (-lba, k);
    } else if (lba > 0) {
      lb = safe_pow (lba, k);
      ub = safe_pow (uba, k);
    } else {
      lb = 0;
      ub = safe_pow (CoinMax (-lba, uba), k);
    }

  } else if (k > 0) { // monotone increasing: x^{2h+1} with h integer, or x^h with h real, or signpower

    lb = safe_pow (lba, k, issignpower_);
    ub = safe_pow (uba, k, issignpower_);

  } else if (isInt && !(intk % 2 || issignpower_)) { // x^{-2h} or x^{-1/2h} with h integer

    if (uba < 0) {
      lb = safe_pow (-lba, k);
      ub = safe_pow (-uba, k);
    } else if (lba > 0) {
      lb = safe_pow (uba, k);
      ub = safe_pow (lba, k);
    } else {
      lb = safe_pow (CoinMax (-lba, uba), k);
      ub = COUENNE_INFINITY;
    }

  } else { // x^k, k<0
    assert(!issignpower_);
    if (uba < 0) {
      lb = safe_pow (uba, k);
      ub = safe_pow (lba, k);
    } else if (lba > 0) {
      lb = safe_pow (uba, k);
      ub = safe_pow (lba, k);
    } else {
      lb = -COIN_DBL_MAX; // !!! may not be reached
      ub =  COIN_DBL_MAX;
    }
  }
}
