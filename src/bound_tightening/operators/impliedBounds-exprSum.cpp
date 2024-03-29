/*
 *
 * Name:    impliedBounds-exprSum.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bound enforcing for exprSum and exprGroup
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprSum.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

/// vector operation to find bound to variable in a sum
static CouNumber scanBounds (int, int, int *, CouNumber *, CouNumber *, int *);


/// implied bound processing for expression w = x+y+t+..., upon change
/// in lower- and/or upper bound of w, whose index is wind

bool exprSum::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  /**
   *  An expression
   *
   *  \f$w = a0 + \sum_{i\in I1} a_i x_i + \sum_{i\in I2} a_i x_i\f$
   *
   *  is given such that all $a_i$ are positive for $i \in I1$ and
   *  negative for $i in I2$. If the bounds on $w \in [l,b]$, implied
   *  bounds on all $x_i, i\in I1 \cup I2$ are as follows:
   *
   *  \f$\forall i\in I1\f$
   *
   *    \f$x_i \ge (l - a0 - sum_{j\in I1 | j\neq i} a_j u_j - sum_{j\in I2}        a_j l_j) / a_i\f$
   *    \f$x_i \le (u - a0 - sum_{j\in I1 | j\neq i} a_j l_j - sum_{j\in I2}        a_j u_j) / a_i\f$
   *
   *  \f$\forall i\in I2\f$
   *
   *    \f$x_i \ge (u - a0 - sum_{j\in I1} a_j u_j        - sum_{j\in I2 | j\neq i} a_j l_j) / a_i\f$
   *    \f$x_i \le (l - a0 - sum_{j\in I1} a_j l_j        - sum_{j\in I2 | j\neq i} a_j u_j) / a_i\f$,
   *
   *  where \f$l_i\f$ and \f$u_i\f$ are lower and upper bound,
   *  respectively, of \f$x_i\f$. We also have to check if some of
   *  these bounds are infinite.
   */

  // compute number of pos/neg coefficients and sum up constants

  CouNumber
    a0 = 0.,   // constant term in the sum
    wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
    wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  // quick check: if both are infinite, nothing is implied...

  if ((wl < -COUENNE_INFINITY) &&
      (wu >  COUENNE_INFINITY))
    return false;

  int
    nterms = nargs_, // # nonconstant terms
    nlin   = 0;      // # linear terms

  exprGroup *eg = NULL;

  bool isExprGroup = (code () == COU_EXPRGROUP);

  if (isExprGroup) { // if exprGroup, compute no. linear terms

    eg = dynamic_cast <exprGroup *> (this);

    a0   += eg -> getc0 ();
    nlin += eg -> lcoeff (). size ();
  }

  nterms += nlin;

  // Coefficients and indices of the positive and the negative
  // non-constant parts of the sum (at most nlin are negative, as all
  // "nonlinear" terms have coefficient 1)

  CouNumber *C1 = (CouNumber *) malloc (nterms * sizeof (CouNumber)),
            *C2 = (CouNumber *) malloc (nlin   * sizeof (CouNumber));
  int       *I1 = (int       *) malloc (nterms * sizeof (int)),
            *I2 = (int       *) malloc (nlin   * sizeof (int));

  int ipos, ineg = ipos = 0; // #pos and #neg terms

  std::set <int> intSet; // contains indices of integer variables

  // classify terms as positive or constant for the exprSum

  for (int i = nargs_; i--;) {

    int index = arglist_ [i] -> Index ();

    if (index < 0) // assume a non variable is a constant
      a0 += arglist_ [i] -> Value ();
    else {
      I1 [ipos]   = index;
      C1 [ipos++] = 1.;

      // add entry to integer variable
      if (arglist_ [i] -> isDefinedInteger ())
	intSet.insert (index);
    }
  }

  // classify terms as pos/neg or constant for the remaining of exprGroup

  if (isExprGroup) {

    exprGroup::lincoeff &lcoe = eg -> lcoeff ();

    for (exprGroup::lincoeff::iterator el = lcoe.begin ();
	 el != lcoe.end (); ++el) {

      CouNumber coe = el -> second;
      int       ind = el -> first -> Index ();

      if      (coe >  0.) {I1 [ipos] = ind; C1 [ipos++] = coe;}
      else if (coe < -0.) {I2 [ineg] = ind; C2 [ineg++] = coe;}

      if (el -> first -> isDefinedInteger ())
	intSet.insert (ind);
    }
  }

  // realloc to save some memory

  /*printf ("  a0 = %g, w in [%g, %g]\n", a0, wl, wu);
  printf ("  pos: "); for (int i=0; i<ipos; i++) printf ("%g x%d [%g,%g] ", C1 [i], I1 [i], l [I1 [i]], u [I1 [i]]); printf ("\n");
  printf ("  neg: "); for (int i=0; i<ineg; i++) printf ("%g x%d [%g,%g] ", C2 [i], I2 [i], l [I2 [i]], u [I2 [i]]); printf ("\n");*/

  // now we have correct  I1, I2, C1, C2, ipos, ineg, and a0

  // indices of the variable in I1 or I2 with infinite lower or upper
  // bound. If more than one is found, it is set to -2

  int infLo1 = -1, infLo2 = -1,
      infUp1 = -1, infUp2 = -1;

  // Upper bound of the sum, considering lower/upper bounds of the
  // variables but negliging the infinite ones:
  //
  // lower = $a0 + \sum_{i in I_1: a_i <  \infinity} a_i l_i
  //             + \sum_{i in I_2: a_i > -\infinity} a_i u_i$
  //
  // upper = $a0 + \sum_{i in I_1: a_i <  \infinity} a_i u_i
  //             + \sum_{i in I_2: a_i > -\infinity} a_i l_i$

  CouNumber

    lower = a0 + scanBounds (ipos, -1, I1, C1, l, &infLo1)
               + scanBounds (ineg, +1, I2, C2, u, &infUp2),

    upper = a0 + scanBounds (ipos, +1, I1, C1, u, &infUp1)
               + scanBounds (ineg, -1, I2, C2, l, &infLo2);

  // Compute lower bound for all or for some of the variables: There
  // is a bound for all variables only if both infUp1 and infLo2 are
  // -1, otherwise there is a bound only for one variable if one is -1
  // and the other is nonnegative.

  bool tighter = false;

  // check if there is room for improvement

  {
    CouNumber
      slackL = lower - wl, // each is positive if no implication
      slackU = wu - upper; //

    // if lower < wl or upper > wu, some bounds can be tightened.
    //
    // otherwise, there is no implication, but if lower
    //
    // steal some work to bound propagation...

    if ((slackL > 0.) &&
	(infLo1 == -1) &&
	(infUp2 == -1)) {   // no implication on lower

      // propagate lower bound to w
      if ((sign != expression::AUX_LEQ) && (updateBound (-1, l + wind, lower))) {
	chg [wind].setLower(t_chg_bounds::CHANGED);
	tighter = true;
	if (intSet.find (wind)!= intSet.end ())
	  l [wind] = ceil (l [wind] - COUENNE_EPS);
      }

      if ((slackU > 0.) &&
	  (infLo2 == -1) &&
	  (infUp1 == -1)) { // no implication on upper

	// propagate upper bound to w
	if ((sign != expression::AUX_GEQ) && (updateBound (+1, u + wind, upper))) {
	  chg [wind].setUpper(t_chg_bounds::CHANGED);
	  tighter = true;
	  if (intSet.find (wind)!= intSet.end ())
	    u [wind] = floor (u [wind] + COUENNE_EPS);
	}

	free (I1); free (I2);
	free (C1); free (C2);

	return false; // both bounds were weak, no implication possible
      }
    }
    else if ((slackU > 0.) &&
	     (infLo2 == -1) &&
	     (infUp1 == -1)) {

      // propagate upper bound to w
      if ((sign != expression::AUX_GEQ) && (updateBound (+1, u + wind, upper))) {
	tighter = true;
	chg [wind].setUpper(t_chg_bounds::CHANGED);
	if (intSet.find (wind)!= intSet.end ())
	  u [wind] = floor (u [wind] + COUENNE_EPS);
      }
    }
  }

  // Subtle... make two copies of lower and upper to avoid updating
  // bounds with some previously updated bounds

  // first of all, find maximum index in I1 and I2

  int maxind = -1;

  for (int i=ipos; i--; I1++) if (*I1 > maxind) maxind = *I1;
  for (int i=ineg; i--; I2++) if (*I2 > maxind) maxind = *I2;

  I1 -= ipos;
  I2 -= ineg;

  CouNumber *lc = (CouNumber *) malloc (++maxind * sizeof (CouNumber));
  CouNumber *uc = (CouNumber *) malloc (maxind   * sizeof (CouNumber));

  CoinCopyN (l, maxind, lc);
  CoinCopyN (u, maxind, uc);

  // Update lowers in I1 and uppers in I2

  if ((infLo1 == -1) && (infUp2 == -1) && (wu < COUENNE_INFINITY / 1e10)) {
    // All finite bounds. All var. bounds can be tightened.

    // tighten upper bound of variables in I1
    for (int i=ipos; i--;) {
      int ind = I1 [i];
      if ((tighter = (updateBound (+1, u + ind, (wu - lower) / C1 [i] + lc [ind]) || tighter))) {
	chg [ind].setUpper(t_chg_bounds::CHANGED);
	if (intSet.find (ind)!= intSet.end ())
	  u [ind] = floor (u [ind] + COUENNE_EPS);
      }
    }

    // tighten lower bound of variables in I2
    for (int i=ineg; i--;) {
      int ind = I2 [i];
      if ((tighter = (updateBound (-1, l + ind, (wu - lower) / C2 [i] + uc [ind]) || tighter))) {
	chg [ind].setLower(t_chg_bounds::CHANGED);
	if (intSet.find (ind)!= intSet.end ())
	  l [ind] = ceil (l [ind] - COUENNE_EPS);
      }
    }
  } else

    if ((infLo1 >= 0) && (infUp2 == -1)) {    // There is one infinite lower bound in I1
      int ind = I1 [infLo1];
      if ((tighter = (updateBound (+1, u + ind, (wu - lower) / C1 [infLo1]) || tighter))) {
	chg [ind].setUpper(t_chg_bounds::CHANGED);
	if (intSet.find (ind)!= intSet.end ())
	  u [ind] = floor (u [ind] + COUENNE_EPS);
      }
    }
    else
      if ((infLo1 == -1) && (infUp2 >= 0)) {  // There is one infinite upper bound in I2
	int ind = I2 [infUp2];
	if ((tighter = (updateBound (-1, l + ind, (wu - lower) / C2 [infUp2]) || tighter))) {
	  chg [ind].setLower(t_chg_bounds::CHANGED);
	  if (intSet.find (ind)!= intSet.end ())
	    l [ind] = ceil (l [ind] - COUENNE_EPS);
	}
      }

  // Update uppers in I1 and lowers in I2

  if ((infUp1 == -1) && (infLo2 == -1) && (wl > -COUENNE_INFINITY / 1e10)) {
    // All finite bounds. All var. bounds can be tightened.

    for (int i=ipos; i--;) {
      int ind = I1 [i];
      if ((tighter = (updateBound (-1, l + ind, (wl - upper) / C1 [i] + uc [ind]) || tighter))) {
	chg [ind].setLower(t_chg_bounds::CHANGED);
	if (intSet.find (ind) != intSet.end ())
	  l [ind] = ceil (l [ind] - COUENNE_EPS);
      }
    }

    for (int i=ineg; i--;) {
      int ind = I2 [i];
      if ((tighter = (updateBound (+1, u + ind, (wl - upper) / C2 [i] + lc [ind]) || tighter))) {
	chg [ind].setUpper(t_chg_bounds::CHANGED);
	if (intSet.find (ind) != intSet.end ())
	  u [ind] = floor (u [ind] + COUENNE_EPS);
      }
    }
  } else

    if ((infUp1 >= 0) && (infLo2 == -1)) { // There is one infinite lower bound in I2
      int ind = I1 [infUp1];
      if ((tighter = (updateBound (-1, l + ind, (wl - upper) / C1 [infUp1]) || tighter))) {
	chg [ind].setLower(t_chg_bounds::CHANGED);
	if (intSet.find (ind) != intSet.end ())
	  l [ind] = ceil (l [ind] - COUENNE_EPS);
      }
    }
    else
      if ((infUp1 == -1) && (infLo2 >= 0)) {  // There is one infinite upper bound in I1
	int ind = I2 [infLo2];
	if ((tighter = (updateBound (+1, u + ind, (wl - upper) / C2 [infLo2]) || tighter))) {
	  chg [ind].setUpper(t_chg_bounds::CHANGED);
	  if (intSet.find (ind) != intSet.end ())
	    u [ind] = floor (u [ind] + COUENNE_EPS);
	}
      }

  // ...phew!

  free (I1); free (I2);
  free (C1); free (C2);
  free (lc); free (uc);

  return tighter;
}


/// sum bounds depending on coefficients

static CouNumber scanBounds (int        num,      /// cardinality of the set (I1 or I2)
			     int        sign,     /// +1: check for +inf, -1: check for -inf
			     int       *indices,  /// indices in the set, $i\in I_k$
			     CouNumber *coeff,    /// coefficients in the set
			     CouNumber *bounds,   /// variable bounds to check
			     int       *infnum) { /// return value: index of variable with infinite
                                                  /// bound, or -2 if more than one exist
  CouNumber bound = 0.;

  for (int i = num; i--;) {

    CouNumber bd = bounds [indices [i]];

    // be sensitive here, check for bounds a little within the finite realm

    if (((sign > 0) ? bd : -bd) > COUENNE_INFINITY / 1e10) {

      bounds [indices [i]] = (sign > 0) ? COIN_DBL_MAX : -COIN_DBL_MAX;
      //bounds [indices [i]] = (sign > 0) ? 1e300 : -1e300;

      // this variable has an infinite bound, mark it
      if      (*infnum == -1) *infnum =  i; // first variable with infinite bound, so far
      else if (*infnum >=  0) *infnum = -2; // oops... more than one found, no finite bound
    }
    else bound += coeff [i] * bd; // no infinity detected, sum a weighted, finite bound
  }

  return bound;
}
