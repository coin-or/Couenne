/*
 *
 * Name:    createCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: a standard cut creator for use with convexification
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiRowCut.hpp"

#include "CouenneTypes.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

using namespace Ipopt;
using namespace Couenne;

/// checks if very large or very small nonzero
bool badCoeff (CouNumber coe) {

  coe = fabs (coe);
  return ((coe > COU_MAX_COEFF) || ((coe < COU_MIN_COEFF) && (coe > 0.)));
}


/// general procedure for inserting a linear cut with up to three
/// variables. Return 1 if cut inserted, 0 if none, <0 if error

int CouenneCutGenerator::createCut (OsiCuts &cs,
				    CouNumber lb, CouNumber ub,
				    int i1, CouNumber c1,
				    int i2, CouNumber c2,
				    int i3, CouNumber c3,
				    bool is_global) const {
  bool numerics = false;

  // a maximum of three terms allowed here. If index -1 (default) the
  // term is not considered

  int nterms = 0;

  // first check: duplicate index.

  if ((i1 != -1) && (i1 == i3)) {i3 = -1; c1 += c3; c3 = 0;}
  if ((i2 != -1) && (i2 == i3)) {i3 = -1; c2 += c3; c3 = 0;}
  if ((i1 != -1) && (i1 == i2)) {i2 = -1; c1 += c2; c2 = 0;}

  // CAUTION: this can make the problem infeasible...
  if (fabs (c3) <= 1.e-21) {                                    i3 = -1;} // shift coeff/index to
  if (fabs (c2) <= 1.e-21) {                  c2 = c3; i2 = i3; i3 = -1;} // keep consistency
  if (fabs (c1) <= 1.e-21) {c1 = c2; i1 = i2; c2 = c3; i2 = i3; i3 = -1;}
  // why 1.e-21? Look at CoinPackedMatrix.cpp:2273

#if 0
  if (i1 >= 0) {if (fabs (c1) > COU_MAX_COEFF) numerics = true; nterms++;} else c1 = 0;
  if (i2 >= 0) {if (fabs (c2) > COU_MAX_COEFF) numerics = true; nterms++;} else c2 = 0;
  if (i3 >= 0) {if (fabs (c3) > COU_MAX_COEFF) numerics = true; nterms++;} else c3 = 0;
#else
  if (i1 >= 0) {if (badCoeff (c1)) numerics = true; nterms++;} else c1 = 0;
  if (i2 >= 0) {if (badCoeff (c2)) numerics = true; nterms++;} else c2 = 0;
  if (i3 >= 0) {if (badCoeff (c3)) numerics = true; nterms++;} else c3 = 0;
#endif

  if (!nterms) // nonsense cut
    return 0;

  // cut has large coefficients/rhs, bail out
  if (numerics
      //|| ((fabs (lb) < COU_MIN_COEFF) ||
      //(fabs (ub) < COU_MIN_COEFF))
      || ((fabs (lb) > COU_MAX_COEFF) &&
	  (fabs (ub) > COU_MAX_COEFF))) {

    jnlst_->Printf(J_STRONGWARNING, J_CONVEXIFYING,
		   "### Discarding cut, large coeff/rhs: %g (%d), %g (%d), %g (%d); [%g,%g]\n",
		   c1, i1, c2, i2, c3, i3, lb, ub);
    return 0;
  }

  if (!firstcall_ && addviolated_) { // need to check violation

    const CouNumber *x = problem_ -> X ();

    // compute violation
    CouNumber violation = 0.;

    if (i1 >= 0) violation += c1 * x [i1];
    if (i2 >= 0) violation += c2 * x [i2];
    if (i3 >= 0) violation += c3 * x [i3];

    // quit if not violated

    if ((violation < ub + 0 * COUENNE_EPS) &&
	(violation > lb - 0 * COUENNE_EPS))
      return 0;
  }

  // check if cut violates optimal solution (irrespective of the
  // branching rules applied, so handle with care)

  CouNumber *best = problem_ -> bestSol ();

  if (best &&
      ((i1 < 0) || ((best [i1] >= problem_ -> Lb (i1)) && (best [i1] <= problem_ -> Ub (i1)))) &&
      ((i2 < 0) || ((best [i2] >= problem_ -> Lb (i2)) && (best [i2] <= problem_ -> Ub (i2)))) &&
      ((i3 < 0) || ((best [i3] >= problem_ -> Lb (i3)) && (best [i3] <= problem_ -> Ub (i3))))) {

    CouNumber lhs = 0.;

    if (i1 >= 0) lhs += c1 * best [i1];
    if (i2 >= 0) lhs += c2 * best [i2];
    if (i3 >= 0) lhs += c3 * best [i3];

    if (lhs > ub + COUENNE_EPS)
      jnlst_->Printf(J_STRONGWARNING, J_CONVEXIFYING,
		     "### cut (%d,%d,%d) (%g,%g,%g) cuts optimum: %g >= %g [%g]\n",
		     i1,i2,i3, c1,c2,c3, lhs, ub, lhs - ub);

    if (lhs < lb - COUENNE_EPS)
      jnlst_->Printf(J_STRONGWARNING, J_CONVEXIFYING,
		     "### cut (%d,%d,%d) (%g,%g,%g) cuts optimum: %g <= %g [%g]\n",
		     i1,i2,i3, c1,c2,c3, lhs, lb, lb - lhs);
  }

  // You are here if:
  //
  // 1) this is the first call to CouenneCutGenerator::generateCuts(), or
  // 2) you also want unviolated cuts, or
  // 3) the cut is violated

  // two cases: cut is of the form w1 [<|>]= alpha, hence a column
  // cut, or it is of the form (a w1 + b w2 + c w3 [<|>]= alpha), a
  // row cut

  if ((i2 < 0) && (i3 < 0)) { // column cut /////////////////////////////////////////

    if (   (fabs (c1) < COUENNE_EPS)
	&& (fabs (lb) > COU_MAX_COEFF * COUENNE_EPS)
	&& (fabs (ub) > COU_MAX_COEFF * COUENNE_EPS)) {

      jnlst_->Printf(J_STRONGWARNING, J_CONVEXIFYING,
		     "#### nonsense column cut: %e <= %e w_%d <= %e\n",
		     lb, c1, i1, ub);
      return 0;
    }

    OsiColCut *cut = new OsiColCut;

    CouNumber
      ll = lb / c1,
      uu = ub / c1;

    if (c1 < 0) {
      CouNumber tmp = ll;
      ll = uu;
      uu = tmp;
    }

    CouNumber &curL = problem_ -> Lb (i1),
              &curU = problem_ -> Ub (i1);

    if ((uu < COUENNE_INFINITY) &&
	(uu < curU - COUENNE_EPS)) {

      cut -> setUbs (1, &i1, &uu);
      curU = uu; // TODO: chg_bds
    }

    if ((ll > -COUENNE_INFINITY) &&
	(ll > curL + COUENNE_EPS)) {
      cut -> setLbs (1, &i1, &ll);
      curL = ll; // idem
    }

    cut -> setGloballyValid (is_global); // global?

    cs.insert (cut);
    delete cut;

  } else {

    // row cut //////////////////////////////////////////////////////////////////////

    CouNumber *coeff = new CouNumber [nterms];
    int       *index = new int       [nterms];
    OsiRowCut *cut   = new OsiRowCut;

    int nt = 0;

    if (i1 >= 0) {coeff [nt] = c1; index [nt++] = i1;}
    if (i2 >= 0) {coeff [nt] = c2; index [nt++] = i2;}
    if (i3 >= 0) {coeff [nt] = c3; index [nt++] = i3;}

    if (lb > -COUENNE_INFINITY) cut -> setLb (lb);
    if (ub <  COUENNE_INFINITY) cut -> setUb (ub);

    cut -> setRow (nterms, index, coeff);

    delete [] coeff;
    delete [] index;

    cut -> setGloballyValid (is_global); // global?

    cs.insert (cut);
    delete cut;
  }

  return 1;
}


/// general procedure for inserting a linear cut with up to three
/// variables. Return 1 if cut inserted, 0 if none, <0 if error

int CouenneCutGenerator::createCut (OsiCuts &cs,
				    CouNumber rhs, int sign,
				    int i1, CouNumber c1,
				    int i2, CouNumber c2,
				    int i3, CouNumber c3,
				    bool is_global)       const {

  return createCut (cs, (sign >= 0) ? rhs : - COIN_DBL_MAX,
		        (sign <= 0) ? rhs :   COIN_DBL_MAX,
		    i1, c1, i2, c2, i3, c3, is_global);
}
