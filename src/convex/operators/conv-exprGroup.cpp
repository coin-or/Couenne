/*
 *
 * Name:    conv-exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of convexification methods for exprGroup
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"

#include "CouenneCutGenerator.hpp"

#include "CouenneExprGroup.hpp"
#include "CouenneExprBound.hpp"
#include "CouenneExprMul.hpp"

#include "CouenneProblem.hpp"

using namespace Couenne;

/// Get lower and upper bound of an expression (if any)
void exprGroup::getBounds (expression *&lb, expression *&ub) {

  expression
    **lbnl = new expression * [1],
    **ubnl = new expression * [1];

  // TODO: do not aggregate members of exprSum

  // compute lower/upper bound of nonlinear part
  exprSum::getBounds (*lbnl, *ubnl);

  lincoeff
    *coeL = new lincoeff,
    *coeU = new lincoeff;

  // derive linear part (obtain constant)
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    std::pair <exprVar *, CouNumber> pairL, pairU;

    CouNumber coeff = el -> second;
    int         ind = el -> first -> Index ();

    pairL .second = pairU .second = coeff;

    if (coeff < 0.) {
      pairL.first = new exprUpperBound (ind, el -> first -> domain ());
      pairU.first = new exprLowerBound (ind, el -> first -> domain ());
    } else {
      pairL.first = new exprLowerBound (ind, el -> first -> domain ());
      pairU.first = new exprUpperBound (ind, el -> first -> domain ());
    }

    coeL -> push_back (pairL);
    coeU -> push_back (pairU);
  }

  lb = new exprGroup (c0_, *coeL, lbnl, 1);
  ub = new exprGroup (c0_, *coeU, ubnl, 1);

  delete coeL;
  delete coeU;
}


// /// Get expressions of lower and upper bound of an expression (if any)
// void exprGroup::getBounds (expression *&lb, expression *&ub) {

//   expression *lbnl, *ubnl;

//   // TODO: do not aggregate members of exprSum

//   // compute lower/upper bound of nonlinear part
//   exprSum::getBounds (lbnl, ubnl);

//   // count linear and constant terms
//   int nlin = lcoeff_.size();
//   if (fabs (c0_) > COUENNE_EPS) nlin++;
//   //  for (int *ind = index_; *ind++>=0; nlin++);

//   expression
//     **linall = new expression * [nlin + 1], // linear arglist for lower bound
//     **linalu = new expression * [nlin + 1]; //                    upper

//   // add constant to bounds
//   if (fabs (c0_) > COUENNE_EPS) {
//     *linall++ = new exprConst (c0_);
//     *linalu++ = new exprConst (c0_);
//   }

//   // TODO: make it another exprGroup!

//   // derive linear part (obtain constant)
//   for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

//     //    c0 += el -> second;

//     CouNumber coeff = el -> second;
//     int         ind = el -> first -> Index ();

//     expression *l = new exprLowerBound (ind, el -> first -> domain ()),
//                *u = new exprUpperBound (ind, el -> first -> domain ());

//     if (fabs (coeff - 1.) < COUENNE_EPS) {
//       *linall++ = l;
//       *linalu++ = u;
//     } else {

//       expression *c1 = new exprConst (coeff),
//                  *c2 = new exprConst (coeff);

//       if (coeff < 0) {
// 	*linall++ = new exprMul (c1, u);
// 	*linalu++ = new exprMul (c2, l);
//       } else {
// 	*linall++ = new exprMul (c1, l);
// 	*linalu++ = new exprMul (c2, u);
//       }
//     }
//   }

//   *linall = lbnl;
//   *linalu = ubnl;

//   lb = new exprSum (linall - nlin, nlin + 1);
//   ub = new exprSum (linalu - nlin, nlin + 1);
// }



/// Get values of  lower and upper bound of an expression (if any)
void exprGroup::getBounds (CouNumber &lb, CouNumber &ub) {

  // compute lower/upper bound of nonlinear part

  exprSum::getBounds (lb, ub);
  /*{
    expression *le, *ue;
    getBounds (le, ue);

    lb = (*le) ();
    ub = (*ue) ();

    delete le;
    delete ue;
    }*/

  lb += c0_;
  ub += c0_;

  // derive linear part (obtain constant)
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    exprVar    *var = el -> first;
    CouNumber coeff = el -> second, vlb, vub;

    bool
      inf_lb = false,
      inf_ub = false;

    if ((vlb = var -> lb ()) < -COUENNE_INFINITY) {
      if (coeff > 0) inf_lb = true;
      else           inf_ub = true;
    } else {
      if (coeff > 0) lb += vlb * coeff;
      else           ub += vlb * coeff;
    }

    if ((vub = var -> ub ()) >  COUENNE_INFINITY) {
      if (coeff > 0) inf_ub = true;
      else           inf_lb = true;
    } else {
      if (coeff > 0) ub += vub * coeff;
      else           lb += vub * coeff;
    }

    if (inf_lb)
      lb = -COUENNE_INFINITY;

    if (inf_ub) {
      ub =  COUENNE_INFINITY;
      if (inf_lb)
	break;
    }
  }
}


// generate equality between *this and *w
void exprGroup::generateCuts (expression *w,
			      OsiCuts &cs, const CouenneCutGenerator *cg,
			      t_chg_bounds *chg,
			      int wind, CouNumber lb, CouNumber ub) {

  // very similar to exprSum::generateCuts. First of all, this has
  // been standardized into a sum, so it only gets a cut in the
  // initial relaxation
  if (!(cg -> isFirst ()))
    return;

  // there is one linear term so far: -w
  int nterms = lcoeff_.size ();

  OsiRowCut *cut = new OsiRowCut;

  // If this aux is fixed, don't write
  //
  // "- w + ax = -b" but just
  //
  // "ax = -b+ w0"
  //
  // with w0 its constant value

  CouNumber vlb, vub;
  w -> getBounds (vlb, vub);
  bool uselessAux = (vub < vlb + COUENNE_EPS);

  // TODO: generalize to sign!= ::EQ

  int displacement = (wind < 0 && !uselessAux) ? 1: 0;

  CouNumber *coeff = new CouNumber [nargs_ + nterms + displacement];
  int       *index = new int       [nargs_ + nterms + displacement];

  if (wind < 0 && !uselessAux) {
    // first, make room for aux variable
    coeff [0] = -1.;
    index [0] = w -> Index ();
    lb = ub = 0.;
  }

  if (uselessAux)
    lb = ub = vlb;

  lb -= c0_;
  ub -= c0_;

  // now add linear terms
  lincoeff::iterator el = lcoeff_.begin ();

  nterms = displacement;

  for (; el != lcoeff_.end (); ++el)

    if (fabs (el -> second) > 1.e-21) {
      // why 1.0e-21? Look at CoinPackedMatrix.cpp:2237

      coeff [nterms]   = el -> second;
      index [nterms++] = el -> first -> Index ();
    }

  // scan arglist for (aux) variables and constants
  for (int i=0; i<nargs_; i++) {

    expression *curr = arglist_ [i];

    if (curr -> Type () == CONST) {// constant term in sum
      lb -= curr -> Value ();
      ub -= curr -> Value ();
    }
    else {                        // variable
      coeff [nterms]   = 1.;
      index [nterms++] = curr -> Index ();
    }
  }

  cut -> setRow (nterms, index, coeff);

  delete [] index;
  delete [] coeff;

  enum auxSign sign = expression::AUX_EQ;

  if (wind < 0)
    sign = cg -> Problem () -> Var (w -> Index ()) -> sign ();

  if (lb > -COUENNE_INFINITY && (wind >= 0 || sign != expression::AUX_GEQ)) cut -> setLb (lb);
  if (ub <  COUENNE_INFINITY && (wind >= 0 || sign != expression::AUX_LEQ)) cut -> setUb (ub);

  cut -> setGloballyValid (); // added only once, it is global
  cs.insert (cut);
  delete cut;
}
