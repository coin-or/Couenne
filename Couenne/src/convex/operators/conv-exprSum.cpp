/* $Id$
 *
 * Name:    conv-exprSum.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to standardize/convexify sum expressions
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExprAux.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprConst.hpp"

#include "CouenneProblem.hpp"

using namespace Couenne;

// generate equality between *this and *w
void exprSum::generateCuts (expression *w, //const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, 
			    int wind, CouNumber lb, CouNumber ub) {
  if (!(cg -> isFirst ()))
    return;

  CouNumber *coeff = new CouNumber [nargs_ + 1];
  int       *index = new int       [nargs_ + 1];
  OsiRowCut *cut   = new OsiRowCut;

  /// first, make room for aux variable

  int nv = 0;

  bool uselessAux = (ub < lb + COUENNE_EPS); 

  if (wind < 0 && !uselessAux) {
    coeff [0] = -1.; index [0] = w -> Index ();
    nv++;
    lb = ub = 0;
  }

  /// scan arglist for (aux) variables and constants
  for (int i=0; i<nargs_; i++) {

    if (arglist_ [i] -> Type () == CONST) {

      CouNumber val = arglist_ [i] -> Value ();

      lb -= val;
      ub -= val;
    }
    else {
      coeff [nv]   = 1.; 
      index [nv++] = arglist_ [i] -> Index ();
    }
  }

  cut -> setRow (nv, index, coeff);

  delete [] index;
  delete [] coeff;

  enum auxSign sign = cg -> Problem () -> Var (w -> Index ()) -> sign ();

  if (lb > -COUENNE_INFINITY && (sign != expression::AUX_GEQ)) cut -> setLb (lb);
  if (ub <  COUENNE_INFINITY && (sign != expression::AUX_LEQ)) cut -> setUb (ub);

  /// added only once, it is global
  cut -> setGloballyValid ();

  cs.insert (cut);
  delete cut;
}
