/* */
/*
 * Name:    exprMinMax.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of min and max operators
 *
 * (C) Carnegie-Mellon University, 2006.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiSolverInterface.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneTypes.hpp"
#include "CouenneExprMax.hpp"
#include "CouenneExprMin.hpp"
#include "CouenneExprConst.hpp"

using namespace Couenne;

// Get lower and upper bound of an expression (if any)
void exprMin::getBounds (expression *&lower, expression *&upper) {
  lower = new exprConst (-COIN_DBL_MAX);
  upper = new exprConst ( COIN_DBL_MAX);
}


// Get lower and upper bound of an expression (if any)
void exprMax::getBounds (expression *&lower, expression *&upper) {
  lower = new exprConst (-COIN_DBL_MAX);
  upper = new exprConst ( COIN_DBL_MAX);
}


void exprMin::generateCuts (expression *w, //const OsiSolverInterface &si,
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber)
{}


void exprMax::generateCuts (expression *w, //const OsiSolverInterface &si,
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber)
{}
