/* $Id$
 *
 * Name:    tightenBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: bound tightening for current linear relaxation
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

using namespace Ipopt;
using namespace Couenne;

/// Bound propagation for auxiliary variables

int CouenneProblem::tightenBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. 

  // check all auxiliary variables for changes in their upper,
  // lower bound, depending on the bound changes of the variables
  // they depend on

  CouNumber *knownOptimum = optimum_;

  if (optimum_) {

    for (int i=nVars(); i--; knownOptimum++)

      if (*knownOptimum < Lb (i) || 
	  *knownOptimum > Ub (i)) {

	knownOptimum = NULL;
	break;
      }

    if (knownOptimum) 
      knownOptimum -= nVars ();
  }

  bool dbgOutput = Jnlst () -> ProduceOutput (J_DETAILED, J_BOUNDTIGHTENING);

  if (dbgOutput) {
    // ToDo: Pipe all output through journalist
    Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING, "  forward  =====================\n  ");
    int j=0;
    for (int i=0; i < nVars (); i++) 
      if (variables_ [i] -> Multiplicity () > 0) {
	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"x_%03d [%+10g %+10g] ", i, 
			domain_. lb (i),
			domain_. ub (i));
	if (!(++j % 6)) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n  ");
    }
    if (j % 6) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
  }

  for (int ii = 0, j = nVars (); j--; ii++) {

    int i = numbering_ [ii];

    exprVar *var = Var (i);

    if (var -> Multiplicity () <= 0) 
      continue;

    CouNumber &lower_i = Lb (i);
    CouNumber &upper_i = Ub (i);

    // early test to avoid a loop

    if ((fabs (upper_i) < COIN_DBL_MAX / 10) &&
	(fabs (lower_i) < COIN_DBL_MAX / 10) &&
	(lower_i > upper_i + COUENNE_BOUND_PREC * (1 + CoinMin (fabs (lower_i), fabs (upper_i))))) {
	// (upper_i < - MAX_BOUND) ||
	// (lower_i >   MAX_BOUND)) {

      if (Jnlst()->ProduceOutput(J_ITERSUMMARY, J_BOUNDTIGHTENING)) {

	Jnlst()->Printf(J_ITERSUMMARY, J_BOUNDTIGHTENING,
			"pre-check: w_%d has infeasible bounds [%.10e,%.10e]. ", i, lower_i, upper_i);

	if (dbgOutput) {

	  var -> Lb () -> print (std::cout);
	  Jnlst () -> Printf (J_DETAILED, J_BOUNDTIGHTENING," --- ");
	  var -> Ub () -> print (std::cout);
	}

	Jnlst()->Printf(J_ITERSUMMARY, J_BOUNDTIGHTENING,"\n");
      }

      return -1; // declare this node infeasible
    }

    /*if ((var -> Type () == VAR) &&
	(var -> isInteger ())) {
      lower_i = ceil  (lower_i - COUENNE_EPS);
      upper_i = floor (upper_i + COUENNE_EPS);
      }*/

    if (var -> Type () != AUX)
      continue;

    // TODO: also test if any indep variable of this expression
    // have changed. If not, skip

    CouNumber ll, uu; 
    //ll = (*(variables_ [i] -> Lb ())) (),
    //uu = (*(variables_ [i] -> Ub ())) ();

    var -> Image () -> getBounds (ll, uu);

    if (var -> isInteger ()) {
      if (var -> sign () != expression::AUX_LEQ) ll = ceil  (ll - COUENNE_EPS);
      if (var -> sign () != expression::AUX_GEQ) uu = floor (uu + COUENNE_EPS);
    }

    if      (var -> sign () == expression::AUX_LEQ) ll = (*(var -> Lb ())) ();
    else if (var -> sign () == expression::AUX_GEQ) uu = (*(var -> Ub ())) ();

    if (ll - uu > COUENNE_BOUND_PREC * (1 + CoinMin (fabs (ll), fabs (uu)))) {

      //if (Jnlst()->ProduceOutput(J_ITERSUMMARY, J_BOUNDTIGHTENING)) {

      Jnlst()->Printf(J_ITERSUMMARY, J_BOUNDTIGHTENING,
		      "w_%d has infeasible bounds [%g,%g]: ", i, ll, uu);

      if (dbgOutput) {
	var -> Lb () -> print (std::cout);
	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING," --- ");
	var -> Ub () -> print (std::cout);
      }

      Jnlst()->Printf(J_ITERSUMMARY, J_BOUNDTIGHTENING,"\n");
      //}

      return -1; // declare this node infeasible
    }

    // Enter the sign of this auxiliary: if defined as w <= f(x),
    // then the lower bound on f(x) shouldn't be propagated to
    // w. similarly, if defined as w >= f(x), then the same should
    // hold for the upper bound.

    if (var -> sign () == expression::AUX_LEQ) ll = -COUENNE_INFINITY;
    if (var -> sign () == expression::AUX_GEQ) uu =  COUENNE_INFINITY;

    // check if lower bound has tightened
    if ((ll > - COUENNE_INFINITY) && 
	(ll >= lower_i + COUENNE_EPS) &&
	((fabs (ll)        < COUENNE_EPS) || 
	 (fabs (lower_i) < COUENNE_EPS) ||
	 (fabs (ll / (lower_i) - 1) > COUENNE_EPS)) ) {

      if (dbgOutput) {

	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
			"  prop L %2d [%g,(%g)] -> [%g,(%g)] (%g) ", 
			i, lower_i, upper_i, ll, uu, lower_i - ll);
	var -> print (std::cout);

	if (Jnlst()->ProduceOutput(J_MOREDETAILED, J_BOUNDTIGHTENING)) {
	  Jnlst()->Printf(J_MOREDETAILED, J_BOUNDTIGHTENING," := ");
	  var -> Image () -> print (std::cout);
	}

	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"\n");

	if (knownOptimum && 
	    (knownOptimum [i] >= lower_i) && 
	    (knownOptimum [i] <= ll - COUENNE_EPS)) {

	  Jnlst()->Printf(J_STRONGWARNING, J_BOUNDTIGHTENING,
			  "Couenne: propagating l_%d cuts optimum: [%g --> %g -X-> %g] :: ", 
			  i, lower_i, knownOptimum [i], ll);
	  var -> Lb () -> print (std::cout);
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING," --- ");
	  var -> Ub () -> print (std::cout);
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"\n");
	}
      }

      lower_i = ll;

      if (ll > upper_i + COUENNE_BOUND_PREC * (1. + CoinMin (fabs (ll), fabs (upper_i)))) {
	Jnlst () -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING,
			    "just-check: w_%d has infeasible bounds [%g,%g]. ", i, lower_i, upper_i);
	return -1;
      }

      chg_bds [i].setLower (t_chg_bounds::CHANGED);
      nchg++;
    }

    // check if upper bound has tightened
    if ((uu < COUENNE_INFINITY) && 
	(uu <= upper_i - COUENNE_EPS) &&
	((fabs (uu)      < COUENNE_EPS) || 
	 (fabs (upper_i) < COUENNE_EPS) ||
	 (fabs (uu / (upper_i) - 1) > COUENNE_EPS)) ) {
      //      if ((uu < COUENNE_INFINITY) && (uu <= ub_ [i+j] - COUENNE_EPS)) {

      /*printf ("update ubound %d: %.10e <= %.10e - %.12e (%.12e)\n", 
	i+j, uu, ub_ [i+j], COUENNE_EPS, uu - ub_ [i+j]);*/
      /*printf ("update ubound %d: %g >= %g\n", 
	i+j, uu, ub_ [i+j]);*/

      if (dbgOutput) {

	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
			"  prop U %2d [(%g),%g] -> [(%g),%g] (%g) ", 
			i, lower_i, upper_i, ll, uu, upper_i - uu);
	var -> print (std::cout);

	if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," := ");
	  var -> Image () -> print (std::cout);
	}

	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"\n");

	if (knownOptimum && 
	    (knownOptimum [i] <= upper_i) && 
	    (knownOptimum [i] >= uu + COUENNE_EPS)) {

	  Jnlst()->Printf(J_STRONGWARNING, J_BOUNDTIGHTENING,
			  "Couenne: propagating u_%d cuts optimum: [%g <-X- %g <-- %g] :: ", 
			  i, uu, knownOptimum [i], upper_i);
	  var -> Lb () -> print (std::cout);
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING," --- ");
	  var -> Ub () -> print (std::cout);
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"\n");
	}
      }

      upper_i = uu;

      if (uu < lower_i - COUENNE_BOUND_PREC) {
	Jnlst () -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING,
			    "just-check: w_%d has infeasible bounds [%g,%g]. ", i, lower_i, upper_i);
	return -1;
      }

      chg_bds [i].setUpper(t_chg_bounds::CHANGED);
      nchg++;
    }

    // useless if assume expression::[lu]b_ etc already point to
    // problem::[lu]b_
  }

  if (nchg)
    Jnlst () -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING,
			"  forward tightening %d changes\n", nchg);

  return nchg;
}
