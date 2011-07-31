/* $Id$
 *
 * Name:    fake_tightening.cpp
 * Author:  Pietro Belotti
 * Purpose: fake single bounds in variables to exclude parts of the solution space 
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"

#include "BonBabInfos.hpp"

using namespace Couenne;

#define MAX_ITER  3 // max # fake tightening (inner) iterations 
#define AGGR_MUL  2 // the larger,  the more conservative. Must be > 0
#define AGGR_DIV  2 // the smaller, the more conservative. Must be > 1

// golden ratio, used to find the ideal bound
const CouNumber phi = 0.5 * (1. + sqrt (5.));

// create fictitious bounds to tighten current interval
CouNumber fictBounds (char direction,
		      CouNumber  x,
		      CouNumber  lb,   
		      CouNumber  ub) {

#define LARGE_BOUND 1e10

  if   (lb < -LARGE_BOUND) {
    if (ub >  LARGE_BOUND) { // ]-inf,+inf[

      return (!direction ? -sqrt (-lb) : sqrt (ub));

      //if (fabs (x) < COUENNE_EPS) return (direction ? AGGR_MUL : - AGGR_MUL);
      //else                        return (direction ? AGGR_MUL : - AGGR_MUL) * fabs (x);

    } else { // ]-inf,u]

      if (!direction)
	return -sqrt (-lb); // try to tighten interval from a very large value

      if      (x < -COUENNE_EPS) return (CoinMin (0., (x+ub)/2));
      else if (x >  COUENNE_EPS) return ((x + (ub-x)/AGGR_DIV));
      else                       return ((ub/AGGR_DIV));

      // if      (x < -COUENNE_EPS) return (direction ? CoinMin (0., (x+ub)/2) : AGGR_MUL * x);
      // else if (x >  COUENNE_EPS) return (direction ? (x + (ub-x)/AGGR_DIV) : 0);
      // else                       return (direction ? (ub/AGGR_DIV) : -AGGR_MUL);
    }
  }
  else {
    if (ub >  LARGE_BOUND) { // [l,+inf[

      if (direction)
	return sqrt (ub);

      if      (x < -COUENNE_EPS) return ((x - (x-lb) / AGGR_DIV));
      else if (x >  COUENNE_EPS) return (CoinMax (0.,(x+lb)/2));
      else                       return (lb/AGGR_DIV);

      // if      (x < -COUENNE_EPS) return (direction ? 0 : (x - (x-lb) / AGGR_DIV));
      // else if (x >  COUENNE_EPS) return (direction ? (AGGR_MUL * x) : CoinMax (0.,(x+lb)/2));
      // else                       return (direction ? AGGR_MUL : lb/AGGR_DIV);

    } else // [l,u]
      return (direction ? 
	      (x + (ub-x) / AGGR_DIV) : 
	      (x - (x-lb) / AGGR_DIV));
  }
}


// Single fake tightening. Return
//
// -1   if infeasible
//  0   if no improvement
// +1   if improved
int CouenneProblem::
fake_tighten (char direction,  // 0: left, 1: right
	      int index,       // index of the variable tested
	      const double *X, // point round which tightening is done
	      CouNumber *olb,  // cur. lower bound
	      CouNumber *oub,  //      upper
	      t_chg_bounds *chg_bds,
	      t_chg_bounds *f_chg) const {
  int
    ncols    = nVars (),
    objind   = Obj (0) -> Body  () -> Index ();

  assert (objind >= 0);

  bool 
    tightened = false,
    intvar    = variables_ [index] -> isInteger ();

  CouNumber 
    xcur      = X [index],
    inner     = xcur,                                                 // point closest to current x
    outer     = (direction ? oub : olb) [index],                      // point closest to bound
    fb        = fictBounds (direction, xcur, Lb (index), Ub (index)); // starting point

  // This is a one-dimensional optimization problem between inner and
  // outer, on a monotone function of which we can compute the value
  // (with relative expense) but not the derivative.

  jnlst_ -> Printf (Ipopt::J_ERROR, J_BOUNDTIGHTENING, 
		    "  x_%d.  x = %10g, lb = %g, cutoff = %g-----------------\n", 
		    index,xcur,Lb (objind),getCutOff());

  /*if (index == objind)
    printf ("  x_%d [%g,%g].  x = %10g, break at %g, cutoff = %g-----------------\n", 
    index, Lb (index), Ub (index), xcur, fb, getCutOff());*/

  for (int iter = 0; iter < MAX_ITER; iter++) {

    if (intvar) {

      if (!direction) {inner = floor (inner + COUENNE_EPS); outer = ceil  (outer - COUENNE_EPS);}
      else            {inner = ceil  (inner - COUENNE_EPS); outer = floor (outer + COUENNE_EPS);}

      if ( (direction && (inner > outer + .5)) || // more robust check on integer-valued doubles
	  (!direction && (inner < outer - .5))) {

	// fictitious interval is empty, hence useless to check. 

	// apply new (valid, tightened) bound
	if (direction) {oub[index] = Ub (index) = fb; chg_bds[index].setUpper(t_chg_bounds::CHANGED);}
	else           {olb[index] = Lb (index) = fb; chg_bds[index].setLower(t_chg_bounds::CHANGED);}

	tightened = true;

	if (!(btCore (f_chg))) 
	  return -1;

	CoinCopyN (Lb (), ncols, olb);
	CoinCopyN (Ub (), ncols, oub);

	// restore initial bound. 
	CoinCopyN (chg_bds, ncols, f_chg);
	//CoinCopyN (olb, ncols, Lb ());
	//CoinCopyN (oub, ncols, Ub ());

	break;
      }

      if ( (direction && ((fb < inner) || (fb > outer))) ||
	  (!direction && ((fb > inner) || (fb < outer))))
	fb = 0.5 * (inner + outer);
    }

    if (direction) {
      Lb (index) = intvar ? ceil (fb - COUENNE_EPS)  : fb; 
      f_chg [index].setLower (t_chg_bounds::CHANGED);
    } else {
      Ub (index) = intvar ? floor (fb + COUENNE_EPS) : fb; 
      f_chg [index].setUpper (t_chg_bounds::CHANGED);
    }

    //    (direction ? lb_ : ub_) [index] = fb; 

    if (jnlst_ -> ProduceOutput (Ipopt::J_ERROR, J_BOUNDTIGHTENING)) {
      char c1 = direction ? '-' : '>', c2 = direction ? '<' : '-';
      printf ("    # x%d = %g iter %3d: [%+10g -%c %+10g %c- %+10g] /\\/\\ ",index,xcur,iter,olb[index],c1,fb,c2, oub [index]);
      printf (" [%10g,%10g]<%g,%g>=> ",Lb (index),Ub (index),CoinMin(inner,outer),CoinMax(inner,outer));
    }

    bool
      feasible  = btCore (f_chg),                           // true if feasible with fake bound
      betterbds = Lb (objind) > getCutOff () + COUENNE_EPS; // true if over cutoff

    jnlst_ -> Printf (Ipopt::J_ERROR, J_BOUNDTIGHTENING,
		      " [%10g,%10g] lb = %g {fea=%d,btr=%d} ",
		      Lb (index), Ub (index), Lb (objind),feasible,betterbds);

    if (feasible && !betterbds) {

      // case 1: too tight, move inner out
      inner = fb;

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, Lb ());
      CoinCopyN (oub, ncols, Ub ());

    } else {

      // Here, !feasible || betterbds
      //
      // If !feasible OR
      //    (betterbds and the new lb is above the cutoff)
      //
      // then there is a tightening

      // case 2: tightening valid, apply and move outer in

      //printf (" --> %cbound [x_%d]: %g --> %g",direction?'U':'L',index,(direction?oub:olb)[index],fb);

      if (optimum_ && 
	  ((!direction &&
	    (optimum_ [index] >= olb [index]) && 
	    (optimum_ [index] <= Lb (index) - COUENNE_EPS)) ||
	   (direction &&
	    (optimum_ [index] <= oub [index]) && 
	    (optimum_ [index] >= Ub (index) + COUENNE_EPS)))) {

	jnlst_ -> Printf (Ipopt::J_ERROR, J_BOUNDTIGHTENING, 
			  "fake tightening CUTS optimum: x%d=%g in [%g,%g] but not in [%g,%g]\n",
			  index, olb [index], oub [index], Lb (index), Ub (index));
      }

      /*bool do_not_tighten = false;

      // check if cut best known solution
      if (optimum_) {
	if (direction) {
	  if ((oub [index] > optimum_ [index]) && 
	      (fb          < optimum_ [index])) {
	    printf ("aggressive bt cuts optimum ub %d: %g < %g < %g\n", 
		    index, fb, optimum_ [index], oub [index]);
	    do_not_tighten = true;
	  }
	} else {
	  if ((olb [index] < optimum_ [index]) && 
	      (fb          > optimum_ [index])) {
	    printf ("aggressive bt cuts optimum lb %d: %g < %g < %g\n", 
		    index, olb [index], optimum_ [index], fb);
	    do_not_tighten = true;
	  }
	}
	}*/

      //if (!do_not_tighten) {

      // apply bound
      if (direction) {

	oub [index] = Ub (index) = intvar ? floor (fb + COUENNE_EPS) : fb; 
	chg_bds [index]. setUpper (t_chg_bounds::CHANGED);

      } else {

	olb [index] = Lb (index) = intvar ? ceil  (fb - COUENNE_EPS) : fb; 
	chg_bds [index]. setLower (t_chg_bounds::CHANGED);
      }

      outer = fb; // we have at least a tightened bound, save it 

      tightened = true;
      //}

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, Lb ());
      CoinCopyN (oub, ncols, Ub ());

      //#if BR_TEST_LOG < 0 // for fair testing
      // check tightened problem for feasibility
      if (!(btCore (chg_bds))) {

	jnlst_ -> Printf (Ipopt::J_ERROR, J_BOUNDTIGHTENING, 
			  "\n    pruned by Probing\n");
	return -1;

      } else {

	// bounds further tightened should be saved
	
	CoinCopyN (Lb (), ncols, olb);
	CoinCopyN (Ub (), ncols, oub);
      }
      //#endif
    }

    // TODO: compute golden section
    //fb = (inner + outer) / 2;

    //fb = fictBounds (direction, fb, CoinMin (inner, outer), CoinMax (inner, outer));

    // inner and outer might have to change. Update 

    CouNumber 
      lb = Lb (index),
      ub = Ub (index);

    if ((!direction && ((inner > ub) || (outer < lb))) ||
	( direction && ((inner < lb) || (outer > ub)))) {

      // keep it simple

      inner = direction ? lb : ub;
      outer = direction ? ub : lb;
    }

    CouNumber diff = fabs (inner - outer);

    if (diff <= COUENNE_EPS) break;

    if (diff > 1.) {

      CouNumber L = log (diff) / log (10.);

      if (direction) fb = inner + exp (log (10.) * L/2);
      else           fb = inner - exp (log (10.) * L/2);

    } else fb = (inner + outer)/2;

    //    if () fb = (          inner + (phi-1) * outer) / phi;
    //    else  fb = ((phi-1) * inner +           outer) / phi;

    //	if (!feasible)       
    //    fb = fictBounds (direction, xcur, 
    //	     direction ? lb [index] : outer,
    //	     direction ? outer      : ub [index]);

    jnlst_ -> Printf (Ipopt::J_ERROR, J_BOUNDTIGHTENING, "\n");
  }

  Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING, "\n");
  if (tightened) 
    Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING, 
		    "  [x%2d] pruned %s [%g, %g] -- lb = %g cutoff = %g\n", 
		    index,direction?"right":"left",
		    olb[index],oub[index], Lb (objind), getCutOff ());

  return tightened ? 1 : 0;
}
