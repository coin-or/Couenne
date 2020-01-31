/* $Id$
 *
 * Name:    genRowCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate Row Cuts for current convexification
 *
 * (C) Carnegie-Mellon University, 2006-09. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

using namespace Couenne;

/// generate OsiRowCuts for current convexification
void CouenneCutGenerator::genRowCuts (const OsiSolverInterface &si,
				      OsiCuts &cs,
				      int nchanged, 
				      int *changed,
				      t_chg_bounds *chg) const {

  // TODO: make nchanged and changed useful
  // TODO: pass have_NLP to all

  // For each auxiliary variable, create convexification cut (or set
  // of cuts) and add it to cs

  // do NOT project current point into new bounding box! it prevents
  // convexification cuts that are violated by current point

  /*for (int i=0, j = problem_ -> nVars (); j--; i++) {

    CouNumber &x = problem_ -> X (i),
      lb = problem_ -> Lb (i),
      ub = problem_ -> Ub (i);

    if (x < lb) x = lb;
    else if (x > ub) x = ub;
    }*/

  if (firstcall_)
    for (int i=0, j = problem_ -> nVars (); j--; i++) {

      if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	break;

      exprVar *var = problem_ -> Var (i);

      if ((var -> Multiplicity () > 0) && 
	  (var -> Type () == AUX)) {

	var -> generateCuts (cs, this, chg);
      }
    }
  else { // chg_bds contains the indices of the variables whose bounds
	 // have changed (a -1 follows the last element)

    for (int i = 0, j = problem_ -> nVars (); j--; i++) {

      // TODO: check if list contains all and only aux's to cut

      /*expression * image = problem_ -> Aux (i) -> Image ();
      
      if ((image -> dependsOn (changed, nchanged)) && 
	  (image -> Linearity () > LINEAR)) {
	printf ("         ");
	problem_ -> Aux (i) -> print ();
	printf (" : = ");
	image -> print ();
	printf ("\n");
      }
      */
      // cut only if:

      /*if (   (image -> Linearity () > LINEAR)    // 1) expression is non linear
	&& (image -> dependsOn (changed, nchanged) // 2) it depends on changed variables
	|| have_NLP
	|| info.pass > 0)) {
      */

      exprVar *var = problem_ -> Var (problem_ -> evalOrder (i));

      if ((var -> Type () == AUX) &&
	  (var -> Multiplicity () > 0) &&
	  (var -> Image () -> Linearity () > LINEAR)) {

	if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	  break;

	var -> generateCuts (cs, this, chg);
      }
    }
  }
}
