/* $Id$
 *
 * Name:    impliedBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: backward implied bound search
 *
 * (C) Carnegie-Mellon University, 2006-10. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

using namespace Couenne;

/// Bound tightening for auxiliary variables

int CouenneProblem::impliedBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

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

  if (Jnlst()->ProduceOutput(Ipopt::J_DETAILED, J_BOUNDTIGHTENING)) {  
    Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING,"  backward =====================\n  ");
    int j=0;
    for (int i=0; i < nVars (); i++) 
      if (variables_ [i] -> Multiplicity () > 0) {
	Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING,
			"x_%03d [%+10g %+10g] ", i, 
			domain_.lb (i),
			domain_.ub (i));
	if (!(++j % 6)) Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING,"\n  ");
      }
    if (j % 6) Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING,"\n");
  }

  for (int ii = nVars (); ii--;) {

    int i = numbering_ [ii];

    if (Lb (i) > Ub (i) &&
	(Lb (i) < Ub (i) + COUENNE_BOUND_PREC * (1 + CoinMin (fabs (Lb (i)), fabs (Ub (i)))))) {

      // This is to prevent very tiny infeasibilities to propagate
      // down and make the problem infeasible. Example pointed out in
      // http://list.coin-or.org/pipermail/couenne/2010-October/000145.html

      CouNumber tmp = Lb (i);
      Lb (i)        = Ub (i);
      Ub (i)        = tmp;
    }

    if ((variables_ [i] -> Type () == AUX) &&
	(variables_ [i] -> Multiplicity () > 0)) {

      if (Lb (i) > Ub (i) + COUENNE_BOUND_PREC * (1 + CoinMin (fabs (Lb (i)), fabs (Ub (i))))) {
	Jnlst () -> Printf (Ipopt::J_DETAILED, J_BOUNDTIGHTENING,
			    "  implied bounds: w_%d has infeasible bounds [%g,%g]\n", 
			    i, Lb (i), Ub (i));
	return -1;
      }

      // TODO: also test if this expression, or any of its indep
      // variables, have changed. If not, skip

      CouNumber 
	l0 = Lb (i), 
	u0 = Ub (i);

      if (variables_ [i] -> Image () -> impliedBound 
	  (variables_ [i] -> Index (), Lb (), Ub (), chg_bds, variables_ [i] -> sign ())) {

	// conservative check for integer variables. 
	/*if (Var (i) -> isInteger ()) {
	  Lb (i) = ceil  (Lb (i) - COUENNE_EPS);
	  Ub (i) = floor (Ub (i) + COUENNE_EPS);
	  }*/

	if (Jnlst()->ProduceOutput(Ipopt::J_VECTOR, J_BOUNDTIGHTENING)) {
	  // todo: send all output through journalist
	  Jnlst()->Printf(Ipopt::J_VECTOR, J_BOUNDTIGHTENING,
			  "  impli %2d [%15.8g, %15.8g] -> [%15.8g, %15.8g]: ",
			  i, l0, u0, Lb (i), Ub (i));

	  variables_ [i] -> print (std::cout);

	  if (Jnlst()->ProduceOutput(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING)) {
	    Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_BOUNDTIGHTENING," := ");
	    variables_ [i] -> Image () -> print (std::cout);
	  }

	  Jnlst()->Printf(Ipopt::J_VECTOR, J_BOUNDTIGHTENING,"\n");
	}

	if (knownOptimum && 
	    ((knownOptimum [i] < Lb (i) - COUENNE_EPS) ||
	     (knownOptimum [i] > Ub (i) + COUENNE_EPS)))

	  Jnlst () -> Printf (Ipopt::J_DETAILED, J_BOUNDTIGHTENING,
			      "#### implied b_%d [%g,%g] cuts optimum %g\n",
			      i, Lb (i), Ub (i), 
			      knownOptimum [i]);

	//printf ("impli %2d ", nvar+i);

	/*if (variables_ [i] -> Image () -> Argument () || 
	    variables_ [i] -> Image () -> ArgList ()) {

	  expression *arg = variables_ [i] -> Image () -> Argument ();

	  if (!arg) {
	    for (int k=0; k < variables_ [i] -> Image () -> nArgs (); k++) {
	      arg =  variables_ [i] -> Image () -> ArgList () [k];
	      Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING," ");
	      arg -> print (std::cout);
	      if (arg -> Index () >= 0) {
		int ind = arg -> Index ();
		Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING,
				" in [%g,%g]", 
				expression::Lbound (ind), 
				expression::Ubound (ind));
	      }	    
	    }
	  } else {
	    Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING," ");
	    arg -> print (std::cout);
	    if (arg -> Index () >= 0) {
	      int ind = arg -> Index ();
	      Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING," in [%g,%g]", 
		      expression::Lbound (ind), 
		      expression::Ubound (ind));
	    }
	  }
	} else Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING," [no args]");
	Jnlst()->Printf(Ipopt::J_DETAILED, J_BOUNDTIGHTENING,"\n");*/

	nchg++;
      }
    }
  }

  if (nchg)
    Jnlst () -> Printf (Ipopt::J_DETAILED, J_BOUNDTIGHTENING, "  implied bounds: %d changes\n", nchg);

  return nchg;
}
