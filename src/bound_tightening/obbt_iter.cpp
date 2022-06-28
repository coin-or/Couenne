/*
 *
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"

using namespace Ipopt;
using namespace Couenne;

#define OBBT_EPS 1e-3

// TODO: seems like Clp doesn't like large bounds and crashes on
// explicit bounds around 1e200 or so. For now simply use fictitious
// bounds around 1e14. Fix.

int obbt_supplement (const OsiSolverInterface *csi, /// interface to use as a solver
		     int index,                     /// variable being looked at
		     int sense);                    /// 1: minimize, -1: maximize

/// reoptimize and change bound of a variable if needed
static bool obbt_updateBound (OsiSolverInterface *csi, /// interface to use as a solver
			      int sense,               /// 1: minimize, -1: maximize
			      CouNumber &bound,        /// bound to be updated
			      bool isint) {            /// is this variable integer


  // NEW: TODO: save min x^\star_i to check at every iteration
  // (comparison each iteration is less advtgs) when checking if x_i
  // needs minz for obbt

  //csi -> deleteScaleFactors ();
  csi -> setDblParam (OsiDualObjectiveLimit, COIN_DBL_MAX);
  csi -> setDblParam (OsiPrimalObjectiveLimit, (sense==1) ? bound : -bound);
  //csi -> setObjSense (1); // always minimize, just change the sign of the variable // done in caller

  ////////////////////////////////////////////////////////////////////////

  //csi -> resolve_nobt (); // this is a time-expensive part, be considerate...
  csi -> resolve (); // this is a time-expensive part, be considerate...

  ////////////////////////////////////////////////////////////////////////

  if (csi -> isProvenOptimal ()) {

    double opt = csi -> getObjValue ();

    if (sense > 0)
         {if (opt        > bound + OBBT_EPS) {bound = (isint ? ceil  (opt - COUENNE_EPS) : opt); return true;}}
    else {if ((opt=-opt) < bound - OBBT_EPS) {bound = (isint ? floor (opt + COUENNE_EPS) : opt); return true;}}
  }

  return false;
}


/// Iteration on one variable

int CouenneProblem::obbt_iter (OsiSolverInterface *csi,
			       t_chg_bounds *chg_bds,
			       const CoinWarmStart *warmstart,
			       Bonmin::BabInfo *babInfo,
			       double *objcoe,
			       int sense,
			       int index) const {

  // TODO: do NOT apply OBBT if this is a variable of the form
  // w2=c*w1, as it suffices to multiply result. More in general, do
  // not apply if w2 is a unary monotone function of w1. Even more in
  // general, if w2 is a unary function of w1, apply bound propagation
  // from w1 to w2 and mark it as exact (depending on whether it is
  // non-decreasing or non-increasing

  //static int iter = 0;

  // exclude checking known optimal solution if initial bounding box
  // already excludes it

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

  std::set <int> deplist;
  int deplistsize;

  bool issimple = false;

  exprVar *var = Var (index);

  int
    objind  = Obj (0) -> Body () -> Index (),
    ncols   = csi -> getNumCols (),
    nImprov = 0;

  if ((var -> Type  () == AUX) &&
      ((deplistsize = var -> Image () -> DepList (deplist, STOP_AT_AUX)) <= 1)) {

    if (!deplistsize) { // funny, the expression is constant...

      CouNumber value = (*(var -> Image ())) ();

      if (csi -> getColLower () [index] < value - COUENNE_EPS) {
	csi -> setColLower (index, value);
	chg_bds    [index].setLowerBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
      }
      else chg_bds [index].setLowerBits(t_chg_bounds::EXACT);

      if (csi -> getColUpper () [index] > value + COUENNE_EPS) {
	csi -> setColUpper (index, value);
	chg_bds    [index].setUpperBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
      }
      else chg_bds [index].setUpperBits(t_chg_bounds::EXACT);

      issimple = true;

    } else { // the expression only depends on one variable, meaning
	     // that bound propagation should be sufficient

      int indInd = *(deplist.begin ());

      //      expression *image = var -> Image ();
      // TODO: write code for monotone functions...

      if // independent variable is exactly bounded in both ways
	((((chg_bds [indInd].lower() & t_chg_bounds::EXACT) &&
	   (chg_bds [indInd].upper() & t_chg_bounds::EXACT)) ||
	  // or if this expression is of the form w=cx+d, that is, it
	  // depends on one variable only and it is linear
	  (var -> Image () -> Linearity () <= LINEAR)) &&
	 (var -> sign () == expression::AUX_EQ)) {

	issimple = true;

	CouNumber lb, ub;
	var -> Image () -> getBounds (lb, ub);

	if (lb < Lb (index)) lb = Lb (index);
	if (ub > Ub (index)) ub = Ub (index);

	if (var -> isInteger ()) {
	  lb = ceil  (lb - COUENNE_EPS);
	  ub = floor (ub + COUENNE_EPS);
	}

	if (lb > csi -> getColLower () [index] + COUENNE_EPS) {
	  csi -> setColLower (index, lb);
	  Lb (index) = lb;
	  chg_bds      [index].setLowerBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
	} else chg_bds [index].setLowerBits(t_chg_bounds::EXACT);

	if (ub < csi -> getColUpper () [index] - COUENNE_EPS) {
	  csi -> setColUpper (index, ub);
	  Ub (index) = ub;
	  chg_bds      [index].setUpperBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
	} else chg_bds [index].setUpperBits(t_chg_bounds::EXACT);
      }
    }
  }

  // only improve bounds if
  if (!issimple &&
      ((Var (index) -> Type () == VAR) ||        // it is an original variable
       (Var (index) -> Multiplicity () > 0)) &&  // or its multiplicity is at least 1
      (Lb (index) < Ub (index) - COUENNE_EPS) && // in any case, bounds are not equal

      ((index != objind) // this is not the objective
       // or it is, so we use it for re-solving // TODO: check!
       || ((sense ==  1) && !(chg_bds [index].lower() & t_chg_bounds::EXACT))
       )) {
       //((sense==-1) && (psense == MAXIMIZE) && !(chg_bds [index].upper() & t_chg_bounds::EXACT)))) {

    bool isInt = (Var (index) -> isInteger ());

    objcoe [index] = sense;

    csi -> setObjective (objcoe);
    csi -> setObjSense (1); // minimization

    // TODO: Use something else!
#if 0
    for (int iv=0; iv<csi->getNumCols (); iv++) {
      if (fabs (csi -> getColLower () [iv]) > 1e7) csi -> setColLower (iv, -1e14);
      if (fabs (csi -> getColUpper () [iv]) > 1e7) csi -> setColUpper (iv,  1e14);
    }
#endif

    CouNumber &bound =
      (sense == 1) ?
      (Lb (index)) :
      (Ub (index));

    // m{in,ax}imize xi on csi

    /*if (Jnlst()->ProduceOutput(J_MOREVECTOR, J_BOUNDTIGHTENING)) {
      Jnlst()->Printf(J_MOREVECTOR, J_BOUNDTIGHTENING,
		      "m%simizing x%d [%g,%g] %c= %g\n",
	    (sense==1) ? "in" : "ax", index, Lb (index), Ub (index),
	    (sense==1) ? '>'  : '<',  bound); fflush (stdout);
      if (Jnlst()->ProduceOutput(J_MOREMATRIX, J_BOUNDTIGHTENING)) {
	char fname [20];
	sprintf (fname, "m%s_w%03d_%03d", (sense == 1) ? "in" : "ax", index, iter);
	printf ("saving in %s\n", fname);
	//Jnlst()->Printf(J_MOREVECTOR, J_BOUNDTIGHTENING,"writing %s\n", fname);
	csi -> writeLp (fname);
      }
      }*/

    csi -> setWarmStart (warmstart);
    //csi -> continuousModel () -> setPerturbation (50);

    /* From ClpSimplex.cpp:

       If you are re-using the same matrix again and again then the
       setup time to do scaling may be significant.  Also you may not
       want to initialize all values or return all values (especially
       if infeasible).  While an auxiliary model exists it will be
       faster.  If options -1 then model is switched off.  Otherwise
       switched on with following options.

       1 - rhs is constant
       2 - bounds are constant
       4 - objective is constant
       8 - solution in by basis and no djs etc in
       16 - no duals out (but reduced costs)
       32 - no output if infeasible
    */

    //csi -> continuousModel () -> auxiliaryModel (1|8|16|32);

    //Jnlst () -> Printf (J_MATRIX, J_BOUNDTIGHTENING,
    //"obbt___ index = %d [sen=%d,bd=%g,int=%d]\n",
    //index, sense, bound, isInt);

    bool has_updated = false;

    if (obbt_updateBound (csi, sense, bound, isInt)) {

      has_updated = true;

      if (knownOptimum) {
	if (sense == 1) {
	  if (bound       > COUENNE_EPS + knownOptimum [index])
	    Jnlst()->Printf(J_STRONGWARNING, J_BOUNDTIGHTENING,
			    "#### OBBT cuts optimum at x%d: lb = %g, opt = %g, new lb = %g\n",
			    index, Lb (index), knownOptimum [index], bound);
	} else {
	  if (bound       < -COUENNE_EPS + knownOptimum [index])
	    Jnlst()->Printf(J_STRONGWARNING, J_BOUNDTIGHTENING,
			    "#### OBBT cuts optimum at x%d: ub = %g, opt = %g, new ub = %g\n",
			    index, Ub (index), knownOptimum [index], bound);
	}
      }

      if (sense == 1) {if (bound > Lb (index)) Lb (index) = bound;}
      else            {if (bound < Ub (index)) Ub (index) = bound;}

      // more conservative, only change (and set CHANGED) if improve

      if (sense==1)
	if (csi -> getColLower () [index] < bound - COUENNE_EPS) {
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"l_%d: %g --> %g\n",
			  index, csi -> getColLower () [index], bound);
	  csi -> setColLower (index, bound);
	  chg_bds      [index].setLowerBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
	} else chg_bds [index].setLowerBits(t_chg_bounds::EXACT);
      else
	if (csi -> getColUpper () [index] > bound + COUENNE_EPS) {
	  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,"u_%d: %g --> %g\n",
			  index, csi -> getColUpper () [index], bound);
	  csi -> setColUpper (index, bound);
	  chg_bds      [index].setUpperBits(t_chg_bounds::CHANGED | t_chg_bounds::EXACT);
	} else chg_bds [index].setUpperBits(t_chg_bounds::EXACT);

      /*
      if (sense==1) {csi -> setColLower (index, bound); chg_bds [index].lower |= CHANGED | EXACT;}
      else          {csi -> setColUpper (index, bound); chg_bds [index].upper |= CHANGED | EXACT;}
      */

#if 0
      // re-check considering reduced costs (more expensive)

      CouNumber *redcost = NULL;

      // first, compute reduced cost when c = c - e_i, where e_i is
      // a vector with all zero except a one in position i. This
      // serves as a base to compute modified reduced cost below.

      for (int j=0; j<ncols; j++)
	if ((j!=index) && (j!=objind)) {

	  // fake a change in the objective function and compute
	  // reduced cost. If resulting vector is all positive
	  // (negative), this solution is also optimal for the
	  // minimization (maximization) of x_j and the corresponding
	  // chg_bds[j].lower (.upper) can be set to EXACT.

	  if (!(chg_bds [j].lower & EXACT)) {
	  }

	  if (!(chg_bds [j].upper & EXACT)) {
	  }
	}
#endif	

      // re-apply bound tightening -- here WITHOUT reduced cost
      // (first argument =NULL is pointer to solverInterface) as csi
      // is not our problem

      Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
		      "  OBBT: x_%d: [%g, %g]\n", index,
		      csi -> getColLower () [index],
		      csi -> getColUpper () [index]);

      // should be faster

      //if (doFBBT_ && !(boundTightening (chg_bds, babInfo))) {
      if (doFBBT_ && !(btCore (chg_bds))) {
	Jnlst () -> Printf (J_DETAILED, J_BOUNDTIGHTENING,
			"node is infeasible after post-OBBT tightening\n");
	return -1; // tell caller this is infeasible
      }

      nImprov++;
    }

    // Check value and bounds of other variables. Do this regardless
    // of the i-th variable being tightened in this iteration

    const double *sol = csi -> getColSolution ();

    for (int j=0; j<ncols; j++)
      if ((j!=index) && (j!=objind)) {

	if (sol [j] <= Lb (j) + COUENNE_EPS) {

	  // if (!(chg_bds [j].lower() & t_chg_bounds::EXACT) && (!has_updated))
	  //   printf ("told you it would be useful: l_i: %g\n", j, Lb (j));
	
	  chg_bds [j].setLowerBits(t_chg_bounds::EXACT);
	}

	if (sol [j] >= Ub (j) - COUENNE_EPS) {

	  // if (!(chg_bds [j].upper() & t_chg_bounds::EXACT) && (!has_updated))
	  //   printf ("told you it would be useful: u_i: %g\n", j, Ub (j));
	
	  chg_bds [j].setUpperBits(t_chg_bounds::EXACT);
	}
      }

    // check for bounds using dual info

    int result = obbt_supplement (csi, index, sense);

    // if we solved the problem on the objective function's
    // auxiliary variable (that is, we re-solved the extended
    // problem), it is worth updating the current point (it will be
    // used later to generate new cuts).

    // TODO: is it, really? And shouldn't we check the opt sense too?
    /*
    if ((objind == index) && (csi -> isProvenOptimal ()) && (sense == 1))
      update (csi -> getColSolution (), NULL, NULL);
    */
    // restore null obj fun
    objcoe [index] = 0;
  }

  if (nImprov && jnlst_ -> ProduceOutput (J_ITERSUMMARY, J_BOUNDTIGHTENING)) {
    Jnlst () -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING, "OBBT: tightened ", nImprov);
    Var (index) -> print ();
    Jnlst () -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING, "\n");
  }

  return nImprov;
}
