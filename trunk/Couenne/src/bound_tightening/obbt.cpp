/* $Id$
 *
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CglCutGenerator.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CouenneExpression.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneInfeasCut.hpp"

using namespace Ipopt;
using namespace Couenne;

#define THRESH_OBBT_AUX 50 // if more than this originals, don't do OBBT on auxs
#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100
#define MAX_OBBT_ATTEMPTS 1 // number of OBBT iterations at root node
			    // -- fixed at one as for some instance it
			    // doesn't seem to do anything after first run

// minimum #bound changed in obbt to generate further cuts
#define THRES_NBD_CHANGED 1

// maximum number of obbt iterations
#define MAX_OBBT_ITER 1

// defined in generateCuts.cpp
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);


// OBBT for one sense (max/min) and one class of variables (orig/aux)
int CouenneProblem::call_iter (OsiSolverInterface *csi, 
			       t_chg_bounds *chg_bds, 
			       const CoinWarmStart *warmstart, 
			       Bonmin::BabInfo *babInfo,
			       double *objcoe,
			       enum nodeType type,
			       int sense) const {

  int ncols   = csi -> getNumCols (),
      nimprov = 0;

  for (int ii=0; ii<ncols; ii++) {

    if (CoinCpuTime () > maxCpuTime_)
      break;

    int i = evalOrder (ii);

    enum expression::auxSign aSign = Var (i) -> sign ();

    if ((Var (i) -> Type () == type)     &&
	(Var (i) -> Multiplicity () > 0) &&
	((type == VAR)                               || 
	 (aSign  == expression::AUX_EQ) ||
	 ((aSign == expression::AUX_LEQ) && (sense > 0)) ||
	 ((aSign == expression::AUX_GEQ) && (sense < 0)))) {

      int ni = obbt_iter (csi, chg_bds, warmstart, babInfo, objcoe, sense, i);

//       {
// 	// ToDo: Pipe all output through journalist
// 	Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING, 
// 			"  bounds after obbt step  =====================\n  ");
// 	int j=0;
// 	for (int i=0; i < nVars (); i++) 
// 	  if (variables_ [i] -> Multiplicity () > 0) {
// 	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
// 			    "x_%03d [%+10g %+10g] ", i, 
// 			    domain_. lb (i),
// 			    domain_. ub (i));
// 	    if (!(++j % 6)) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n  ");
// 	  }
// 	if (j % 6) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
//       }

      if (ni < 0) return ni;
      nimprov += ni;
    }
  }

  return nimprov;
}


/// Optimality based bound tightening -- inner loop

int CouenneProblem::obbtInner (OsiSolverInterface *csi,
			       OsiCuts &cs,
			       t_chg_bounds *chg_bds,
			       Bonmin::BabInfo * babInfo) const {

  // set large bounds to infinity (as per suggestion by JJF)

  int ncols = csi -> getNumCols ();
  const double *lb = csi -> getColLower (),
               *ub = csi -> getColUpper ();

  double inf = csi -> getInfinity ();

  for (int i=ncols; i--;) {
    if (lb [i] < - COUENNE_INFINITY) csi -> setColLower (i, -inf);
    if (ub [i] >   COUENNE_INFINITY) csi -> setColUpper (i,  inf);
  }

  //  csi -> setHintParam (OsiDoDualInResolve, false);

  // setup cloned interface for later use
  csi -> setObjSense (1); // minimization
  csi -> setIntParam (OsiMaxNumIteration, MAX_OBBT_LP_ITERATION);
  csi -> applyCuts (cs);   // apply all (row+column) cuts to formulation
  csi -> initialSolve ();

  const CoinWarmStart *warmstart = csi -> getWarmStart ();

  // improve each bound

  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  csi -> setObjective (objcoe);
  csi -> setObjSense (1);        // minimization

  int nimprov = 0;
 
  const int Infeasible = 1;

  try {

    int ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, VAR,  1)) < 0) throw Infeasible;
    nimprov += ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, VAR, -1)) < 0) throw Infeasible;
    nimprov += ni;

    if (nVars () < THRESH_OBBT_AUX) {

      if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, AUX,  1)) < 0) throw Infeasible;
      nimprov += ni;

      if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, AUX, -1)) < 0) throw Infeasible;
      nimprov += ni;
    }
  }

  catch (int exception) {

    if (exception == Infeasible)
      nimprov = -1;
  }

  free (objcoe);
  delete warmstart;

  return (nimprov);
}


// Optimality based bound tightening -- main loop
int CouenneProblem::obbt (const CouenneCutGenerator *cg,
			  const OsiSolverInterface &si,
			  OsiCuts &cs,
			  const CglTreeInfo &info,
			  Bonmin::BabInfo * babInfo,
			  t_chg_bounds *chg_bds) {

  // TODO: set up list of hopeless variables and do different OBBT
  // saving lps

  // Check if cs contains only one cut and if it is of the form 1 <=
  // x0 <= -1. That means a previous cut generator has determined that
  // this node is infeasible and we shouldn't take the pain of running
  // this CGL.

  if (isWiped (cs) || info.pass >= MAX_OBBT_ATTEMPTS)
    return 0;

  int nTotImproved = 0;

  // Do OBBT if:
  if (doOBBT_ &&                        // flag is checked, AND
      ((logObbtLev_ != 0) ||               // (parameter is nonzero OR
       (info.level == 0)) &&               //  we are at root node), AND
      (info.pass == 0) &&               // at first round of cuts, AND 
      ((logObbtLev_ < 0) ||               // (logObbtLev = -1, OR
       (info.level <= logObbtLev_) ||     //  depth is lower than COU_OBBT_CUTOFF_LEVEL, OR
                                          //  probability inversely proportional to the level)
       (CoinDrand48 () < pow (2., (double) logObbtLev_ - (info.level + 1))))) {

    if ((info.level <= 0 && !(info.inTree)) || 
    	jnlst_ -> ProduceOutput (J_STRONGWARNING, J_COUENNE))  {

      jnlst_ -> Printf (J_ERROR, J_COUENNE, "Optimality Based BT: "); 
      //nVars () > THRESH_OBBT_AUX ? nOrigVars_ : nVars (), info.pass); 
      fflush (stdout);
    }

    jnlst_ -> Printf (J_ITERSUMMARY, J_BOUNDTIGHTENING, "----- OBBT\n");

    // TODO: why check info.pass==0? Why not more than one pass? It
    // should be anyway checked that info.level be >= 0 as <0 means
    // first call at root node

    OsiSolverInterface *csi = si.clone (true);

    csi -> messageHandler () -> setLogLevel (0);
    //dynamic_cast <CouenneSolverInterface<T> *> 

    OsiClpSolverInterface *clpcsi = dynamic_cast <OsiClpSolverInterface *> (csi);

    if (clpcsi)
      clpcsi -> setupForRepeatedUse ();
    //csi -> doingResolve () = false;

    //csi -> setHintParam (OsiDoDualInResolve, false);

    int nImprov, nIter = 0;

    bool notImproved = false;

    while (!notImproved && 
	   (nIter++ < MAX_OBBT_ITER) &&
	   ((nImprov = obbtInner (csi, cs, chg_bds, babInfo)) > 0) &&
	   (CoinCpuTime () < maxCpuTime_)) {

      int nchanged, *changed = NULL;

      /// OBBT has tightened, add improved bounds
      sparse2dense (nVars (), chg_bds, changed, nchanged);
      cg -> genColCuts (*csi, cs, nchanged, changed);

      nTotImproved += nImprov;

      if ((nIter < MAX_OBBT_ITER) && 
	  (nImprov >= THRES_NBD_CHANGED)) {

	// only generate new row cuts if improvents are enough
	int nCurCuts = cs.sizeRowCuts ();
	cg -> genRowCuts (*csi, cs, nchanged, changed, chg_bds);

	if (nCurCuts == cs.sizeRowCuts ())
	  notImproved = true; // repeat only if new cuts available

      } else notImproved = true;

      if (changed) 
	free (changed);
    }

    //csi -> doingResolve () = true;

    delete csi;

    if ((info.level <= 0 && !(info.inTree)) ||
    	jnlst_ -> ProduceOutput (J_STRONGWARNING, J_COUENNE))
      jnlst_ -> Printf (J_ERROR, J_COUENNE, "%d improved bounds\n", nTotImproved);

    if (nImprov < 0) {
      jnlst_->Printf(J_ITERSUMMARY, J_BOUNDTIGHTENING, "  Couenne: infeasible node after OBBT\n");
      return -1;
    }
  }

  return 0;
}
