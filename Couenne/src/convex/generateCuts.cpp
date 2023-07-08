/* $Id$
 *
 * Name:    generateCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: the generateCuts() method of the convexification class CouenneCutGenerator
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "BonCbc.hpp"
#include "BonBabInfos.hpp"
#include "CglCutGenerator.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneInfeasCut.hpp"

#include "CouenneRecordBestSol.hpp"

//#define FM_PRINT_INFO

#ifdef COIN_HAS_NTY
#include "Nauty.h"
#endif

using namespace Ipopt;

namespace Couenne {

#define Couenne_large_bound2 9.99e12

// checks bad cuts against known optimum
bool isOptimumCut (const CouNumber *opt, OsiCuts &cs, CouenneProblem *p);

// set and lift bound for auxiliary variable associated with objective
// function
void fictitiousBound (OsiCuts &cs,
		      CouenneProblem *p, 
		      bool action) {     // true before convexifying, false afterwards

  // fictitious bound for initial unbounded lp relaxations
  const CouNumber large_tol = (Couenne_large_bound2 / 1e6);

  // set trivial dual bound to objective function, if there is none

  int ind_obj = p -> Obj (0) -> Body () -> Index ();

  if (ind_obj < 0) return;

  // we have a single variable objective function

  //int sense = -1; //(p -> Obj (0) -> Sense () == MINIMIZE) ? -1 : 1;

  if (action)
    //if (sense<0) 
      {if (p -> Lb (ind_obj) < - Couenne_large_bound2) p -> Lb (ind_obj) = - Couenne_large_bound2;}
  //else         {if (p -> Ub (ind_obj) >   large_bound2) p -> Ub (ind_obj) =   large_bound2;}
  else
    //if (sense>0) {if (fabs (p->Ub(ind_obj)-large_bound2)<large_tol) p->Ub(ind_obj)=COUENNE_INFINITY;}
    //else         
      {if (fabs (p->Lb(ind_obj)+Couenne_large_bound2)<large_tol) p->Lb(ind_obj) =-COUENNE_INFINITY;}
}


// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged) {

  // convert sparse chg_bds in something handier

  changed  = (int *) realloc (changed, ncols * sizeof (int));
  nchanged = 0;

  for (int i=ncols, j=0; i--; j++, chg_bds++)
    if (chg_bds -> lower() != t_chg_bounds::UNCHANGED ||
	chg_bds -> upper() != t_chg_bounds::UNCHANGED ) {
      *changed++ = j;
      nchanged++;
    }

  changed -= nchanged;
  //changed = (int *) realloc (changed, nchanged * sizeof (int));
}


/// get new bounds from parents' bounds + branching rules
void updateBranchInfo (const OsiSolverInterface &si, CouenneProblem *p, 
		       t_chg_bounds *chg, const CglTreeInfo &info);

/// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					const CglTreeInfo info)
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
  const
#endif
  {

  // if si.lb(objInd) > cutoff,
  //   return infeasCut

  int indObj = problem_ -> Obj (0) -> Body () -> Index ();
  
  if ((indObj >= 0) && 
      (si.getColLower () [indObj] > problem_ -> getCutOff () + COUENNE_EPS)) {

    WipeMakeInfeas (cs);
    return;
  }

  // check if out of time or if an infeasibility cut (iis of type 0)
  // was added as a result of, e.g., pruning on BT. If so, no need to
  // run this.

  if (isWiped (cs) || 
     (CoinCpuTime () > problem_ -> getMaxCpuTime ()))
    return;

#ifdef FM_TRACE_OPTSOL
  double currCutOff = problem_->getCutOff();
  double bestVal = 1e50;
  CouenneRecordBestSol *rs = problem_->getRecordBestSol();
  if(rs->getHasSol()) {
    bestVal = rs->getVal(); 
  }
  if(currCutOff > bestVal) {
    //problem_ -> setCutOff (bestVal - 1e-6); // FIXME: don't add numerical constants
    problem_ -> setCutOff (bestVal);

    if ((indObj >= 0) && (si. getColUpper () [indObj] > bestVal)) {
      OsiColCut *objCut = new OsiColCut;
      objCut->setUbs(1, &indObj, &bestVal);
      cs.insert(objCut);
      delete objCut;
    }
  }
#endif

#ifdef FM_PRINT_INFO
  if((BabPtr_ != NULL) && (info.level >= 0) && (info.pass == 0) && 
     (BabPtr_->model().getNodeCount() > lastPrintLine)) {
    printLineInfo();
    lastPrintLine += 1;
  }
#endif

  const int infeasible = 1;

  int nInitCuts = cs.sizeRowCuts ();

  CouNumber
    *&realOpt = problem_ -> bestSol (),
    *saveOptimum = realOpt;

  if (!firstcall_ && realOpt) { 

    // have a debug optimal solution. Check if current bounds
    // contain it, otherwise pretend it does not exist

    CouNumber *opt = realOpt;

    const CouNumber 
      *sol = si.getColSolution (),
      *lb  = si.getColLower (),
      *ub  = si.getColUpper ();

    int objind = problem_ -> Obj (0) -> Body () -> Index ();

    for (int j=0, i=problem_ -> nVars (); i--; j++, opt++, lb++, ub++)
      if ((j != objind) && 
	  ((*opt < *lb - COUENNE_EPS * (1 + CoinMin (fabs (*opt), fabs (*lb)))) || 
	   (*opt > *ub + COUENNE_EPS * (1 + CoinMin (fabs (*opt), fabs (*ub)))))) {
	
	jnlst_ -> Printf (J_VECTOR, J_CONVEXIFYING, 
			  "out of bounds, ignore x%d = %g [%g,%g] opt = %g\n", 
			  problem_ -> nVars () - i - 1, *sol, *lb, *ub, *opt);

	// optimal point is not in current bounding box,
	// pretend realOpt is NULL until we return from this procedure
	realOpt = NULL;
	break;
      }
  }

  /*static int count = 0;
  char fname [20];
  sprintf (fname, "relax_%d", count++);
  si.writeLp (fname);
  printf ("writing %s\n", fname);*/

  jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,
		    "generateCuts: level = %d, pass = %d, intree = %d\n",
		    info.level, info.pass, info.inTree);

  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (si.getAuxiliaryInfo ());

  if (babInfo)
    babInfo -> setFeasibleNode ();

  double now   = CoinCpuTime ();
  int    ncols = problem_ -> nVars ();

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  t_chg_bounds *chg_bds = new t_chg_bounds [ncols];

  /*for (int i=0; i < ncols; i++) 
    if (problem_ -> Var (i) -> Multiplicity () <= 0) {
      chg_bds [i].setLower (t_chg_bounds::UNCHANGED);
      chg_bds [i].setUpper (t_chg_bounds::UNCHANGED);
      }*/

  problem_ -> installCutOff (); // install upper bound

  if (firstcall_) {

    // First convexification //////////////////////////////////////

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to update the auxiliary variables and bounds with
    // their current value.

    for (int i=0; i < ncols; i++) 
      if (problem_ -> Var (i) -> Multiplicity () > 0) {
	chg_bds [i].setLower (t_chg_bounds::CHANGED);
	chg_bds [i].setUpper (t_chg_bounds::CHANGED);
      }

    // start with FBBT, should take advantage of cutoff found by NLP
    // run AFTER initial FBBT...
    if (problem_ -> doFBBT () &&
	(! (problem_ -> boundTightening (chg_bds, info, babInfo))))
          jnlst_ -> Printf (J_STRONGWARNING, J_CONVEXIFYING,
            "Couenne: WARNING, first convexification is infeasible\n");

    // For each auxiliary variable replacing the original (nonlinear)
    // constraints, check if corresponding bounds are violated, and
    // add cut to cs

    int nnlc = problem_ -> nCons ();

    for (int i=0; i<nnlc; i++) {

      if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	break;

      // for each constraint
      CouenneConstraint *con = problem_ -> Con (i);

      // (which has an aux as its body)
      int objindex = con -> Body () -> Index ();

      if ((objindex >= 0) && 
	  ((con -> Body () -> Type () == AUX) ||
	   (con -> Body () -> Type () == VAR))) {

	// get the auxiliary that is at the lhs
	exprVar *conaux = problem_ -> Var (objindex);

	if (conaux &&
	    (conaux -> Type () == AUX) &&
	    (conaux -> Image ()) && 
	    (conaux -> Image () -> Linearity () <= LINEAR)) {

	  // reduce density of problem by adding w >= l rather than
	  // ax + b >= l for any linear auxiliary defined as w := ax+b

	  double 
	    lb = (*(con -> Lb ())) (), 
	    ub = (*(con -> Ub ())) ();

	  OsiColCut newBound;
	  if (lb > -COUENNE_INFINITY) newBound.setLbs (1, &objindex, &lb);
	  if (ub <  COUENNE_INFINITY) newBound.setUbs (1, &objindex, &ub);

	  cs.insert (newBound);

	  // the auxiliary w of constraint w <= b is associated with a
	  // linear expression w = ax: add constraint ax <= b
	  /*conaux -> Image () -> generateCuts (conaux, si, cs, this, chg_bds, 
					      conaux -> Index (), 
					      (*(con -> Lb ())) (), 
					      (*(con -> Ub ())) ());*/

	  // take it from the list of the variables to be linearized
	  // 
	  // DO NOT decrease multiplicity. Even if it is a linear
	  // term, its bounds can still be used in implied bounds
	  //
	  // Are we sure? That will happen only if its multiplicity is
	  // nonzero, for otherwise this aux is only used here, and is
	  // useless elsewhere
	  //
	  //conaux -> decreaseMult (); // !!!
	}

	// also, add constraint w <= b

	// not now, do it later

// 	// if there exists violation, add constraint
// 	CouNumber l = con -> Lb () -> Value (),	
// 	          u = con -> Ub () -> Value ();

// 	// tighten bounds in Couenne's problem representation
// 	problem_ -> Lb (index) = CoinMax (l, problem_ -> Lb (index));
// 	problem_ -> Ub (index) = CoinMin (u, problem_ -> Ub (index));

      } else { // body is more than just a variable, but it should be
	       // linear. If so, generate equivalent linear cut

	assert (false);	// TODO
      }
    }

    if (jnlst_ -> ProduceOutput (J_ITERSUMMARY, J_CONVEXIFYING)) {
      if (cs.sizeRowCuts ()) {
	jnlst_ -> Printf (J_ITERSUMMARY, J_CONVEXIFYING,"Couenne: %d constraint row cuts\n",
			  cs.sizeRowCuts ());
	for (int i=0; i<cs.sizeRowCuts (); i++) 
	  cs.rowCutPtr (i) -> print ();
      }
      if (cs.sizeColCuts ()) {
	jnlst_ -> Printf (J_ITERSUMMARY, J_CONVEXIFYING,"Couenne: %d constraint col cuts\n",
			  cs.sizeColCuts ());
	for (int i=0; i<cs.sizeColCuts (); i++) 
	  cs.colCutPtr (i) -> print ();
      }
    }
  } else {

    // use new optimum as lower bound for variable associated w/objective

    // transmit solution from OsiSolverInterface to problem
    problem_ -> domain () -> push (&si, &cs);

    if (indObj >= 0) {

      // Use current value of objvalue's x as a lower bound for bound
      // tightening
      double lp_bound = problem_ -> domain () -> x (indObj);

      //if (problem_ -> Obj (0) -> Sense () == MINIMIZE) 
      {if (lp_bound > problem_ -> Lb (indObj)) problem_ -> Lb (indObj) = lp_bound;}
	   //else {if (lp_bound < problem_ -> Ub (indObj)) problem_ -> Ub (indObj) = lp_bound;}
    }

    updateBranchInfo (si, problem_, chg_bds, info); // info.depth >= 0 || info.pass >= 0
  }

  // restore constraint bounds before tightening and cut generation
  for (int i = problem_ -> nCons (); i--;) {

    // for each constraint
    CouenneConstraint *con = problem_ -> Con (i);

    // (which has an aux as its body)
    int objindex = con -> Body () -> Index ();

    if ((objindex >= 0) && 
	((con -> Body () -> Type () == AUX) ||
	 (con -> Body () -> Type () == VAR))) {

      // if there exists violation, add constraint
      CouNumber 
	l = con -> Lb () -> Value (),	
	u = con -> Ub () -> Value ();

      // tighten bounds in Couenne's problem representation
      problem_ -> Lb (objindex) = CoinMax (l, problem_ -> Lb (objindex));
      problem_ -> Ub (objindex) = CoinMin (u, problem_ -> Ub (objindex));
    }
  }

  problem_ -> installCutOff (); // install upper bound

  fictitiousBound (cs, problem_, false); // install finite lower bound, if currently -inf

  int *changed = NULL, nchanged;

  // Bound tightening ///////////////////////////////////////////

  // do bound tightening only at first pass of cutting plane in a node
  // of BB tree (info.pass == 0) or if first call (creation of RLT,
  // info.pass == -1)

  try {

    // Before bound tightening, compute symmetry group. After bound
    // tightening is done, we can apply further tightening using orbit
    // information.

// #ifdef COIN_HAS_NTY
//     //    ChangeBounds (psi -> getColLower (),  
//     //		  psi -> getColUpper (), 
//     //		  psi -> getNumCols ());

//     if (problem_ -> orbitalBranching ()){

//       problem_ -> ChangeBounds (problem_ -> Lb (),
// 				problem_ -> Ub (),
// 				problem_ -> nVars ());

//       problem_ -> Compute_Symmetry ();
//     }
// #endif

    // Bound tightening ////////////////////////////////////

    //bool is_feas = p -> btCore (chg_bds);

    // Reduced Cost BT -- to be done first to use rcost correctly
    if (!firstcall_  &&                         // have a linearization already
	problem_ -> doRCBT () &&                // authorized to do reduced cost tightening
	problem_ -> redCostBT (&si, chg_bds) && // some variables were tightened with reduced cost
	!(problem_ -> btCore (chg_bds)))        // in this case, do another round of FBBT
      throw infeasible;

    // FBBT
    if (problem_ -> doFBBT () && 
	//(info.pass <= 0) && // do it in subsequent rounds too
	(! (problem_ -> boundTightening (chg_bds, info, babInfo))))
      throw infeasible;

    // OBBT
    if (!firstcall_ && // no obbt if first call (there is no LP to work with)
	problem_ -> obbt (this, si, cs, info, babInfo, chg_bds) < 0)
      throw infeasible;

    // Bound tightening done /////////////////////////////

    if ((problem_ -> doFBBT () ||
	 problem_ -> doOBBT () ||
	 problem_ -> doABT  ()) &&
	(jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING))) {

      jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== after bt =============\n");
      for (int i = 0; i < problem_ -> nVars (); i++)
	if (problem_ -> Var (i) -> Multiplicity () > 0)
	  jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			 problem_ -> X  (i), problem_ -> Lb (i), problem_ -> Ub (i));
      jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
    }

    // Use orbit info to tighten bounds

#ifdef COIN_HAS_NTY

    // TODO: when independent bound tightener, can get original bounds
    // through si.getCol{Low,Upp}er()

    if (problem_ -> orbitalBranching () && !firstcall_) {

      CouNumber 
	*lb = problem_ -> Lb (),
	*ub = problem_ -> Ub ();

      std::vector<std::vector<int> > *new_orbits = problem_ -> getNtyInfo () -> getOrbits();

      for (int i=0, ii = problem_ -> getNtyInfo () -> getNumOrbits (); ii--; i++){

	CouNumber
	  ll = -COUENNE_INFINITY,
	  uu =  COUENNE_INFINITY; 
	
	std::vector <int> orbit = (*new_orbits)[i];

	if (orbit.size () <= 1)
	  continue; // not much to do when only one variable in this orbit

	if (jnlst_ -> ProduceOutput (J_VECTOR, J_BOUNDTIGHTENING)) {
	  printf ("orbit bounds: "); fflush (stdout);
	  for(int j = 0; j < orbit.size (); j++) {
	    printf ("x_%d [%g,%g] ", orbit[j], lb [orbit [j]], ub [orbit [j]]);
	    fflush (stdout);
	  }
	}

	for (int j = 0; j < orbit.size (); j++) {
 
	  int indOrb = orbit [j];

	  if (indOrb < problem_ -> nVars ()) {

	    if (lb [indOrb] > ll) ll = lb [indOrb];
	    if (ub [indOrb] < uu) uu = ub [indOrb];
	  }
	}

	jnlst_ -> Printf (J_VECTOR, J_BOUNDTIGHTENING, 
			  " --> new common bounds: [%g,%g]\n", ll, uu);

	for(int j = 0; j < orbit.size (); j++) {

	  int indOrb = orbit [j];

	  if (indOrb < problem_ -> nVars ()){

	    lb [indOrb] = ll;
	    ub [indOrb] = uu;
	  }
	}
      }

      delete new_orbits;
    }

#endif

    // Generate convexification cuts //////////////////////////////

    sparse2dense (ncols, chg_bds, changed, nchanged);

    double *nlpSol = NULL;

    //--------------------------------------------

    if (true) {

      if (babInfo) 
	nlpSol = const_cast <double *> (babInfo -> nlpSolution ());

      // Aggressive Bound Tightening ////////////////////////////////

      int logAbtLev = problem_ -> logAbtLev ();

      if (problem_ -> doABT () &&             // flag is checked, AND
	  ((logAbtLev != 0) ||                // (parameter is nonzero OR
	   (info.level == 0)) &&              //  we are at root node), AND
	  (info.pass == 0) &&                 // at first round of cuts, AND 
	  ((logAbtLev < 0) ||                 // (logAbtLev = -1, OR
	   (info.level <= logAbtLev) ||       //  depth is lower than COU_OBBT_CUTOFF_LEVEL, OR
	   (CoinDrand48 () <                  //  probability inversely proportional to the level)
	    pow (2., (double) logAbtLev - (info.level + 1))))) {

	jnlst_ -> Printf(J_VECTOR, J_BOUNDTIGHTENING,"  performing ABT\n");
	if (! (problem_ -> aggressiveBT (nlp_, chg_bds, info, babInfo)))
	  throw infeasible;

	sparse2dense (ncols, chg_bds, changed, nchanged);
      }

      // obtain solution just found by nlp solver

      // Auxiliaries should be correct. solution should be the one found
      // at the node even if not as good as the best known.

      // save violation flag and disregard it while adding cut at NLP
      // point (which are not violated by the current, NLP, solution)
      bool save_av = addviolated_;
      addviolated_ = false;

      // save values
      problem_ -> domain () -> push 
	(problem_ -> nVars (), 
	 problem_ -> domain () -> x  (), 
	 problem_ -> domain () -> lb (), 
	 problem_ -> domain () -> ub (), false);

      // fill originals with nlp values
      if (nlpSol) {
	CoinCopyN (nlpSol, problem_ -> nOrigVars (), problem_ -> domain () -> x ());
      //problem_ -> initAuxs ();

      problem_ -> getAuxs (problem_ -> domain () -> x ());
      }

      if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {
	jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== genrowcuts on NLP =============\n");
	for (int i = 0; i < problem_ -> nVars (); i++)
	  if (problem_ -> Var (i) -> Multiplicity () > 0)
	    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			   problem_ -> X  (i),
			   problem_ -> Lb (i),
			   problem_ -> Ub (i));
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
      }

      problem_ -> domain () -> current () -> isNlp () = true;
      genRowCuts (si, cs, nchanged, changed, chg_bds);  // add cuts

      problem_ -> domain () -> pop (); // restore point

      addviolated_ = save_av;     // restore previous value

      //    if (!firstcall_) // keep solution if called from extractLinearRelaxation()
      if (babInfo) 
	babInfo -> setHasNlpSolution (false); // reset it after use //AW HERE

    } else {

      if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {
	jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== genrowcuts on LP =============\n");
	for (int i = 0; i < problem_ -> nVars (); i++)
	  if (problem_ -> Var (i) -> Multiplicity () > 0)
	    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			   problem_ -> X  (i),
			   problem_ -> Lb (i),
			   problem_ -> Ub (i));
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
      }

      genRowCuts (si, cs, nchanged, changed, chg_bds);
    }

    // change tightened bounds through OsiCuts
    if (nchanged)
      genColCuts (si, cs, nchanged, changed);

    if (firstcall_ && (cs.sizeRowCuts () >= 1))
      jnlst_->Printf(J_ITERSUMMARY, J_CONVEXIFYING,
		     "Couenne: %d initial row cuts\n", cs.sizeRowCuts ());

    if (realOpt && // this is a good time to check if we have cut the optimal solution
	isOptimumCut (realOpt, cs, problem_))
      jnlst_->Printf(J_ITERSUMMARY, J_CONVEXIFYING,
		     "Warning: Optimal solution was cut\n");
  }

  catch (int exception) {

    if ((exception == infeasible) && (!firstcall_)) {

      jnlst_ -> Printf (J_ITERSUMMARY, J_CONVEXIFYING,
			"Couenne: Infeasible node\n");

      WipeMakeInfeas (cs);
    }

    if (babInfo) // set infeasibility to true in order to skip NLP heuristic
      babInfo -> setInfeasibleNode ();
  }

  delete [] chg_bds;

  if (changed) 
    free (changed);

  if (firstcall_) {

    jnlst_ -> Printf (J_SUMMARY, J_CONVEXIFYING, 
		      "Couenne: %d cuts (%d row, %d col) for linearization\n", 
		      cs.sizeRowCuts () + cs.sizeColCuts (),
		      cs.sizeRowCuts (),  cs.sizeColCuts ());

    fictitiousBound (cs, problem_, true);
    firstcall_  = false;
    ntotalcuts_ = nrootcuts_ = cs.sizeRowCuts ();

  } else { 

    problem_ -> domain () -> pop ();

    ntotalcuts_ += (cs.sizeRowCuts () - nInitCuts);

    if (saveOptimum)
      realOpt = saveOptimum; // restore debug optimum
  }

  septime_ += CoinCpuTime () - now;

  if (jnlst_ -> ProduceOutput (J_ITERSUMMARY, J_CONVEXIFYING)) {

    if (cs.sizeColCuts ()) {
      jnlst_ -> Printf (J_ITERSUMMARY, J_CONVEXIFYING,"Couenne col cuts:\n");
      for (int i=0; i<cs.sizeColCuts (); i++) 
	cs.colCutPtr (i) -> print ();
    }
  }

  if (!(info.inTree)) 
    rootTime_ = CoinCpuTime ();
}

}
