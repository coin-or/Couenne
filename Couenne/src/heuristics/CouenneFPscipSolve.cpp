/* $Id$
 *
 * Name:    CouenneFPscipSolve.cpp
 * Authors: Timo Berthold, ZIB Berlin
 * Purpose: call SCIP
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneExprVar.hpp"

#include "cons_rowcuts.h"

using namespace Couenne;

#ifdef COIN_HAS_SCIP

/* general SCIP includes */
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_bounddisjunction.h"

SCIP_RETCODE CouenneFeasPump::ScipSolve (const double *nSol, double* &sol, int niter, int* nsuciter, CouNumber &obj) {

  static int currentmilpmethod = 0;

  int depth = (model_ -> currentNode ()) ? model_ -> currentNode () -> depth () : 0;

  SCIP* scip;

  SCIP_VAR** vars;
  SCIP_Real* lbs;
  SCIP_Real* ubs;
  const SCIP_Real* objs;
  const char* vartypes;
  const CoinPackedMatrix * matrix;
  const CoinBigIndex* rowstarts;
  const int* rowlengths;
  const SCIP_Real* coeffs;
  const SCIP_Real* lhss;
  const SCIP_Real* rhss;
  const int* indices; 

  SCIP_Real timelimit;

  SCIP_VAR      **tabuvars;
  SCIP_Real      *tabubounds;
  SCIP_BOUNDTYPE *tabuboundtypes;

  double infinity;
  int nvars;
  int nconss;
  int nscipsols;

  bool solveagain;

  // COUENNE_INFINITY , getInfinity()

  // get problem data
  nvars    = milp_ -> getNumCols ();
  nconss   = milp_ -> getNumRows ();
  infinity = milp_ -> getInfinity ();

  // get variable data
  lbs =      CoinCopyOfArray (milp_ -> getColLower (), milp_ -> getNumCols ());
  ubs =      CoinCopyOfArray (milp_ -> getColUpper (), milp_ -> getNumCols ());

  objs =     milp_ -> getObjCoefficients ();
  vartypes = milp_ -> getColType ();

  // tighten bounds with data from problem_ (which may have run btCore())

  for (int i = 0; i < problem_->nVars(); ++i)
    if (problem_ -> Var (i) -> Multiplicity () > 0) {
      if (problem_ -> Lb (i) > lbs [i]) lbs [i] = problem_ -> Lb (i);
      if (problem_ -> Ub (i) < ubs [i]) ubs [i] = problem_ -> Ub (i);
    }

  // get row data
  lhss = milp_ -> getRowLower ();
  rhss = milp_ -> getRowUpper ();

  // get matrix data
  matrix     = milp_ -> getMatrixByRow();
  rowstarts  = matrix -> getVectorStarts();
  rowlengths = matrix -> getVectorLengths();
  coeffs     = matrix -> getElements();
  indices    = matrix -> getIndices();
     
  // access tabuPool_ in this class

  if (problem_ -> Jnlst () -> ProduceOutput (Ipopt::J_ERROR, J_NLPHEURISTIC)) {
    SCIPdebugMessage("create SCIP problem instance with %d variables and %d constraints.\n", nvars, nconss);
  }

  // initialize SCIP
  SCIP_CALL( SCIPcreate (&scip) );
  assert(scip != NULL);
   
  // include default SCIP plugins
  SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

  // include row cut constraint hanlder
  if( milpCuttingPlane_ == FP_CUT_INTEGRATED )
    { 
      SCIP_CALL( SCIPincludeConshdlrRowcuts(scip, couenneCG_, milp_) );
    }

  // create problem instance in SCIP
  SCIP_CALL( SCIPcreateProb(scip, "auxiliary FeasPump MILP", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

  // allocate local memory for SCIP variable array
  SCIP_CALL( SCIPallocMemoryArray(scip, &vars, nvars) );

  // Allocating space for tabu stuff

  SCIP_CALL( SCIPallocMemoryArray(scip, &tabuvars      , 2*nvars) ); // fix to nvars + nIntvars
  SCIP_CALL( SCIPallocMemoryArray(scip, &tabubounds    , 2*nvars) );
  SCIP_CALL( SCIPallocMemoryArray(scip, &tabuboundtypes, 2*nvars) );
     
  // one variable for objective !!!!!!!!!

  // create variables 
  for (int i=0; i<nvars; i++) {

    char varname[SCIP_MAXSTRLEN];  

    bool neglect = 
      (i < problem_ -> nVars ()) && 
      (problem_ -> Var (i) -> Multiplicity () <= 0);

    // check that all data is in valid ranges
    assert( 0 <= vartypes[i] && vartypes[i] <= 2);

    if (!neglect) {

      checkInfinity(scip, lbs[i], infinity);
      checkInfinity(scip, ubs[i], infinity);
    }

    // all variables are named x_i
    (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", i);
    SCIP_CALL( SCIPcreateVar(scip, &vars[i], varname, 
			     CoinMin (lbs [i], ubs [i]),
			     CoinMax (lbs [i], ubs [i]),
			     objs [i], 
			     vartypes[i] == 0 ? SCIP_VARTYPE_CONTINUOUS : (vartypes[i] == 1 ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER),
			     TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

    // add the variable to SCIP
    SCIP_CALL( SCIPaddVar(scip, vars[i]) );
  }

  for (std::set <CouenneFPsolution, compareSol>:: iterator i = tabuPool_ . begin (); 
       i != tabuPool_ . end (); ++i) {

    const double *x = i -> x ();

    int nEntries = 0;

    SCIP_CONS *tabucons = NULL;

    for (int j = 0; j < i -> n (); ++j) {

      if (problem_ -> Var (j) -> isInteger () && 
	  (problem_ -> Var (j) -> Multiplicity () > 0) &&
	  (fabs (ubs [j] - lbs [j]) > .5)) {

	// problem_ -> Ub  (j) - problem_ -> Lb (j) > .5) {

	// if (fabs (x [j] - floor (x [j] + .5)) >= SCIPfeastol (scip)) {
	//   printf ("integer var x%d not really integer: %e\n", j, x [j]);
	// }

	//assert (fabs (lbs [j] - problem_ -> Lb (j)) < SCIPfeastol (scip));
	//assert (fabs (ubs [j] - problem_ -> Ub (j)) < SCIPfeastol (scip));
	assert (fabs (x [j] - floor (x [j] + .5))   < SCIPfeastol (scip) * 1.e3);

	assert (nEntries <= 2*nvars - 2);

	double x_rounded = floor (x [j] + .5);

	if        (x [j] >= ubs [j] - COUENNE_EPS) {

	  tabuvars       [nEntries]   = vars [j];
	  tabubounds     [nEntries]   = x_rounded - 1.;
	  tabuboundtypes [nEntries++] = SCIP_BOUNDTYPE_UPPER;

	} else if (x [j] <= lbs [j] + COUENNE_EPS) {

	  tabuvars       [nEntries]   = vars [j];
	  tabubounds     [nEntries]   = x_rounded + 1.;
	  tabuboundtypes [nEntries++] = SCIP_BOUNDTYPE_LOWER;

	} else {

	  tabuvars       [nEntries]   = vars [j];
	  tabubounds     [nEntries]   = x_rounded - 1.;
	  tabuboundtypes [nEntries++] = SCIP_BOUNDTYPE_UPPER;

	  tabuvars       [nEntries]   = vars [j];
	  tabubounds     [nEntries]   = x_rounded + 1.;
	  tabuboundtypes [nEntries++] = SCIP_BOUNDTYPE_LOWER;
	}
      } 
    }

    if (nEntries != 0) {

      SCIP_CALL (SCIPcreateConsBounddisjunction (scip, &tabucons, "Tabu Solution", nEntries,
						 tabuvars, tabuboundtypes, tabubounds,  
						 TRUE,  TRUE,  TRUE,  TRUE,  TRUE, 
						 FALSE, FALSE, FALSE, FALSE, FALSE));
      
      SCIP_CALL( SCIPaddCons(scip, tabucons) );
      //SCIP_CALL (SCIPprintCons (scip, tabucons, NULL));
      SCIP_CALL (SCIPreleaseCons (scip, &tabucons));
    }
  }

  // create constraints
  for (int i=0; i<nconss; i++) {

    SCIP_CONS* cons;

    char consname[SCIP_MAXSTRLEN];
    (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "row_%d", i);

    // check that all data is in valid ranges
    checkInfinity(scip, lhss[i], infinity);
    checkInfinity(scip, rhss[i], infinity);

    // create an empty linear constraint
    SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhss[i], rhss[i], 
				    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    // add variables to constraint
    for(int j=rowstarts[i]; j<rowstarts[i]+rowlengths[i]; j++)        
      {
	checkInfinity(scip, coeffs[j], infinity);
	SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[indices[j]], coeffs[j]) );
      }

    // add constraint to SCIP
    SCIP_CALL( SCIPaddCons(scip, cons) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons) );        
  }
  
  // determine the method to solve the MILP
  if (milpMethod_ == 0 && niter == 0)
    {
      // initialize currentmilpmethod: at the root node we solve the MILP with SCIP default
      // deeper in the tree, we will solve the MILP by a FP heuristic
      if( depth == 0 )
	currentmilpmethod = 2;
      else 
	currentmilpmethod = 4;
    }
  else if (milpMethod_ != 0)
    currentmilpmethod = milpMethod_; // use a fixed method to solve the MILP

  // MILP solving loop. If the MILP terminates without a solution, it might get resolved with a more expensive atrategy
  do {
    solveagain = false;

    // reset parameters if MILP is solved agian
    SCIP_CALL( SCIPresetParams(scip) );
        
    // deactivate SCIP output
    if (!(problem_ -> Jnlst () -> ProduceOutput (Ipopt::J_WARNING, J_NLPHEURISTIC))) {
      SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
    }
        
    // do not abort subproblem on CTRL-C
    SCIP_CALL( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );

    // set time limit
    timelimit = problem_ -> getMaxCpuTime () - CoinCpuTime ();

    if (timelimit < 0) 
      break;

    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );        

    if (problem_ -> Jnlst () -> ProduceOutput (Ipopt::J_WARNING, J_NLPHEURISTIC)) {
      SCIPinfoMessage(scip, NULL, "using MILP method: %d\n",currentmilpmethod);
    }
        
    // tune SCIP differently, depending on the chosen method to solve the MILP
    /// -1. MILP and stop at first solution found (similar to other FP implementations and mainly done for comparison with them)
    ///  0. Let Couenne choose
    ///  1. Partially solve the MILP with emphasis on good solutions
    ///  2. Solve the MILP relaxation partially, up to a certain node limit
    ///  3. Apply RENS to 1
    ///  4. Use Objective FP 2.0 for MILPs
    switch(currentmilpmethod)
      {
      case -1: // solve the MILP completely. SCIP's default setting should be best for this
	if( milpCuttingPlane_ == FP_CUT_INTEGRATED )
	  {
	    SCIP_CALL( SCIPsetLongintParam(scip, "constraints/rowcuts/maxcuttingrounds", 0) );
	  }

	// As this is expensive, stop as soon as a solution is found

	SCIP_CALL( SCIPsetIntParam(scip, "limits/bestsol", 1) );
	break;

      case 1: // Be aggressive in finding feasible solutions, but lazy about the dual bound. 
	// Set limits on overall nodes and stall nodes (nodes without incumbent improvement).
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/stallnodes", 1000) );
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 10000) );
	SCIP_CALL( SCIPsetRealParam   (scip, "limits/gap", .001) );

	// disable expensive cutting plane separation 
	//SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_FAST, TRUE) );
        
	// disable expensive presolving 
	//SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

	// use aggressive primal heuristics 
	SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
        
	// use best estimate node selection
	/*
	  if( SCIPfindNodesel(scip, "estimate") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/estimate/stdpriority", INT_MAX/4) ); 
	  }
	*/

	// use inference branching 
	/*
	  if( SCIPfindBranchrule(scip, "inference") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
	  }
	*/

	// disable conflict analysis 
	/*
	  SCIP_CALL( SCIPsetBoolParam(scip, "conflict/useprop", FALSE) );
	  SCIP_CALL( SCIPsetBoolParam(scip, "conflict/useinflp", FALSE) );
	  SCIP_CALL( SCIPsetBoolParam(scip, "conflict/useboundlp", FALSE) );
	  SCIP_CALL( SCIPsetBoolParam(scip, "conflict/usesb", FALSE) );
	  SCIP_CALL( SCIPsetBoolParam(scip, "conflict/usepseudo", FALSE) );
	*/

	break;

      case 2: // default SCIP with node limits
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/stallnodes", 500) );
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 5000) );
	SCIP_CALL( SCIPsetRealParam   (scip, "limits/gap", .005) );

	// disable expensive dual techniques

	SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_FAST, TRUE) );
	SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );
	SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", INT_MAX/4) );

	break;
      case 3: // solve the MILP with RENS. Disable most other features, enable RENS
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );

	// disable cutting plane separation 
	SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
        
	// disable expensive presolving 
	SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

	// besides RENS, only use cheap heuristics 
	SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );

	// use inference branching 
	if( SCIPfindBranchrule(scip, "inference") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
	  }

	// ensure that RENS is called
	if( SCIPfindHeur(scip, "rens") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", 0) );
	    SCIP_CALL( SCIPsetRealParam(scip, "heuristics/rens/minfixingrate", 0.0) );
	  }
	if( milpCuttingPlane_ == FP_CUT_INTEGRATED )
	  { 
	    SCIP_CALL( SCIPsetLongintParam(scip, "constraints/rowcuts/maxcuttingrounds", 0) );
	  }
	break;

      case 4: // solve the MILP with Feasibility Pump. Disable most other features, enable stage 3 for feaspump
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );

	// disable cutting plane separation 
	SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
        
	// disable expensive presolving 
	SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

	// besides feaspump, only use cheap heuristics 
	SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );

	// use inference branching 
	if( SCIPfindBranchrule(scip, "inference") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
	  }

	// ensure that feasibility pump is called
	if( SCIPfindHeur(scip, "feaspump") != NULL )
	  {
	    SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", 0) );
	    SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/maxsols", -1) );

	    // Comment next 8 lines as Workaround for mac/conv/batch
	    // to avoid assertion in stage 3 on nsolutions==0

	    if( SCIPversion() <= 2.01 )
	      {
	    	SCIP_CALL( SCIPsetBoolParam(scip, "heuristics/feaspump2/stage3", TRUE) );
	      }
	    else
	      {
	    	SCIP_CALL( SCIPsetBoolParam(scip, "heuristics/feaspump/stage3", TRUE) );
	      }
	  }
	if( milpCuttingPlane_ == FP_CUT_INTEGRATED )
	  { 
	    SCIP_CALL( SCIPsetLongintParam(scip, "constraints/rowcuts/maxcuttingrounds", 0) );
	  }
	break;

      case 5: // use rounding heuristic in Couenne

	SCIP_CALL( SCIPsetLongintParam(scip, "limits/stallnodes", 500) );
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 5000) );
	SCIP_CALL( SCIPsetRealParam   (scip, "limits/gap", .05) );

	{
	  if (!sol)
	    sol = new CouNumber [problem_ -> nVars ()];

	  if (!nSol)
	    nSol = milp_ -> getColSolution ();

	  // whatever nSol, pick solution through getIntCand

	  problem_ -> getIntegerCandidate (nSol, sol, problem_ -> Lb (), problem_ -> Ub ());

	  SCIP_SOL* scipSol;
	  SCIP_CALL( SCIPcreateSol(scip, &scipSol, NULL ));
	  SCIP_CALL(SCIPsetSolVals(scip, scipSol, nvars, vars,(double*) sol));
	}

	break;

      case 6: // round; TODO: but perturb first if we are cycling

	SCIP_CALL( SCIPsetLongintParam(scip, "limits/stallnodes", 500) );
	SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 5000) );
	SCIP_CALL( SCIPsetRealParam   (scip, "limits/gap", .05) );

	if (!sol)
	  sol = new CouNumber [problem_ -> nVars ()];

	if (!nSol)
	  nSol = milp_ -> getColSolution ();

	for (int j = 0; j < problem_ -> nVars (); ++j)
	  sol [j] =
	    (problem_ -> Var (j) -> Multiplicity () <= 0) ? 0. :
	    (problem_ -> Var (j) -> isInteger    ())      ? COUENNE_round (nSol [j]) : nSol [j];

	  SCIP_SOL* scipSol;
	  SCIP_CALL( SCIPcreateSol(scip, &scipSol, NULL ));
	  SCIP_CALL(SCIPsetSolVals(scip, scipSol, nvars, vars, (double*) sol));

	break;


      default:
	printf("Invalid MILP method in feasibility pump: %d\n", currentmilpmethod);
	assert(false);
	break;
      }

    ////////////////////////////////////////////////////////////////////////

    SCIP_CALL(SCIPtransformProb(scip));

    {
      SCIP_CONSHDLR*  conshdlr_bounddisjunction = SCIPfindConshdlr(scip, "bounddisjunction");
      int nbdconss= SCIPconshdlrGetNConss(conshdlr_bounddisjunction);
      SCIP_CONS** bdconss = SCIPconshdlrGetConss(conshdlr_bounddisjunction);

      for(int i=0; i < nbdconss; i++)
	{
          SCIP_CONS* checkcons = bdconss[i];
          SCIP_VAR** checkvars = SCIPgetVarsBounddisjunction(scip, checkcons);
          SCIP_Real* checkbounds = SCIPgetBoundsBounddisjunction(scip,checkcons);
          SCIP_BOUNDTYPE* checkboundtypes =SCIPgetBoundtypesBounddisjunction(scip, checkcons);
          int checknvars = SCIPgetNVarsBounddisjunction(scip, checkcons);

          for(int j=i+1; j < nbdconss; j++)
	    {
	      SCIP_CONS* tmpcons =SCIPconshdlrGetConss(conshdlr_bounddisjunction)[j];
	      SCIP_VAR      **tmpvars = SCIPgetVarsBounddisjunction(scip,tmpcons);
	      SCIP_Real      *tmpbounds =SCIPgetBoundsBounddisjunction(scip, tmpcons);
	      SCIP_BOUNDTYPE *tmpboundtypes =SCIPgetBoundtypesBounddisjunction(scip, tmpcons);
	      int tmpnvars = SCIPgetNVarsBounddisjunction(scip, tmpcons);
	      int k;

	      if( checknvars !=  tmpnvars )
                continue;

	      for(k=0; k < tmpnvars; k++)
                if(tmpvars[k] != checkvars[k] || tmpbounds[k] !=checkbounds[k] || tmpboundtypes[k] != checkboundtypes[k])
		  break;

	      if (k == tmpnvars)
		problem_ -> Jnlst () -> Printf (Ipopt::J_WARNING, J_NLPHEURISTIC, "ZZZ identical bounddisjunction constraints\n");
	    }
	}
    }

    //////////////////////////////////////////////////////////////////////////

#if 0
    // writes MILP problem and SCIP settings into a file
    SCIP_CALL( SCIPwriteOrigProblem(scip, "debug.lp", NULL, FALSE) );
    SCIP_CALL( SCIPwriteParams(scip, "debug.set", FALSE,TRUE) );
#endif
     
    // solve the MILP
                                                //   /|
                                                //  / +---------+
    SCIP_RETCODE retcode = SCIPsolve (scip);    // <     MILP   |
                                                //  \ +---------+
                                                //   \|
#if 0
    if (SCIPgetStatus (scip) == SCIP_STATUS_INFEASIBLE) {

      // writes MILP problem and SCIP settings into a file
      SCIP_CALL( SCIPwriteOrigProblem(scip, "debug.lp", NULL, FALSE) );
      SCIP_CALL( SCIPwriteParams(scip, "debug.set", FALSE, TRUE));

      printf ("SCIP found that the problem is infeasible, exiting\n");

      exit (-1);
    }
#endif

    if (problem_ -> Jnlst () -> ProduceOutput (Ipopt::J_WARNING, J_NLPHEURISTIC))
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );

    // Greppable line with condensed info

    problem_ -> Jnlst () -> Printf (Ipopt::J_ERROR, J_NLPHEURISTIC, "[FeasPump-SCIP] %5d %5d %7.2f\n", SCIPgetNVars(scip), SCIPgetNConss(scip), SCIPgetSolvingTime(scip));

    if (retcode != SCIP_OKAY) {
      problem_ -> Jnlst () -> Printf (Ipopt::J_WARNING, J_NLPHEURISTIC, "Couenne FP: SCIPsolve did not succeed\n");
      goto TERMINATION;
    }

    nscipsols =  SCIPgetNSols(scip);
     
    // copy the solution
    if( nscipsols)
      {
	SCIP_SOL** scipsols;
	SCIP_SOL* bestsol;
	SCIP_Real cutoffbound;

	int nstoredsols;

	double objval = 0;

	/* get incumbent solution */
	bestsol = SCIPgetBestSol(scip);
	assert(bestsol != NULL);

	/* get SCIP solution pool */
	scipsols = SCIPgetSols(scip);
	assert(scipsols != NULL);

	if (!sol)
	  sol = new CouNumber [problem_ -> nVars ()];

	// get solution values and objective of incumbent
	SCIP_CALL( SCIPgetSolVals(scip, bestsol, problem_ -> nVars (), vars, sol) );
	obj = SCIPgetSolOrigObj(scip, bestsol);

	nstoredsols = 0;

	// if we use an objective feaspump, the objective function might be negative
	cutoffbound = obj > 0 ? 2*obj : obj / 2.0;

	// insert other SCIP solutions into solution pool
	// do not store too many or too poor solutions 
	for(int i=1; i<nscipsols && nstoredsols < 10 && 
	      (objval = SCIPgetSolOrigObj(scip,scipsols[i])) <= cutoffbound; i++){

	  double* tmpsol = new CouNumber [nvars];
           
	  // get solution values
	  SCIP_CALL( SCIPgetSolVals(scip, scipsols[i], problem_ -> nVars (), vars, tmpsol) );
	  CouenneFPsolution couennesol = CouenneFPsolution (problem_, tmpsol);

	  delete [] tmpsol;

	  // add solutions to the pool if they are not in the tabu list
	  if (   tabuPool_      . find (couennesol) == tabuPool_      . end () 
	      && pool_ -> Set (). find (couennesol) == pool_ -> Set() . end ()
	      ){

	    pool_ -> Set (). insert (couennesol);

	    ++nstoredsols;
	  }
	}

	++(*nsuciter);

	// if we succeeded three times in a row, try a cheaper MILP_ solving method next time
	// TODO: if we want to use time limits, hitting the time limit would be another good reason to switch
	if( milpMethod_ == 0 && *nsuciter >= 3 && currentmilpmethod < 4 )
	  {
	    ++currentmilpmethod;
	    *nsuciter = 0;
	  }          
      }
    //try to use a more aggressive, more expensive way to solve the MILP
    else if( milpMethod_ == 0 && currentmilpmethod > 1 )
      {
	--currentmilpmethod;
	solveagain = true;
	*nsuciter = 0;

	// throw away the current solution process
	SCIP_CALL( SCIPfreeTransform(scip) );
      }
    else {

      obj = COIN_DBL_MAX;  
    }

  } while (solveagain);

  ////////////////////////////////////////////////////////////////

 TERMINATION:
   
  // release variables before freeing them
  for (int i=0; i<nvars; i++) {
    SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
  }

  // free memory
  SCIPfreeMemoryArray(scip, &vars);

  // free tabu stuff
  SCIPfreeMemoryArray(scip, &tabuvars      );
  SCIPfreeMemoryArray(scip, &tabubounds    );
  SCIPfreeMemoryArray(scip, &tabuboundtypes);

  delete [] lbs;
  delete [] ubs;

  SCIP_CALL( SCIPfree(&scip) );
   
  BMScheckEmptyMemory();     

  return SCIP_OKAY;
}

#endif
