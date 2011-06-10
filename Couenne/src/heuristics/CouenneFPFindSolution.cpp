/* $Id$
 *
 * Name:    CouenneFPFindSolution.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Find solution by looping through MILP solvers/heuristics
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

//#define SCIP_DEBUG 

#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CoinTime.hpp"

using namespace Couenne;

/// find a feasible or optimal solution of MILP
double CouenneFeasPump::findSolution (double* &sol) {

  /// as found on the notes, these methods can be used, from the most
  /// expensive and accurate (exact) method to a cheap, inexact one:
  ///
  /// 1. Solve a MILP relaxation with Manhattan distance as objective
  /// 2. Apply RENS to 1
  /// 3. Use Objective FP 2.0 for MILPs
  /// 4. round-and-propagate
  /// 5. choose from pool, see 4
  /// 6. random perturbation

  // What order should we use? I suggest we use priorities, assigned
  // at the beginning but changeable in the event of multiple failures
  // (or successes) of a given method.
  //
  // Rule of thumb: 
  //
  // 1) Assign all methods i a number p[i] (for instance those in the
  //    list above)
  //
  // 2) Call each in the order define by p[i], return a solution if
  //    found, otherwise proceed to next method
  //
  // 3) If K consecutive successes at finding new solution (not
  //    necessarily new best feasible), --p[i]
  //
  // 4) if H consecutive failutes, ++p[i]

   double obj;

  /// solve MILP 

#ifdef COIN_HAS_SCIP
   if (useSCIP_) {

     SCIP* scip;

     SCIP_VAR** vars;
     const SCIP_Real* lbs;
     const SCIP_Real* ubs;
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

     double infinity;
     int nvars;
     int nconss;

     int currentmilpmethod;

     // COUENNE_INFINITY , getInfinity()

     // get problem data
     nvars    = milp_ -> getNumCols ();
     nconss   = milp_ -> getNumRows ();
     infinity = milp_ -> getInfinity ();

     // get variable data
     lbs =      milp_ -> getColLower ();
     ubs =      milp_ -> getColUpper ();
     objs =     milp_ -> getObjCoefficients ();
     vartypes = milp_ -> getColType ();

     // get row data
     lhss = milp_ -> getRowLower ();
     rhss = milp_ -> getRowUpper ();

     // get matrix data
     matrix     = milp_ -> getMatrixByRow();
     rowstarts  = matrix -> getVectorStarts();
     rowlengths = matrix -> getVectorLengths();
     coeffs     = matrix -> getElements();
     indices    = matrix -> getIndices();
     
     SCIPdebugMessage("create SCIP problem instance with %d variables and %d constraints.\n", nvars, nconss);

     // initialize SCIP
     SCIP_CALL_ABORT( SCIPcreate(&scip) );
     assert(scip != NULL);
   
     // include default SCIP plugins
     SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );
     
#ifndef SCIP_DEBUG
     SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
#endif

     // do not abort subproblem on CTRL-C
     SCIP_CALL( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );

     // set time limit
     timelimit = problem_ -> getMaxCpuTime () - CoinCpuTime ();
     SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );

     // create problem instance in SCIP
     SCIP_CALL_ABORT( SCIPcreateProb(scip, "auxiliary FeasPump MILP", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

     // allocate local memory for SCIP variable array
     SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &vars, nvars) );
   
     // one variable for objective !!!!!!!!!

     // create variables 
     for (int i=0; i<nvars; i++) {
        char varname[SCIP_MAXSTRLEN];  

        // check that all data is in valid ranges
        assert( 0 <= vartypes[i] && vartypes[i] <= 2);
        checkInfinity(scip, lbs[i], infinity);
        checkInfinity(scip, ubs[i], infinity);
        
        // all variables are named x_i
        (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", i);
        SCIP_CALL_ABORT( SCIPcreateVar(scip, &vars[i], varname, lbs[i], ubs[i], objs[i], 
              vartypes[i] == 0 ? SCIP_VARTYPE_CONTINUOUS : (vartypes[i] == 1 ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER),
              TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

        // add the variable to SCIP
        SCIP_CALL_ABORT( SCIPaddVar(scip, vars[i]) );
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
        SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhss[i], rhss[i], 
              TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
             
        // add variables to constraint
        for(int j=rowstarts[i]; j<rowstarts[i]+rowlengths[i]; j++)        
        {
           checkInfinity(scip, coeffs[j], infinity);
           SCIP_CALL_ABORT( SCIPaddCoefLinear(scip, cons, vars[indices[j]], coeffs[j]) );
        }

        // add constraint to SCIP
        SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
        SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );        
     }
     
     // determine the method to solve the MILP
     if (milpMethod_ == 0 )
     {
        // TODO: add rule for automatic choice of currentmilpmethod
        currentmilpmethod = 1;
     }
     else 
        currentmilpmethod = milpMethod_; // use a fixed method to solve the MILP
     
     SCIPdebugMessage("using MILP method: %d\n",currentmilpmethod);

     // tune SCIP differently, depending on the chosen method to solve the MILP
     switch(currentmilpmethod)
     {
     case 1: // solve the MILP completely. SCIP's default setting should be best for this
        break;

     case 2: // solve the MILP quickly. Set limits on overall nodes and stall nodes (nodes without incumbent improvement)
             // disable or cut down all procedures which are merely used for improving the dual bound, e.g., cuts
        SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/stallnodes", 100) );
        SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/nodes", 1000) );

        // disable cutting plane separation 
        SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
        
        // disable expensive presolving 
        SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

        // use aggressuve Heuristics 
        SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
        
        // use best estimate node selection 
        if( SCIPfindNodesel(scip, "estimate") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "nodeselection/estimate/stdpriority", INT_MAX/4) ); 
        }
        
        // use inference branching 
        if( SCIPfindBranchrule(scip, "inference") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
        }
        
        // disable conflict analysis 
        SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "conflict/useprop", FALSE) );
        SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "conflict/useinflp", FALSE) );
        SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "conflict/useboundlp", FALSE) );
        SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "conflict/usesb", FALSE) );
        SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "conflict/usepseudo", FALSE) );

        break;

     case 3: // solve the MILP with RENS. Disable most other features, enable RENS
        SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/nodes", 1) );

        // disable cutting plane separation 
        SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
        
        // disable expensive presolving 
        SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

        // besides RENS, only use cheap heuristics 
        SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );

        // use inference branching 
        if( SCIPfindBranchrule(scip, "inference") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
        }

        // ensure that RENS is called
        if( SCIPfindHeur(scip, "rens") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "heuristics/rens/freq", 0) );
           SCIP_CALL_ABORT( SCIPsetRealParam(scip, "heuristics/rens/minfixingrate", 0.0) );
        }
        break;

     case 4: // solve the MILP with Feasibility Pump. Disable most other features, enable stage 3 for feaspump
        SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/nodes", 1) );

        // disable cutting plane separation 
        SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
        
        // disable expensive presolving 
        SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );

        // besides feaspump, only use cheap heuristics 
        SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );

        // use inference branching 
        if( SCIPfindBranchrule(scip, "inference") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/inference/priority", INT_MAX/4) );
        }

        // ensure that feasibility pump is called
        if( SCIPfindHeur(scip, "feaspump") != NULL )
        {
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "heuristics/feaspump/freq", 0) );
           SCIP_CALL_ABORT( SCIPsetIntParam(scip, "heuristics/feaspump/maxsols", -1) );
           (void) SCIPsetBoolParam(scip, "heuristics/feaspump/stage3", TRUE);
           (void) SCIPsetBoolParam(scip, "heuristics/feaspump2/stage3", TRUE);
        }

        break;

     default:
        break;
     }
  /// 1. Solve a MILP relaxation with Manhattan distance as objective
  /// 2. Partially solve the MILP with emphasis on good solutions
  /// 3. Apply RENS to 1
  /// 4. Use Objective FP 2.0 for MILPs
  /// 5. round-and-propagate
  /// 6. choose from pool, see 4
  /// 7. random perturbation

          
     // solve the MILP
     SCIP_CALL_ABORT( SCIPsolve(scip) );

     obj = SCIPinfinity(scip);
     // copy the solution
     if( SCIPgetNSols(scip) )
     {
        SCIP_SOL* bestsol;
        bestsol = SCIPgetBestSol(scip);
        assert(bestsol != NULL);

	sol = new CouNumber [nvars];

        // get solution values and objective
        SCIP_CALL_ABORT( SCIPgetSolVals(scip, bestsol, nvars, vars, sol) );
        obj = SCIPgetSolOrigObj(scip, bestsol);
     }
     else 
        sol = NULL;
     
     // release variables before freeing them
     for (int i=0; i<nvars; i++) {
        SCIP_CALL_ABORT( SCIPreleaseVar(scip, &vars[i]) );
     }

     // free memory
     SCIPfreeMemoryArray(scip, &vars);
     SCIP_CALL_ABORT( SCIPfree(&scip) );
   
     BMScheckEmptyMemory();     
  }
   else
#endif      
   {
      milp_ -> branchAndBound ();
      sol = CoinCopyOfArray (milp_ -> getColSolution (), problem_ -> nVars ());
      obj = milp_ -> getObjValue ();
   }

  return obj;
}

/// initialize MILP solvers if needed
void CouenneFeasPump::init_MILP () {

}
