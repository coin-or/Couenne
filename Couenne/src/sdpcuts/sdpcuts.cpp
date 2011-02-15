/* $Id$
 *
 * Name:    sdpcuts.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>

#include <sdpcuts.hpp>

#include <misc_util.hpp>
#include <OsiXxxSolverInterface.hpp>
#include <populate.hpp>
#include <CutGen.hpp>
#include <quadratic_cuts_check.hpp>
#include <linquad_cuts.hpp>
#include <orthocut.hpp>
#include <tracer.hpp>

#ifdef USE_INTERIOR_POINT
#include <OsiCpxSolverInterface.hpp>
#include <cplex.h>
#endif


#if 0
int main (int argc, const char **argv) {
	if (argc < 2) {
		printf("Missing argument [mps file]\n");
		return 1;
	}

	//print current configuration
	print_ifdefs(stdout);

	// determine problem name
	char name[256];
	char *name_pos = strrchr(const_cast <char *> (argv[1]), '/');
	if(name_pos != NULL)
		strcpy(name, &(name_pos[1]));
	else
		strcpy(name, argv[1]);
	char *last_dot_pos = strrchr(name, '.');
	if(last_dot_pos !=NULL)
		last_dot_pos[0] = '\0';

	// initialize random generator
	struct timeval tv;
	gettimeofday (&tv, NULL);
	srand48 (tv.tv_usec);

	Timer globaltimer;
	globaltimer.start();

#ifdef F_RES
	FILE *f_res = open_f_res();
#endif /* F_RES */
#ifdef SHORT_F_RES
	FILE *short_f_res = open_short_f_res();
#endif /* SHORT_F_RES */

	OsiXxxSolverInterface si;

#if (!defined TRACE_ALL)
	si.messageHandler()->setLogLevel(0);
#endif

	si.setDblParam(OsiPrimalTolerance,LP_TOLERANCE);
	si.setDblParam(OsiDualTolerance,LP_TOLERANCE);

	int n;			// number of x_i variables
	int t;			// number of y_i variables
	int origCons;		// number of original constraints
	double *b;		// obj coefficients for the x_i variables
	double *c;		// obj coefficients for the y_i variables
	double **Q;		// obj coefficients for the x_i_j variables
	double constant;	// obj function additive constant
	double **origMat;	// original problem matrix
	double *origRhs;	// original constraints' RHS
	char *origSense;	// original constraints' sense
	double *xlb;		// lower bounds on x_i vars
	double *xub;		// upper bounds on x_i vars
	double *ylb;		// lower bounds on y_i vars
	double *yub;		// upper bounds on y_i vars

	double objValue;

	int status;

	// read problem file
	status = populateProblem (argv[1],&n, &t, &origCons, &b, &c, &Q, &constant, &origMat, 
					&origRhs, &origSense, &xlb, &xub, &ylb, &yub,&si);

	if (status) 
		exit(1);


	Tracer *tracer = new Tracer(name,n,n*(n+1)/2,t);

#ifdef LINQUAD_BOUNDS_CUTS
	OsiCuts cs_diagbounds;
	linQuadCutGenOriginalBounds(xlb , xub , n , cs_diagbounds, tracer);
	si.applyCuts(cs_diagbounds);
#endif


	printf("Problem:                     %-20s\tx_card=%-3d Xext_card=%-6d y_card=%-3d\n",name,n,n*(n+1)/2,t);
	printf("Objective function constant: %.8f  (automatically added)\n",constant);

	CutGen  cg (n, t, origCons, constant, 
			b, c, (const double**) Q, 
			(const double**) origMat, origRhs, origSense, xlb, xub, ylb, yub, 
			&si, &globaltimer, tracer);

#ifdef MAX_CUTS_PER_ITER
	cg.set_max_nb_cuts(MAX_CUTS_PER_ITER);
#else
	cg.set_max_nb_cuts(100000);
#endif

#define OBJ_HISTORY_SIZE 10000
	double *obj_history = new double[OBJ_HISTORY_SIZE];
	for (int i=0;i<OBJ_HISTORY_SIZE;i++)
		obj_history[i] = COIN_DBL_MAX;

#ifdef DELETE_INACTIVE_CUTS
	int *del_row_ind = new int[1000000];
#ifdef FANCY_CUT_ACTIVITY
	int card_act = 0;
	int *activity = new int[1000000];
	int *gen_iter = new int[1000000];
#endif
#endif

	int niter = 0, ncuts = 0, tot_gen_cuts = 0, tot_del_cuts = 0;
	bool do_exit = false;


	printf("Initialization time:         %.2f\n",globaltimer.time());
	printf("\n");

	globaltimer.start();

#ifdef USE_INTERIOR_POINT
	si.XxxInitialSolveBaropt();
#else
	si.initialSolve();
#endif
	solver_status(&si);
	cg.updateSol();

	objValue = si.getObjValue() + constant;
	obj_history[0] = objValue;

#ifdef EXIT_IMPROVEMENT_ITER
	double *last_lp_val = new double [EXIT_IMPROVEMENT_ITER];
#ifdef EXIT_TAILING_OFF_VALUE
	double tailing_off_value = EXIT_TAILING_OFF_VALUE; 
	for(int i=0; i<EXIT_IMPROVEMENT_ITER; i++) {
		last_lp_val[i] = objValue + tailing_off_value;
	}
#endif
#ifdef EXIT_TAILING_OFF_PERC
	for(int i=0; i<EXIT_IMPROVEMENT_ITER; i++) {
		last_lp_val[i] = objValue + (fabs(objValue)*EXIT_TAILING_OFF_PERC/100);
	}	
#endif
#endif

#ifdef DELETE_INACTIVE_CUTS
	int del_init_nrows;
	del_init_nrows = si.getNumRows(); //never cut these ones (original and rlt (if separation of the latter is not allowed))
#endif
	
	print_current_sol(niter, globaltimer.time(), tot_gen_cuts, tot_gen_cuts-tot_del_cuts, 
				objValue, cg.bestObj(),cg.currObj());
#ifdef F_RES
	globaltimer.pause();
	print_file_current_sol(f_res, name, niter,
				globaltimer.time(), tot_gen_cuts, tot_gen_cuts-tot_del_cuts, 
				objValue, cg.bestObj(), cg.currObj());
	globaltimer.restore();
#endif

	tracer->setMainBound(objValue);
	tracer->setMainIterationTime(globaltimer.time());
	tracer->setMainLPTime(globaltimer.time());
	double time_at_previous_iteration = globaltimer.time();
	globaltimer.pause();

	//determining active rows
	const double *dualvars = si.getRowPrice();
	int active_rows = 0;
	for (int i=0;i<si.getNumRows();i++)
		if (fabs(dualvars[i]) < 2 * LP_TOLERANCE)
			active_rows++;
	tracer->setMainActiveCuts(active_rows);
	tracer->setMainAddedCuts(0);
	tracer->setMainTotalCuts(si.getNumRows());
	tracer->setMainDeletedCuts(0);
	tracer->setMainTotalEigendecompositions(0);

	tracer->setSDPNumNegativeEV(0);
	tracer->setSDPMostNegativeEV(0.0);
	tracer->setSDPCutsTime(0.0);
	tracer->setSDPCutsTotalCuts(0);


	tracer->setSparsifyTime(0.0);
	tracer->setSparsifyTotalCuts(0);
	tracer->setSparsifyDuplicatedCuts(0);
	tracer->setSparsifyWiseDecompositions(0);
	tracer->addSparsifyNz(0);
	tracer->addSparsifySingleColumnSparsity(0);
	tracer->addSparsifyColumnPairSparsity(0);
	tracer->addSparsifyTop20PercCutsViolation(0.0);

	tracer->setOrthocutTime(0.0);
	tracer->setOrthocutTotalCuts(0);

	// setLinquad[...] in linQuadCutGenOriginalBounds()

	tracer->setDisjunctiveCutsTime(0.0);
	tracer->setDisjunctiveCutsTotalCuts(0);

	// setHeuristics[...] in cg.updateSol()

	globaltimer.restore();

#ifdef CHECK_QUADRATIC_CUTS
	globaltimer.pause();
	QuadraticCuts qc(n,si.getColSolution(),tracer);
	globaltimer.restore();
#endif

	do {
		niter++;
		tot_del_cuts = 0;


#ifdef EXIT_ON_ITER
		if(niter > EXIT_ON_ITER) {
			do_exit = true;
			printf("Exit: Iteration limit [%d] exceeded\n",EXIT_ON_ITER);
		}
#endif

#ifdef EXIT_ON_TIME
		if(globaltimer.time() > EXIT_ON_TIME) {
			do_exit = true;
			printf("Exit: Time limit [%g] exceeded\n",EXIT_ON_TIME);
		}
#endif
		
		if (do_exit)
			break;

		tracer->newIter();

		cg.setIter(niter);

		OsiCuts cs;
		cg.generateCuts (si, cs);

		// the slow way to add cuts:
		// si.applyCuts (cs);

		// the fast way :)
		int size_cs = cs.sizeRowCuts();
		const OsiRowCut **newRowCuts = new const OsiRowCut * [size_cs];
		for(int i=0; i<size_cs; i++) {
			newRowCuts[i] = &cs.rowCut(i);
		}
		si.applyRowCuts(size_cs, newRowCuts);
		delete[] newRowCuts; 




		Timer lp_timer;
		lp_timer.start();
#ifdef USE_INTERIOR_POINT
		si.XxxResolveBaropt();
#else
		si.resolve();
#endif
		solver_status(&si);
		lp_timer.pause();

		objValue = si.getObjValue() + constant;
		cg.updateSol();

#ifdef WRITE_ITER_LP
		globaltimer.pause();

		si.writeLp("lastiter");
		globaltimer.restore();
#endif

#ifdef CHECK
		globaltimer.pause();
		const double *sol  = si.getColSolution ();
		int status = feasibility_check(n,t,origCons,sol,(const double**)origMat,(const double*)origRhs,(const char*)origSense,(const double*)xlb,(const double*)xub,(const double*)ylb,(const double*)yub);
		if (status)
			printf("ERROR: Infeasible Outer-approximation solution [status=%d]\n",status);

		// check that heuristic value is BELOW the current solution value
		if (cg.currObj() - LP_TOLERANCE > objValue)
			printf("ERROR: heuristic solution value > current solution value\n");
		globaltimer.restore();
#endif

		globaltimer.pause();

		tot_gen_cuts += ncuts = cs.sizeRowCuts ();

		print_current_sol(niter,globaltimer.time(), tot_gen_cuts, ncuts,
					objValue, cg.bestObj(), cg.currObj());
#ifdef F_RES
		print_file_current_sol(f_res, name, niter,globaltimer.time(), tot_gen_cuts, ncuts,
					objValue, cg.bestObj(), cg.currObj());
#endif
		tracer->setMainIterationTime(fabs(globaltimer.time() - time_at_previous_iteration));
		time_at_previous_iteration = globaltimer.time();
		tracer->setMainLPTime(lp_timer.time());
		tracer->setMainBound(objValue);
		//determining active rows
		const double *dualvars = si.getRowPrice();
		active_rows = 0;
		for (int i=0;i<si.getNumRows();i++)
			if (fabs(dualvars[i]) < 2 * LP_TOLERANCE)
				active_rows++;
		tracer->setMainActiveCuts(active_rows);
		tracer->setMainAddedCuts(cs.sizeRowCuts ());
		tracer->setMainTotalCuts(si.getNumRows());
		tracer->setMainDeletedCuts(0); // set later
		globaltimer.restore();


#ifdef CHECK_QUADRATIC_CUTS
		globaltimer.pause();
		qc.refresh(si.getColSolution());
		globaltimer.restore();
#endif

		if (cs.sizeRowCuts () == 0) {
			do_exit = true;
			printf("Exit: No cut generated\n");
		}

#ifdef EXIT_IMPROVEMENT_ITER
		int ind_iter = niter % EXIT_IMPROVEMENT_ITER;
#ifdef EXIT_TAILING_OFF_VALUE
		if ((last_lp_val[ind_iter] < objValue  + tailing_off_value) 
			&& (niter > EXIT_IMPROVEMENT_ITER)) {
			do_exit = true;
			printf("Exit: Solution improvement < %g in the last %d iterations\n"
				,EXIT_TAILING_OFF_VALUE,EXIT_IMPROVEMENT_ITER);
		} else
			last_lp_val[ind_iter] = objValue;
#endif
#ifdef EXIT_TAILING_OFF_PERC
		if((last_lp_val[ind_iter] - 
				(fabs(last_lp_val[ind_iter])*EXIT_TAILING_OFF_PERC/100) < objValue)
			&& (niter > EXIT_IMPROVEMENT_ITER)) {
			do_exit = true;
			printf("Exit: Solution improvement < %g%% in the last %d iterations\n"
				,EXIT_TAILING_OFF_PERC,EXIT_IMPROVEMENT_ITER);
		} else
			last_lp_val[ind_iter] = objValue;
#endif
#endif
		if (do_exit)
			break;

#ifdef DELETE_INACTIVE_CUTS
		// consider only cuts from the previous iteration
		int curr_nrows = si.getNumRows() - cs.sizeRowCuts();
		const double *y = si.getRowPrice();
#ifdef FANCY_CUT_ACTIVITY
		int card_old_activity = card_act;
		card_act = 0;
#endif
		int card_del_row_ind = 0;
		for(int del_i = del_init_nrows; del_i<curr_nrows; del_i++) {

#ifdef FANCY_CUT_ACTIVITY
			int old_act = 0, old_gen = niter;
			int old_pos = del_i - del_init_nrows;
		
			if(old_pos < card_old_activity) {
				old_act = activity[old_pos];
				old_gen = gen_iter[old_pos];
			}
				
			if(fabs(y[del_i]) < DELETE_INACTIVE_CUTS) {
				if(old_act > FANCY_CUT_ACTIVITY - 2) {
					del_row_ind[card_del_row_ind] = del_i;
					card_del_row_ind++;
				}
				else {
					activity[card_act] = old_act + 1;
					gen_iter[card_act] = old_gen;
					card_act++;
				}
			}
			else {
				activity[card_act] = 0;
				gen_iter[card_act] = old_gen;
				card_act++;
			}
#else
			if(fabs(y[del_i]) < DELETE_INACTIVE_CUTS) {
				del_row_ind[card_del_row_ind] = del_i;
				card_del_row_ind++;
			}
#endif
		}

#if (defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC) && (defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER)
		int historyIdx = (niter-CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER) % OBJ_HISTORY_SIZE;
		if (historyIdx < 0)
			historyIdx = 0;
		double oldObjValue = 
		  obj_history[historyIdx];

		if ((fabs(objValue - oldObjValue)*100/ fabs(oldObjValue))
				>= CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC)
		{
			si.deleteRows(card_del_row_ind, del_row_ind);
			tot_del_cuts += card_del_row_ind;
		}
#else
		si.deleteRows(card_del_row_ind, del_row_ind);
		tot_del_cuts += card_del_row_ind;
#endif
#ifdef CPLEX
		si.resolve();
#endif
#endif // DELETE_INACTIVE_CUTS

		tracer->setMainDeletedCuts(tot_del_cuts);

		obj_history[niter % OBJ_HISTORY_SIZE] = objValue;

	} while(do_exit == false);

#ifdef SHORT_F_RES
	print_file_short_sol(short_f_res,name,niter,globaltimer.time(),
			 tot_gen_cuts,tot_gen_cuts-tot_del_cuts,objValue,cg.bestObj());
#endif

	tracer->detailedReport();
	tracer->globalReport();



#ifdef F_RES
	fclose(f_res);
#endif

#ifdef SHORT_F_RES
	fclose(short_f_res);
#endif

#ifdef EXIT_IMPROVEMENT_ITER
	delete [] last_lp_val;
#endif
#ifdef DELETE_INACTIVE_CUTS
	delete [] del_row_ind;
#ifdef FANCY_CUT_ACTIVITY
	delete [] activity;
	delete [] gen_iter;
#endif
#endif
	delete [] obj_history;

	// freeing structures allocated by populate()
	delete [] b;
	delete [] c;
	for(int i=0;i<n;i++)
		delete [] Q[i];
	delete [] Q;
	for (int i=0;i<origCons;i++)
		delete [] origMat[i];
	delete [] origMat;
	delete [] origRhs;
	delete [] origSense;
	delete [] xlb;
	delete [] xub;
	delete [] ylb;
	delete [] yub;

	delete tracer;

	return 0;
}
#endif




/***********************************************************************/
void print_ifdefs(FILE *f) {
	fprintf(f, "Solver:        ");
#ifdef CPLEX
	fprintf(f, " CPLEX");
#endif
#ifdef CLP
	fprintf(f, " CLP");
#endif
#ifdef USE_INTERIOR_POINT
	fprintf(f, " USE_INTERIOR_POINT");
#endif
	

	fprintf(f, "\nTrace:         ");
#ifdef TRACE
	fprintf(f, " TRACE");
#endif
#ifdef TRACE_ALL
	fprintf(f, " TRACE_ALL");
#endif
#ifdef TRACE_EIGENVALUES
	fprintf(f, " TRACE_EIGENVALUES");
#endif
#ifdef TRACE_USED_EIGENVECTORS
	fprintf(f, " TRACE_USED_EIGENVECTORS");
#endif
#ifdef TRACE_RLT_MPS
	fprintf(f, " TRACE_RLT_MPS");
#endif
#ifdef TRACE_HEUR
	fprintf(f, " TRACE_HEUR");
#endif
#ifdef TRACE_CUT_TIME
	fprintf(f, " TRACE_CUT_TIME");
#endif
#ifdef TRACE_SPARSIFY_CUT_DEPTH
	fprintf(f, " TRACE_SPARSIFY_CUT_DEPTH");
#endif
#ifdef TRACE_DISJUNCTIVE_CUTS
	fprintf(f, " TRACE_DISJUNCTIVE_CUTS");
#endif
#ifdef TRACE_ORTHOCUT
	fprintf(f, " TRACE_ORTHOCUT");
#endif
#ifdef CHECK_QUADRATIC_CUTS
	fprintf(f, " CHECK_QUADRATIC_CUTS");
#endif
#ifdef F_RES
	fprintf(f, " F_RES");
#endif
#ifdef SHORT_F_RES
	fprintf(f, " SHORT_F_RES");
#endif



	fprintf(f, "\nExit:          ");	
#ifdef EXIT_ON_ITER
	fprintf(f, " EXIT_ON_ITER=%d",EXIT_ON_ITER);
#endif
#ifdef EXIT_ON_TIME
	fprintf(f, " EXIT_ON_TIME=%.2f",EXIT_ON_TIME);
#endif
#if (!defined EXIT_ON_ITER) && (!defined EXIT_ON_TIME)
#error EXIT_ON_TIME or EXIT_ON_ITER must be set
#endif
#ifdef EXIT_IMPROVEMENT_ITER
	fprintf(f, " EXIT_IMPROVEMENT_ITER=%d",EXIT_IMPROVEMENT_ITER);
#ifdef EXIT_TAILING_OFF_VALUE
	fprintf(f, " EXIT_TAILING_OFF_VALUE=%e",EXIT_TAILING_OFF_VALUE);
#endif
#ifdef EXIT_TAILING_OFF_PERC
	fprintf(f, " EXIT_TAILING_OFF_PERC=%g%%",(double)EXIT_TAILING_OFF_PERC);
#endif
#if (!defined EXIT_TAILING_OFF_VALUE) && (!defined EXIT_TAILING_OFF_PERC) 
#error EXIT_TAILING_OFF_VALUE or EXIT_TAILING_OFF_PERC must be set
#endif
#if (defined EXIT_TAILING_OFF_VALUE) && (defined EXIT_TAILING_OFF_PERC) 
#error Only one between EXIT_TAILING_OFF_VALUE or EXIT_TAILING_OFF_PERC can be defined
#endif
#endif

	fprintf(f, "\nAlg. options:  ");
#ifdef RLT_CUTS
	fprintf(f, " RLT_CUTS");
#endif
#ifdef RLT_SEPARATION
#ifndef RLT_CUTS
#error RLT_CUTS must be defined in order to have RLT_SEPARATION
#endif
	fprintf(f, " RLT_SEPARATION");
#endif
#ifdef INCL_ORIG
	fprintf(f, " INCL_ORIG");
#endif
#ifdef SCALE_EIGENV
	fprintf(f, " SCALE_EIGENV");
#endif
#ifdef ONLY_NEG_EIGENV
	fprintf(f, " ONLY_NEG_EIGENV");
#endif
#ifdef ONLY_MOST_NEG
	fprintf(f, " ONLY_MOST_NEG");
#ifdef ONLY_NEG_EIGENV
#error ONLY_NEG_EIGENV and ONLY_MOST_NEG conflicts
#endif
#endif
#ifdef SPARSIFY2
	fprintf(f, " SPARSIFY2");
#endif
#ifdef SPARSIFY
	fprintf(f, " SPARSIFY");
#endif
#ifdef SPARSIFY_MINOR_SDP_CUTS
	fprintf(f," SPARSIFY_MINOR_SDP_CUTS");
#endif
#ifdef WISE_SPARSIFY
	fprintf(f, " WISE_SPARSIFY");
#endif
#ifdef SPARSIFY_REMOVE_DUPLICATES
	fprintf(f, " SPARSIFY_REMOVE_DUPLICATES");
#endif
#ifdef NORMALIZE_SPARSE_CUTS
	fprintf(f, " NORMALIZE_SPARSE_CUTS");
#endif
#ifdef SLIDING
	fprintf(f, " SLIDING");
#endif
#ifdef ORTHOCUT
	fprintf(f, " ORTHOCUT");
#endif
#ifdef DISJUNCTIVE_CUTS
	fprintf(f, " DISJUNCTIVE_CUTS");
#endif
#ifdef ZERO_LB_DIAGONAL
	fprintf(f, " ZERO_LB_DIAGONAL");
#endif
#ifdef LINQUAD_CUTS
	fprintf(f, " LINQUAD_CUTS");
#endif
#ifdef LINQUAD_BOUNDS_CUTS
	fprintf(f, " LINQUAD_BOUNDS_CUTS");
#endif


	fprintf(f, "\nCut management:");
#ifdef MAX_CUTS_PER_ITER
	fprintf(f, " MAX_CUTS_PER_ITER=%d",MAX_CUTS_PER_ITER);
#endif
#ifdef RAND_CUT_ADD
	fprintf(f, " RAND_CUT_ADD=%.2f",RAND_CUT_ADD);
#endif
#ifdef DELETE_INACTIVE_CUTS
	fprintf(f, " DELETE_INACTIVE_CUTS(y_i<%e)",DELETE_INACTIVE_CUTS);
#endif
#ifdef FANCY_CUT_ACTIVITY
#if (!defined DELETE_INACTIVE_CUTS)
#error FANCY_CUT_ACTIVITY must be defined only with DELETE_INACTIVE_CUTS
#endif
	fprintf(f, " FANCY_CUT_ACTIVITY=%d",FANCY_CUT_ACTIVITY);
#endif

#if (defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC) && (!defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER)
#error CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC must be defined together with CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER
#endif
#if (!defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC) && (defined CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER)
#error CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER must be defined together with CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC
#endif
#ifdef CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC
#if (!defined DELETE_INACTIVE_CUTS)
#error CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC & CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER must be defined only with DELETE_INACTIVE_CUTS
#endif
	fprintf(f, " CUT_DELETION_OBJ_PERC_IMPROVEMENT=(%.4f,%d)",
	CUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC,CUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER);
#endif
#ifdef ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS
	fprintf(f, " ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS=%.2f%%",ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS);
#endif
#ifdef ADD_RAND_PERC_SPARSIFY_CUTS
	fprintf(f, " ADD_RAND_PERC_SPARSIFY_CUTS=%.2f%%",ADD_RAND_PERC_SPARSIFY_CUTS);
#endif
#if (defined ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS) && (defined ADD_RAND_PERC_SPARSIFY_CUTS)
#error either ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS or ADD_RAND_PERC_SPARSIFY_CUTS can be defined
#endif

	fprintf(f, "\nHeuristics:    ");
#ifdef HEUR_XXT
	fprintf(f, " HEUR_XXT");
#endif
#ifdef HEUR_MIN_MATRIX_NORM
	fprintf(f, " HEUR_MIN_MATRIX_NORM");
#endif
#ifdef HEUR_IMPROVE_SOLUTION
	fprintf(f, " HEUR_IMPROVE_SOLUTION");
#endif

	fprintf(f, "\nDebug:         ");
#ifdef CHECK
	fprintf(f, " CHECK");
#endif
#ifdef WRITE_ITER_LP
	fprintf(f, " WRITE_ITER_LP");
#endif


	fprintf(f, "\n\n");
	fflush(f);
} /* print_ifdefs */

/***********************************************************************/
void solver_status(OsiSolverInterface *solver) {

	if (solver->isAbandoned()) {
		printf("### ERROR: Numerical difficulties in Solver\n");
		exit(1);
	}
	
	if (solver->isProvenPrimalInfeasible()) {
		printf("### WARNING: Problem is infeasible\n");
	
#ifdef TRACE_ALL
		solver->writeLp("xxxinfeas.lp");
		printf("Problem written in xxxinfeas.lp\n");
#endif
		exit(1);
	}
} /* solver_status */
/***********************************************************************/
FILE *open_f_res() {
	FILE *f_res = NULL;
	f_res = fopen("f_res.xxx", "r");

	if(f_res == NULL) {
		f_res = fopen("f_res.xxx", "w");
		print_ifdefs(f_res);
	}
	else {
		fclose(f_res);
		f_res = fopen("f_res.xxx", "a");
	}
	return f_res;
} /* open_f_res() */
/***********************************************************************/
FILE *open_short_f_res() {
	FILE *short_f_res = NULL;
	short_f_res = fopen("short_f_res.xxx", "r");

	if(short_f_res == NULL) {
		short_f_res = fopen("short_f_res.xxx", "w");
		print_ifdefs(short_f_res);
	}
	else {
		fclose(short_f_res);
		short_f_res = fopen("short_f_res.xxx", "a");
	}
	return short_f_res;
} /* open_short_f_res() */

/***********************************************************************/
void print_current_sol(int niter, double time, int tot_gen_cuts, int current_cuts, double objValue,double bestHeurObj, double currHeurObj) {

	if (niter == 0) 
		printf("    %4s %8s %5s %5s %12s %12s %12s\n",
			"iter","time","tcuts","gcuts","upper bound","curr heur","best heur");
	printf ("::: %4d %8.2f %5d %5d %12.4f ", 
		niter, time, tot_gen_cuts, current_cuts, objValue);
	if (currHeurObj > -1e40)
		printf("%12.4f ",currHeurObj);
	else
		printf("%12s ","N/A");

	if (bestHeurObj > -1e40)
		printf("%12.4f ",bestHeurObj);
	else
		printf("%12s ","N/A");
	printf("\n");
} /* print_current_sol() */
/***********************************************************************/
void print_file_current_sol(FILE* file, char* name, int niter, double time, int tot_gen_cuts, int current_gen_cuts, double objValue,double bestHeurObj, double currHeurObj) {
	if(niter == 0)
		fprintf(file, "%-16s ", name);
	else
		fprintf(file, "                 ");

	fprintf(file, "%12.4f %10d %10.2f %10d %10d ", 
		objValue, niter, time, tot_gen_cuts, current_gen_cuts);

	if (currHeurObj > -1e40)
		fprintf(file,"%12.4f",currHeurObj);
	else
		fprintf(file,"%12s ","N/A");

	if (bestHeurObj > -1e40)
		fprintf(file,"%12.4f ",bestHeurObj);
	else
		fprintf(file,"%12s ","N/A");
	fprintf(file,"\n");
	fflush(file);
} /* print_file_current_sol() */
/***********************************************************************/
void print_file_short_sol(FILE* file,char* name, int niter, double time, int tot_gen_cuts, int current_gen_cuts, double objValue,double bestHeurObj) {

	fprintf(file,"%-16s %12.4f ",
		name,
		objValue);
	if (bestHeurObj > -1e40)
		fprintf(file,"%12.4f ",bestHeurObj);
	else
		fprintf(file,"%12s ","N/A");
	fprintf(file,"%10.2f %10d %10d %10d\n",time,niter,tot_gen_cuts,current_gen_cuts);
	fflush(file);
} /* print_file_short_sol() */
/***********************************************************************/
int feasibility_check(const int n, const int t, const int cons, const double *sol, const double **origMat , const double *origRhs, const char *origSense, const double *xlb,const double *xub,const double *ylb,const double *yub) {
	int N =  n*(n+3)/2;;
	bool violation = false;
	bool violation_hard = false;
/*
for(int i=0;i<n;i++)
	printf("x_%d=%.2f ",i,sol[i]);
for(int i=n;i<N;i++)
	printf("X_%d=%.2f ",i,sol[i]);
for(int i=N;i<N+t;i++)
	printf("y_%d=%.2f ",i,sol[i]);
printf("\n");
*/

	//checking bound feasibility
	for(int i=0;i<n;i++) {
		if (( sol[i] < xlb[i] - LP_TOLERANCE )||( sol[i] > xub[i] + LP_TOLERANCE )) {
			printf("feasibility_check:: x_%d=%.8f violates bounds [%.8f,%.8f]\n",i,sol[i],xlb[i],xub[i]);
		}
	}
	for(int i=0;i<t;i++) {
		if (( sol[N+i] < ylb[i] - LP_TOLERANCE )||( sol[N+i] > yub[i] + LP_TOLERANCE )) {
			printf("feasibility_check:: y_%d=%.8f violates bounds [%.8f,%.8f]\n",i,sol[N+i],ylb[i],yub[i]);
			return FEAS_CHECK_BOUNDS_VIOLATION;
		}
	}
	
	//checking constraint feasibility
	for (int k=0;k<cons;k++) {
		double lhs = 0.0; 
		bool y_terms = false;
		for (int i=0;i<N;i++) {
			lhs += origMat[k][i] * sol[i];
		}
		for (int i=N;i<N+t;i++) {
			lhs += origMat[k][i] * sol[i];
			if (origMat[k][i] != 0.0)
				y_terms = true;
		}


/*
printf("\n");
printf("origMat cons %d: ",k);
for(int i=0;i<n;i++)
	printf("origMat[%d][%d]x=%.2f ",k,i,origMat[k][i]);
printf("\n");
for(int i=n;i<N;i++)
	printf("origMat[%d][%d]X=%.2f ",k,i,origMat[k][i]);
printf("\n");
for(int i=N;i<N+t;i++)
	printf("origMat[%d][%d]y=%.2f ",k,i,origMat[k][i]);
printf("\n");
printf("lhs=%.2f rhs=%.2f sense=%c\n",lhs,origRhs[k],origSense[k]);
*/
		switch(origSense[k]) {
			case 'L': 
				if ((lhs - LP_TOLERANCE) > origRhs[k]) {
					violation = true;
					if (y_terms == false)
						violation_hard = true;
#ifdef TRACE_ALL
					printf("feasibility_check:: violation in the original 'L' constraint %d     [%d]\n",k,violation_hard);
#endif
				}
				break;
			case 'G': 
				if ((lhs + LP_TOLERANCE) < origRhs[k]) {
					violation = true;
					if (y_terms == false)
						violation_hard = true;
#ifdef TRACE_ALL
					printf("feasibility_check:: violation in the original 'G' constraint %d     [%d]\n",k,violation_hard); 
#endif
				}
				break;
			case 'E': 
				if (((lhs + LP_TOLERANCE) < origRhs[k]) || ((lhs - LP_TOLERANCE) > origRhs[k])) {
					violation = true;
					if (y_terms == false)
						violation_hard = true;
#ifdef TRACE_ALL
					printf("feasibility_check:: violation in the original 'E' constraint %d     [%d]\n",k,violation_hard);
#endif
				}
				break;
			default: 
				printf("feasibility_check:: unrecognized constraint sense '%c' (constraint index = %d)\n",origSense[k],k);
				return -999;
		}
	}

	if (violation == true)
		if (violation_hard == true)
			return FEAS_CHECK_CONSTRAINT_VIOLATION_NO_RECOVER;
		else
			return FEAS_CHECK_CONSTRAINT_VIOLATION;
	else 
		return FEAS_CHECK_NO_VIOLATION;
} //easibility_check()

/***********************************************************************/
double evaluateSolution(const int n, const int t,const double *heurSol, const double *b, const double *c , const double **Q, const double objConstant) {
	int N = n*(n+3)/2;
	double heurSolValue = 0.0;
	for(int i=0;i<n;i++) {
		heurSolValue += heurSol[i] * b[i];
	}
	for(int i=0;i<t;i++) {
		heurSolValue += heurSol[i+N] * c[i];
	}
	for(int i=0;i<n;i++)
		for(int j=i;j<n;j++) {
			heurSolValue += heurSol[indexQ(i,j,n)] * Q[i][j];
			heurSolValue += heurSol[indexQ(i,j,n)] * Q[j][i];
		}

	heurSolValue += objConstant;

	return heurSolValue;
} //evaluateSolution()
/***********************************************************************/
