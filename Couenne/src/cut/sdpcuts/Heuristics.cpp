/* $Id$
 *
 * Name:    Heuristics.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>

#include <Heuristics.hpp>
#include <CoinPackedVector.hpp>
#include <OsiSolverInterface.hpp>
#include <OsiXxxSolverInterface.hpp>
#include <sdpcuts.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

// constructor
Heuristics::Heuristics(	const int n, 
			const int t,
			const int cons,
			const double objConst,
			const double *b,
			const double *c,
			const double **Q,
			const double **origMat,
			const double *origRhs,
			const char *origSense,
			const double *xlb,
			const double *xub,
			const double *ylb,
			const double *yub,
			const OsiSolverInterface *si,
			Tracer *tracer
			):

	n_ (n),
	t_ (t),
	cons_ (cons),
	objConst_ (objConst),
	b_ (b),
	c_ (c),
	Q_ (Q),
	origMat_ (origMat),
	origRhs_ (origRhs),
	origSense_ (origSense),
	xlb_ (xlb),
	xub_ (xub),
	ylb_ (ylb),
	yub_ (yub),
	si_ (si),
	currObj_ (-DBL_MAX),
	bestObj_ (-DBL_MAX),
	tracer_ (tracer)
	{

	N_ = n*(n+3)/2;

	bestSol_ = new double[N_+t_];
	xxTSol_  = new double[N_+t_];
	MNSol_   = new double[N_+t_];
	for (int i=0;i<N_+t_;i++) {
		bestSol_[i] = -DBL_MAX;
		xxTSol_ [i] = -DBL_MAX;
		MNSol_  [i] = -DBL_MAX;
	}



	//initialize heurLP solver interface
	heurLPimproveSi_.messageHandler()->setLogLevel(0);
	// maximization problem
	heurLPimproveSi_.setObjSense(-1);

	heurLbRowAdded_ = new bool[cons_];

	for(int i=0;i<cons_;i++)
		heurLbRowAdded_[i] = false;
	for (int i=0;i<t;i++)
		heurLPimproveSi_.addCol (0, NULL, NULL,  ylb[i], yub[i], c[i]);
	for (int k=0;k<cons;k++) {
		double rhs = origRhs[k];
		CoinPackedVector vector;
		vector.clear();
		for (int j=0;j<t;j++) {
			double val = origMat[k][N_+j];
			if (val != 0.0)
				vector.insert(j,val);
		}
		if (vector.getNumElements() > 0) {
			heurLPimproveSi_.addRow(vector,origSense[k],rhs,0);
			heurLbRowAdded_[k] = true;
		}
	}
	heurLPimproveSi_.initialSolve(); // we solve via primal simplex to permit a warm start later on with resolve();



	//initialize min norm LP Heur solver interface

	double *coeff = const_cast <double *> ( si_->getObjCoefficients() );

	temp_row_ = new double[N_+t];  
	for(int i=0;i<N_+t;i++)
		temp_row_[i] = 0.0;
	MNLPSi_.messageHandler()->setLogLevel(0);
	MNLPSi_.setObjSense(1); //minimization
	//add variables in MNLPSi
	for(int i=0;i<n;i++) {
		MNLPSi_.addCol (0, NULL, NULL,  xlb[i], xub[i], 0.0); //x_i
	}
	for(int i=0;i<t;i++) {
		MNLPSi_.addCol (0, NULL, NULL,  ylb[i], yub[i], 0.0); //y_i
	}
	for(int i=0;i<n;i++) {
		double objcoeff;
		objcoeff = fabs(coeff[i]);
		MNLPSi_.addCol (0, NULL, NULL,  0.0, MNLPSi_.getInfinity(), objcoeff); //r^+_f(i,0)
		MNLPSi_.addCol (0, NULL, NULL,  0.0, MNLPSi_.getInfinity(), objcoeff); //r^-_f(i,0)
	}
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++) {
			double objcoeff;
			if (i==j)
 				objcoeff = fabs(coeff[indexQ(i,i,n_)]);
			else
				objcoeff = fabs(coeff[indexQ(i,j,n_)]);
			MNLPSi_.addCol (0, NULL, NULL,  0.0, MNLPSi_.getInfinity(), objcoeff); //r^+_f(i,j)
			MNLPSi_.addCol (0, NULL, NULL,  0.0, MNLPSi_.getInfinity(), objcoeff); //r^-_f(i,j)
		}
	}

	//add constraints in MNLPSi
	for (int k=0;k<cons;k++) {
		// original constraints
		CoinPackedVector vector;
		vector.clear();
		MNLPSi_.addRow(vector,origSense[k],origRhs[k],0);
	}
	int idx;
	idx = n_+t_;
	for (int i=0;i<n;i++) {
		// constraints x_i +r^+_f(i,0) -r^-_f(i,0)= x_i
		CoinPackedVector vector;
		vector.clear();
		vector.insert(i,1.0); //x_i
		vector.insert(idx,1.0);	//r+
		idx++;
		vector.insert(idx,-1.0); //r-
		idx++;
		MNLPSi_.addRow(vector,'E',0.0,0);
	}
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++) {
			// constraints appr.X_ij +r^+_f(i,j) -r^-_f(i,j)= X_ij
			CoinPackedVector vector;
			vector.clear();
			vector.insert(idx,1.0); //r+
			idx++;
			vector.insert(idx,-1.0); //r-
			idx++;
			MNLPSi_.addRow(vector,'E',0.0,0);
		}
	}
	MNLPSi_.initialSolve(); // we solve via primal simplex to permit a warm start later on with resolve();
} // Heuristics()
/***********************************************************************/
// destructor
Heuristics::~Heuristics() {
	delete [] heurLbRowAdded_;
	delete [] bestSol_;
	delete [] xxTSol_;
	delete [] MNSol_;
	delete [] temp_row_;
} // ~Heuristics()
/***********************************************************************/
int Heuristics::run() {
	Timer heur_timer;
	heur_timer.start();

	int status;
	bool feasible_solution = false;
	double *tmpSol;

	currObj_ = -DBL_MAX;

	// xx^T heuristic solution


#ifdef HEUR_XXT
	Timer xxtheur_timer;
	xxtheur_timer.start();
	double xxtheurvalue = -DBL_MAX;
	double xxt_heurlp_value = -DBL_MAX;

	tmpSol = xxTHeur();
#ifdef HEUR_IMPROVE_SOLUTION
	status = processSol(tmpSol,true,&xxtheurvalue,&xxt_heurlp_value);
#else
	status = processSol(tmpSol,false,&xxtheurvalue,&xxt_heurlp_value);
#endif
	if (status == 0)
		feasible_solution = true;

	tracer_->setHeuristicsxxTTime(xxtheur_timer.time());
	if (xxtheurvalue > -DBL_MAX) {
		tracer_->setHeuristicsxxTSolution(xxtheurvalue);
	}
	if (xxt_heurlp_value > -DBL_MAX) {
		tracer_->setHeuristicsxxTSolutionLPHeuristicImprovement
			(fabs(xxt_heurlp_value-xxtheurvalue));
	}
#endif //HEUR_XXT



	// matrix norm minimization heuristic solution
#ifdef HEUR_MIN_MATRIX_NORM
	Timer MNheur_timer;
	MNheur_timer.start();
	double MNheurvalue = -DBL_MAX;
	double MN_heurlp_value = -DBL_MAX;

	tmpSol = MNHeur();
#ifdef HEUR_IMPROVE_SOLUTION
	status = processSol(tmpSol,true,&MNheurvalue,&MN_heurlp_value);
#else
	status = processSol(tmpSol,false,&MNheurvalue,&MN_heurlp_value);
#endif
	if (status == 0)
		feasible_solution = true;

	tracer_->setHeuristicsMNLPTime(MNheur_timer.time());
	if (MNheurvalue > -DBL_MAX)
		tracer_->setHeuristicsMNLPSolution(MNheurvalue);
	if (MN_heurlp_value > -DBL_MAX)
		tracer_->setHeuristicsMNLPSolutionLPHeuristicImprovement
				(fabs(MN_heurlp_value-MNheurvalue));
#endif //HEUR_MIN_MATRIX_NORM

	tracer_->setHeuristicsTime(heur_timer.time());
	if (bestObj_ > -DBL_MAX) {
		tracer_->setHeuristicsBestSolution(bestObj_);
	}
	if (feasible_solution) {
		tracer_->setHeuristicsCurrentSolution(currObj_);
		return 0;
	}
	else {
		return 1;
	}
} // run()
/***********************************************************************/
double* Heuristics::xxTHeur() {

	double *sol = const_cast <double *> ( si_->getColSolution () );
	double *x,*y;
	x = sol;
	if (t_ > 0)
		y = &sol[N_];

	// compute xx^T
	for(int i=0;i<n_;i++)
		xxTSol_[i] = x[i];
	for(int i=0;i<t_;i++)
		xxTSol_[i+N_] = y[i];
	for(int i=0;i<n_;i++)
		for(int j=i;j<n_;j++)
			xxTSol_[indexQ(i,j,n_)] = x[i] * x[j];
	return xxTSol_;
} // xxTHeur()
/***********************************************************************/
double* Heuristics::MNHeur() {
	double *sol = const_cast <double *> ( si_->getColSolution () );

	double *x,*y;
	x = sol;
	if (t_ > 0)
		y = &sol[N_];

	//update original constraints
	for(int k=0;k<cons_;k++) {
		for (int i=0;i<n_+t_;i++) //clear temp_row_
			temp_row_[i] = 0.0;
		for (int i=0;i<n_;i++)
			temp_row_[i] = origMat_[k][i]; //coeff of x_i
		for (int i=0;i<t_;i++)
			temp_row_[N_+i] = origMat_[k][N_+i]; //coeff of y_i
		for (int i=0;i<n_;i++)
			for (int j=i;j<n_;j++) {
				double value = sqrt(fabs(x[indexQ(i,j,n_)]))/2.0;
				temp_row_[i] += value;
				temp_row_[j] += value;
			}

		for(int i=0;i<N_+t_;i++)
			MNLPSi_.XxxModifyCoefficient(k,i,temp_row_[i]);
	}
	//update constraints x_i +r^+_f(i,0) -r^-_f(i,0)= x_i

	int rowidx;
	rowidx = cons_;
	for (int i=0;i<n_;i++) {
		MNLPSi_.setRowBounds(rowidx+i,x[i],x[i]);
	}
	rowidx = cons_ + n_;
	for (int i=0;i<n_;i++) 
		for (int j=i;j<n_;j++) {
			if (i==j) {
				double value = sqrt(fabs(x[indexQ(i,j,n_)]));
				MNLPSi_.XxxModifyCoefficient(rowidx,i,value);
			} else { //i != j
				double value = sqrt(fabs(x[indexQ(i,j,n_)]))/2.0;
				MNLPSi_.XxxModifyCoefficient(rowidx,i,value);
				MNLPSi_.XxxModifyCoefficient(rowidx,j,value);
			}
			MNLPSi_.setRowBounds(rowidx,x[indexQ(i,j,n_)],x[indexQ(i,j,n_)]);
			rowidx++;
		}

	MNLPSi_.resolve();
	if (MNLPSi_.isAbandoned()) {
		printf("Heur:: HeurLP: Numerical difficulties in Solver\n");
		return NULL;
	}
	if (MNLPSi_.isProvenPrimalInfeasible()) {evaluateSolution(n_,t_,MNSol_,b_,c_,Q_,objConst_);
#if (defined TRACE_ALL)||(defined TRACE_HEUR)
		printf("Heur:: MNLPHeur Problem is infeasible\n");
#endif
		return NULL;
	}
	
	const double *lpheurSol = MNLPSi_.getColSolution();

	// compute xx^T
	for(int i=0;i<n_;i++)
		MNSol_[i] = lpheurSol[i];
	for(int i=0;i<t_;i++)
		MNSol_[N_+i] = lpheurSol[n_+i];
	for(int i=0;i<n_;i++)
		for(int j=i;j<n_;j++)
			MNSol_[indexQ(i,j,n_)] = MNSol_[i] * MNSol_[j];



/*
printf("Current Solution:\n");
for (int i=0;i<n_;i++)
	printf("x_%d=%.2f ",i,x[i]);
printf("\n");
for (int i=0;i<n_;i++)
	for (int j=i;j<n_;j++)
		printf("x_%d_%d=%.2f ",i,j,x[indexQ(i,j,n_)]);
if (t_>0) {
printf("\n");
for (int i=0;i<t_;i++)
	printf("y_%d=%.2f ",i,y[i]);
}
printf("\n");
int status = feasibility_check(n_,t_,cons_,MNSol_,origMat_,origRhs_,origSense_,xlb_,xub_,ylb_,yub_);
double value = evaluateSolution(n_,t_,MNSol_,b_,c_,Q_,objConst_);
printf("---->feas:%d  value=%.2f\n",status,value);

printf("curr sol value=%.2f\n",evaluateSolution(n_,t_,const_cast <double *>(si_->getColSolution ()),b_,c_,Q_,objConst_));

for (int i=0;i<N_+t_;i++)
	printf("%.2f ",const_cast <double *>(si_->getColSolution ())[i]);
printf("\n");
for (int i=0;i<n_;i++)
	printf("%.2f ",b_[i]);
for (int i=0;i<n_;i++)
	for (int j=i;j<n_;j++)
		printf("%.2f ",Q_[i][j]);
for (int i=0;i<t_;i++)
	printf("%.2f ",c_[i]);
printf("\n");
MNLPSi_.writeLp("test", "lp",1e-5,10,5,0.0,true);

exit(9);
*/
/*
int status = feasibility_check(n_,t_,cons_,MNSol_,origMat_,origRhs_,origSense_,xlb_,xub_,ylb_,yub_);
double value = evaluateSolution(n_,t_,MNSol_,b_,c_,Q_,objConst_);
printf("---->feas:%d  value=%.2f\n",status,value);
*/
	return MNSol_;
} //MNHeur()
/***********************************************************************/
int Heuristics::processSol(double* sol, bool improveHeurLP, double *origvalue, double *lpheurvalue) {
	int status;
	double value   = -DBL_MAX;
	
	if (sol == NULL)
		return -1;
	status = feasibility_check(n_,t_,cons_,sol,origMat_,origRhs_,origSense_,xlb_,xub_,ylb_,yub_);
	if (status == FEAS_CHECK_NO_VIOLATION)	
		value = evaluateSolution(n_,t_,sol,b_,c_,Q_,objConst_);
	(*origvalue) = value;
#if (defined TRACE_ALL)||(defined TRACE_HEUR)
	if (value > -DBL_MAX)
		printf("Heur:: Heuristic solution value = %.3f",value);
	else 
		printf("Heur:: Heuristic solution infeasible. ");
	switch(status) {
		case FEAS_CHECK_NO_VIOLATION: 
			printf("(feasible)\n");
			break;
		case FEAS_CHECK_CONSTRAINT_VIOLATION:
			printf("Enforcing feasibility...\n");
			break;
		case FEAS_CHECK_CONSTRAINT_VIOLATION_NO_RECOVER:
			printf("Cannot enforce feasibility!\n");
			break;
		case FEAS_CHECK_BOUNDS_VIOLATION:
			printf("\nERROR: variable bounds are violated!!!\n"); //should never happen
			break;
		default:
			printf("\nERROR: feasibility_check() return code unknown!!\n");
	}
#endif
#ifdef TRACE_ALL
	printf("Heur:: Heuristic Solution:\nx = ");
	for(int i=0;i<n_;i++)
		printf("%.2f ",sol[i]);
	if (t_>0) {
		printf("\ny = ");
		for(int i=0;i<t_;i++)
			printf("%.2f ",sol[i+N_]);
		printf("\n");
	}
#endif
	if ( (status != FEAS_CHECK_NO_VIOLATION) && (status != FEAS_CHECK_CONSTRAINT_VIOLATION) )
		return status;

	if (status == FEAS_CHECK_NO_VIOLATION)
		update(sol,value);

	// by changing the y_i variables via LP we can try to
	//	-improve the current solution heurSol if it is feasible
	//	-make the current solution heurSol feasible if it is infeasible and all the violated
	//	 constraints contain a a variable y_i
	// if there are no y_i variables (t_=0), heurLP is useless
	if ((improveHeurLP) && (t_>0) && (heurLP_improveSolution(sol) == 0)) {
		status = feasibility_check(n_,t_,cons_,sol,origMat_,origRhs_,origSense_,xlb_,xub_,ylb_,yub_);
		double newValue = evaluateSolution(n_,t_,sol,b_,c_,Q_,objConst_);
		if (status == FEAS_CHECK_NO_VIOLATION) {
			(*lpheurvalue) = newValue; 
			update(sol,newValue);
#if (defined TRACE_ALL)||(defined TRACE_HEUR)
			if (newValue > value)
				printf("Heur:: HeurLP improved feasible solution found! value=%.2f\n",newValue);
			else
				printf("Heur:: HeurLP could not improve the current heuristic solution\n");
#endif
		}
	}
	return status;
	
} // processSol()
/***********************************************************************/
int Heuristics::heurLP_improveSolution(double* sol) {
	int rowcnt = 0;
	for (int k=0;k<cons_;k++) {
		if (heurLbRowAdded_[k]) {
			double rhs = origRhs_[k];
			for (int i=0;i<N_;i++) {
				rhs -= origMat_[k][i] * sol[i];
			}
			if (origSense_[k] == 'G')
				heurLPimproveSi_.setRowLower(rowcnt,rhs);
			if (origSense_[k] == 'L')
				heurLPimproveSi_.setRowUpper(rowcnt,rhs);
			if (origSense_[k] == 'E') {
				heurLPimproveSi_.setRowLower(rowcnt,rhs);
				heurLPimproveSi_.setRowUpper(rowcnt,rhs);
			}
			rowcnt++;
		}
	}

	heurLPimproveSi_.resolve();

	if (heurLPimproveSi_.isAbandoned()) {
		printf("Heur:: HeurLP: Numerical difficulties in Solver\n");
		return -1;
	}
	if (heurLPimproveSi_.isProvenPrimalInfeasible()) {
#if (defined TRACE_ALL)||(defined TRACE_HEUR)
		printf("Heur:: HeurLP Problem is infeasible\n");
#endif
		return 1;
	}

	const double *lpheurSol = heurLPimproveSi_.getColSolution();
	for(int i=0;i<t_;i++)
		sol[i+N_] = lpheurSol[i];

	return 0;
} // improveSolution()
/***********************************************************************/
int Heuristics::update(double* sol, double value) {
	if (currObj_ < value)
		currObj_ = value;
	if (bestObj_ < value) {
		bestObj_ = value;
		for(int i=0;i<N_+t_;i++)
			bestSol_[i] = sol[i];
		return 0;
	} else
		return 1;
}
/***********************************************************************/



