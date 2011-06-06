/* $Id$
 *
 * Name:    sedumi.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>

#include <misc_util.hpp>
#include <OsiXxxSolverInterface.hpp>
#include <populate.hpp>

#include <CoinPackedVector.hpp>



//#define DEBUG_ORIGINALPROBLEM
//#define DEBUG_ORIGINALPROBLEM_ARRAYS
//#define DEBUG_SDPMODEL
//#define DEBUG_SDPMODEL_ARRAYS
//#define DEBUG_SDPMODEL_ARRAYS_FULLMATRIX

//#define ADD_LB
//#define ADD_UB
//#define ADD_RLT




#if 0
double getCoinPackedVectorElementAt(const CoinPackedVector vector, int index){
	int size = vector.getNumElements();
	const int *indices = vector.getIndices();
	const double *elements = vector.getElements();
	for (int i=0;i<size;i++)
		if (indices[i] == index)
			return elements[i];
	return 0.0;
}
#endif

#if 0
int main (int argc, const char **argv) {
	if (argc < 2) {
		printf("Missing argument [mps file] <[matlab file]>\n");
		return 1;
	}

	FILE *matlabout;
	if (argc >= 3)
		matlabout = fopen(argv[2],"w");
	else
		matlabout = stdout;

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

	OsiXxxSolverInterface si;

	Timer globaltimer;
	globaltimer.start();

	int n;			// number of x_i variables
	int t;			// number of y_i variables
	int origCons_tmp;	// number of original constraints
	double *b_tmp;		// obj coefficients for the x_i variables
	double *c_tmp;		// obj coefficients for the y_i variables
	double **Q_tmp;		// obj coefficients for the x_i_j variables
	double constant;	// obj function additive constant
	double **origMat_tmp;	// original problem matrix
	double *origRhs_tmp;	// original constraints' RHS
	char *origSense_tmp;	// original constraints' sense
	double *xlb_tmp;	// lower bounds on x_i vars
	double *xub_tmp;	// upper bounds on x_i vars
	double *ylb_tmp;	// lower bounds on y_i vars
	double *yub_tmp;	// upper bounds on y_i vars


	int status;

	// read problem file
	status = populateProblem (argv[1],&n, &t, &origCons_tmp, &b_tmp, &c_tmp, 
		&Q_tmp, &constant, &origMat_tmp, &origRhs_tmp, &origSense_tmp, 
		&xlb_tmp, &xub_tmp, &ylb_tmp, &yub_tmp,&si);

	if (status) {
		printf("ERROR! _\n");
		exit(1);
	}

	const double *si_lb = si.getColLower();
	const double *si_ub = si.getColUpper();
	const double *si_objcoeff = si.getObjCoefficients();

#ifdef DEBUG_ORIGINALPROBLEM
printf("[debug original problem]\n");
printf("File               : %s\n",argv[1]);
printf("constant           : %.8f\n",constant);
printf("n (x vars)         : %d\n",n);
printf("t (y vars)         : %d\n",t);
printf("original pb cons   : %d\n",origCons_tmp);
printf("\n");
#ifdef DEBUG_ORIGINALPROBLEM_ARRAYS
printf("b (x coeff)        : ");
for(int i=0;i<n;i++)
printf("%.2f ",b_tmp[i]);
printf("\n");
printf("c (y coeff)        : ");
for(int i=0;i<t;i++)
printf("%.2f ",c_tmp[i]);
printf("\n");
printf("Q (X coeff)        :\n");
for(int i=0;i<n;i++){
for(int j=0;j<n;j++)
printf("%.2f ",Q_tmp[i][j]);
printf("\n");
}
printf("x lb               : ");
for(int i=0;i<n;i++)
printf("%.2f ",xlb_tmp[i]);
printf("\n");
printf("x ub               : ");
for(int i=0;i<n;i++)
printf("%.2f ",xub_tmp[i]);
printf("\n");
printf("y lb               : ");
for(int i=0;i<t;i++)
printf("%.2f ",ylb_tmp[i]);
printf("\n");
printf("y ub               : ");
for(int i=0;i<t;i++)
printf("%.2f ",yub_tmp[i]);
printf("\n");
printf("\n");
printf("sdp+rlt model obj  : ");
for(int i=0;i<si.getNumCols();i++){
printf("[%d]:%.2f ",i,si_objcoeff[i]);
}
printf("\n");
printf("sdp+rlt model lb   : ");
for(int i=0;i<si.getNumCols();i++){
printf("[%d]:%.2f ",i,si_lb[i]);
}
printf("\n");
printf("sdp+rlt model ub   : ");
for(int i=0;i<si.getNumCols();i++){
printf("[%d]:%.2f ",i,si_ub[i]);
}
printf("\n");
printf("\n");
#endif
printf("sdp+rlt model rows : %d  (rlt inequalities = %d)\n",si.getNumRows(),si.getNumRows()-origCons_tmp);
printf("sdp+rlt model cols : %d\n",si.getNumCols());
printf("-----------\n");
#endif


	// determining slack variables
	int slackvars_cnt = 0;

	// determine slack variables for lb ub
	int slacklb_cnt = 0;
	int slackub_cnt = 0;
	// lb_slackvaridx[i] = index of the slack variable associated with the lb constraint for
	//  variable i in the original model si
	int *lb_slackvaridx = new int[si.getNumCols()];
	int *ub_slackvaridx = new int[si.getNumCols()];
	// lb_origidx[i] = index of the variable in the original model si that has the slack
	// variable i associated with the lb constraint for it
	int *lb_origidx = new int[si.getNumCols()];
	int *ub_origidx = new int[si.getNumCols()];
	// rhs of the explicit constraints for the bounds
	double *lb_value = new double[si.getNumCols()];
	double *ub_value = new double[si.getNumCols()];

	for (int i=0;i<si.getNumCols();i++) {
#ifdef ADD_LB
		lb_origidx[slacklb_cnt] = i;
		lb_value[slacklb_cnt] = si_lb[i];
		lb_slackvaridx[slacklb_cnt] = t + slackvars_cnt++;
		slacklb_cnt++;
#endif
#ifdef ADD_UB
		ub_origidx[slackub_cnt] = i;
		ub_value[slackub_cnt] = si_ub[i];
		ub_slackvaridx[slackub_cnt] = t + slackvars_cnt++;
		slackub_cnt++;
#endif
	}

	// determine slack variables for constraints
	int slackcons_cnt = 0;
	const char *si_sense = si.getRowSense();
	int *slackcons_varidx = new int[si.getNumRows()];

	for (int i=0;i<si.getNumRows();i++) {
		slackcons_varidx[i] = -1;
		if (si_sense[i] != 'E') {
#ifndef ADD_RLT
		if (i < origCons_tmp )
#endif
			{
				slackcons_varidx[i] = t +slackvars_cnt++;
				slackcons_cnt++;
			}
		}
	}

	int yvars_first = 0;

	int sdpvars_first = t + slackvars_cnt;


	// variable correspondance
	int *varmap = new int[si.getNumCols()];
	int *varmapsym = new int[si.getNumCols()];
	for(int i=0;i<si.getNumCols();i++) {
		varmap[i] = -1;
		varmapsym[i] = -1;
	}
	for(int i=0;i<n;i++) {
		varmap[i] = sdpvars_first + i+1;
		varmapsym[i] = sdpvars_first + (i+1) * (n+1); 
	}

	for(int i=0;i<n;i++)
		for(int j=i;j<n;j++) {
			varmap[indexQ(i,j,n)] = sdpvars_first + (n+2) + (n+1)*j + i;
			varmapsym[indexQ(i,j,n)] = sdpvars_first + (n+2) + (n+1)*i + j;
		}

	for(int i=0;i<t;i++) {
		varmap[(n*(n+3)/2) + i] = yvars_first + i;
		varmapsym[(n*(n+3)/2) + i] = yvars_first + i;
	}

#ifdef DEBUG_SDPMODEL
printf("[debug sdp model]\n");
printf("slackvars_cnt      : %d\n",slackvars_cnt);
printf("\n");
printf("slackcons_cnt      : %d\n",slackcons_cnt);
printf("slacklb_cnt        : %d\n",slacklb_cnt);
printf("slackub_cnt        : %d\n",slackub_cnt);
printf("\n");
printf("yvars_first        : %d\n",yvars_first);
printf("sdpvars_first      : %d\n",sdpvars_first);
printf("\n");
#ifdef DEBUG_SDPMODEL_ARRAYS
printf("lb_origidx         : ");
for (int i=0;i<slacklb_cnt;i++) {
printf("[%d]:%d  ",i,lb_origidx[i]);
}
printf("\n");
printf("lb_slackvaridx     : ");
for (int i=0;i<slacklb_cnt;i++) {
printf("[%d]:%d  ",i,lb_slackvaridx[i]);
}
printf("\n");
printf("lb_value           : ");
for (int i=0;i<slacklb_cnt;i++) {
printf("[%d]:%.1f  ",i,lb_value[i]);
}
printf("\n");
printf("ub_origidx         : ");
for (int i=0;i<slackub_cnt;i++) {
printf("[%d]:%d  ",i,ub_origidx[i]);
}
printf("\n");
printf("ub_slackvaridx     : ");
for (int i=0;i<slackub_cnt;i++) {
printf("[%d]:%d  ",i,ub_slackvaridx[i]);
}
printf("\n");
printf("ub_value           : ");
for (int i=0;i<slackub_cnt;i++) {
printf("[%d]:%.1f  ",i,ub_value[i]);
}
printf("\n");
printf("varmap             : ");
for (int i=0;i<si.getNumCols();i++)
printf("[%d]>%d ",i,varmap[i]);
printf("\n");
printf("varmapsym          : ");
for (int i=0;i<si.getNumCols();i++)
printf("[%d]>%d ",i,varmapsym[i]);
printf("\n");
#endif
#endif


	CoinPackedMatrix *A = new CoinPackedMatrix();
	CoinPackedVector b(0,(double*)NULL,(int*)NULL);
	CoinPackedVector c(0,(double*)NULL,(int*)NULL);

	int origcons = origCons_tmp; //original constraints (not the RLT)
#ifdef ADD_RLT
	origcons = si.getNumRows(); //original constraints + RLT constraints
#endif
	int A_rows = slacklb_cnt + slackub_cnt + origcons + 1;
	int A_cols = slacklb_cnt + slackub_cnt + slackcons_cnt + t + (n+1)*(n+1);
	
	A->setDimensions(A_rows,A_cols);

#ifdef DEBUG_SDPMODEL
printf("A_rows             : %d\n",A_rows);
printf("A_cols             : %d\n",A_cols);
#endif

	int cons = 0;

	// add lb ub constraints
	for (int i=0;i<slacklb_cnt;i++){
		A->modifyCoefficient(cons,varmap[lb_origidx[i]],
				A->getCoefficient(cons,varmap[lb_origidx[i]]) + 0.5);
		A->modifyCoefficient(cons,varmapsym[lb_origidx[i]],
				A->getCoefficient(cons,varmapsym[lb_origidx[i]]) + 0.5);


		A->modifyCoefficient(cons,lb_slackvaridx[i],-1.0);
		b.insert(cons,lb_value[i]);
		cons++;
	}
	for (int i=0;i<slackub_cnt;i++){
		A->modifyCoefficient(cons,varmap[ub_origidx[i]],
				A->getCoefficient(cons,varmap[ub_origidx[i]]) + 0.5);
		A->modifyCoefficient(cons,varmapsym[ub_origidx[i]],
				A->getCoefficient(cons,varmapsym[ub_origidx[i]]) + 0.5);

		A->modifyCoefficient(cons,ub_slackvaridx[i],1.0);
		b.insert(cons,ub_value[i]);
		cons++;
	}
	
	// add constraints that set first element to 1
	A->modifyCoefficient(cons,sdpvars_first,1.0);
	b.insert(cons,1.0);
	cons++;
	
	// add constraints from si
	const double *rhs = si.getRightHandSide();
	const CoinPackedMatrix *matrixByRow = si.getMatrixByRow();

	for (int i=0;i<si.getNumRows();i++) {
#ifndef ADD_RLT
		if (i>=origCons_tmp)
			break;
#endif
		if (si_sense[i] == 'G') {
			A->modifyCoefficient(cons,slackcons_varidx[i],-1.0);
		}
		if (si_sense[i] == 'L') {
			A->modifyCoefficient(cons,slackcons_varidx[i],1.0);
		}
		const CoinShallowPackedVector row = matrixByRow->getVector(i);
		const int *indices = row.getIndices();
		const double *elements = row.getElements();
		for(int j=0;j<row.getNumElements();j++) {
			// assuming that populate returns a model with the following order of vars
			// x X y

			A->modifyCoefficient(cons,varmap[indices[j]],
				A->getCoefficient(cons,varmap[indices[j]]) + (0.5 * elements[j]));
			A->modifyCoefficient(cons,varmapsym[indices[j]],
				A->getCoefficient(cons,varmapsym[indices[j]]) + (0.5 * elements[j]));
		}
		b.insert(cons,rhs[i]);
		cons++;
	}

	//objective function
	const double *coeff = si.getObjCoefficients();
	double *newCoeff = new double[A_cols];
	for(int i=0;i<A_cols;i++)
		newCoeff[i] = 0.0;
	for (int i=0;i<si.getNumCols();i++) {
		newCoeff[varmap[i]] += -0.5 * coeff[i];
		newCoeff[varmapsym[i]] += -0.5 * coeff[i];
	}
	newCoeff[sdpvars_first] = - constant; // sdpvars_first corresponds to the variable which is always 1
	for(int i=0;i<A_cols;i++)
		if (newCoeff[i] != 0.0)
			c.insert(i,newCoeff[i]);	


#ifdef DEBUG_SDPMODEL_ARRAYS
printf("Sdp model c        : \n");
for(int j=0;j<A_cols;j++) {
printf("[%d]%.1f",j,getCoinPackedVectorElementAt(c,j));
}
printf("\n");
printf("A & RHS            : \n");
for(int i=0;i<A_rows;i++) {
	for(int j=0;j<A_cols;j++) {
#ifdef DEBUG_SDPMODEL_ARRAYS_FULLMATRIX
			printf("%.1f ",A->getCoefficient(i,j));
#else
		if (A->getCoefficient(i,j) != 0.0)
			printf("A[%d][%d]%.1f ",i,j,A->getCoefficient(i,j));
#endif
	}
printf("   RHS=%.1f\n",getCoinPackedVectorElementAt(b,i));
}
#endif
	
	//output mat file
	fprintf(matlabout,"c = [ ");
	
	for (int i=0;i<A_cols;i++){
		fprintf(matlabout,"%f ",getCoinPackedVectorElementAt(c,i));
		if (i!=A_cols-1) {
			fprintf(matlabout,", ");
		}
	}
	fprintf(matlabout,"];\n");

	fprintf(matlabout,"b = [ ");
	for (int i=0;i<A_rows;i++){
		fprintf(matlabout,"%f ",getCoinPackedVectorElementAt(b,i));
		if (i!=A_rows-1) {
			fprintf(matlabout,"; ");
		}
	}
	fprintf(matlabout,"];\n");

	// matrix A input in sparse format
	fprintf(matlabout,"A=[];\n");
	fprintf(matlabout,"A=sparse(A);\n");
	for (int i=0;i<A_rows;i++) {
		for (int j=0;j<A_cols;j++) {
			if (A->getCoefficient(i,j) != 0.0)
				fprintf(matlabout,"A(%d,%d)=%f;\n",i+1,j+1,A->getCoefficient(i,j));
		}
	}

	fprintf(matlabout,"K.f=%d;\n",t);
	fprintf(matlabout,"K.l=%d;\n",slacklb_cnt + slackub_cnt + slackcons_cnt);
	fprintf(matlabout,"K.s=%d;\n",n+1);
	fprintf(matlabout,"%[X,Y,info]=sedumi(A,b,c,K,pars);\n");
	fprintf(matlabout,"[X,Y,Z,info]=csdp(A,b,c,K,pars);\n");
	
	printf("ADD_LB : ");
#ifdef ADD_LB
	printf(" added\n");
#else
	printf(" NOT added\n");
#endif
	printf("ADD_UB : ");
#ifdef ADD_UB
	printf(" added\n");
#else
	printf(" NOT added\n");
#endif
	printf("ADD_RLT: ");
#ifdef ADD_RLT
	printf(" added\n");
#else
	printf(" NOT added\n");
#endif

	delete [] slackcons_varidx;
	delete [] lb_origidx;
	delete [] ub_origidx;
	delete [] lb_slackvaridx;
	delete [] ub_slackvaridx;
	delete [] lb_value;
	delete [] ub_value;

	delete [] varmap;
	delete [] varmapsym;
	delete A;

	return 0;

}
#endif
