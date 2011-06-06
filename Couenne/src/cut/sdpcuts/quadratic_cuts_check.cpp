/* $Id$
 *
 * Name:    quadratic_cuts_check.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <quadratic_cuts_check.hpp>
#include <stdio.h>
#include <dsyevx_wrapper.hpp>
#include <tracer.hpp>

QuadraticCuts::QuadraticCuts(int n, const double *initial_sol, Tracer *tracer) {
	n_ = n;
	L = new double[n*n];
	tracer_ = tracer;
	
	sol = new double[(n*(n+3))/2];
	previous_sol = new double[(n*(n+3))/2];
	for (int i=0;i<(n*(n+3))/2;i++) {
		sol[i] = 0.0;
		previous_sol[i] = 0.0;
	}
	eigenvectors = new double*[n];
	for (int i=0;i<n;i++)
		eigenvectors[i] = new double[n];

#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	Xtilde = new double[(n+1)*(n+1)];
	eigenvectors_Xtilde = new double*[n+1];
	for (int i=0;i<n+1;i++)
		eigenvectors_Xtilde[i] = new double[n+1];
#endif
	updateSolution(initial_sol);

	checkQuadraticDiagonalCutsOnCurrentSolution();

	computeEigenvectorsFromCurrentSolution();
}

QuadraticCuts::~QuadraticCuts() {
	delete [] L;
	delete [] sol;
	delete [] previous_sol;
	for (int i=0;i<n_;i++)
		delete [] eigenvectors[i];
	delete [] eigenvectors;

#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	delete [] Xtilde;
	for (int i=0;i<n_+1;i++)
		delete [] eigenvectors_Xtilde[i];
	delete [] eigenvectors_Xtilde;
#endif
}

void QuadraticCuts::refresh(const double *current_sol) {
	updateSolution(current_sol);

	checkQuadraticDiagonalCutsOnCurrentSolution();
	checkPreviousQuadraticEVCutsOnCurrentSolution();

	computeEigenvectorsFromCurrentSolution();
}

void QuadraticCuts::updateSolution(const double *current_sol){

	// L = X - xxT
	double *X = new double[n_*n_];
	double *xxT = new double[n_*n_];

	for (int i=0;i<(n_*(n_+3))/2;i++)
		previous_sol[i] = sol[i];
	
	for (int i=0;i<(n_*(n_+3))/2;i++)
		sol[i] = current_sol[i];
	
	for (int i=0;i<n_;i++) {
		for (int j=i;j<n_;j++) {
			X   [i*n_ + j] = current_sol[indexQ(i,j,n_)];
			X   [j*n_ + i] = current_sol[indexQ(i,j,n_)];
			xxT [i*n_ + j] = current_sol[i]*current_sol[j];
			xxT [j*n_ + i] = current_sol[i]*current_sol[j];
		}
	}
	for (int i=0;i<n_;i++) {
		for (int j=0;j<n_;j++) {
			L[i*n_ + j] = X[i*n_ + j] - xxT[i*n_ + j];
		}
	}

#ifdef QUADRATIC_CUTS_DEBUG
//	printf("current solution = ");
//	for (int i=0;i<(n_*(n_+3))/2;i++)
//		printf("%4f ",current_sol[i]);
//	printf("\n");
	printf("\n");
	printf("x     :\t");
	for (int i=0;i<n_;i++)
		printf("%8.5f\t",current_sol[i]);
	printf("\n");
	printf("X     :\t");
	for (int i=0;i<n_;i++) {
		for (int j=0;j<n_;j++)
			printf("%8.5f\t",X   [i*n_ + j]);
		printf("\n\t");
	}
	printf("\n");
	printf("xxT   :\t");
	for (int i=0;i<n_;i++) {
		for (int j=0;j<n_;j++)
			printf("%8.5f\t",xxT   [i*n_ + j]);
		printf("\n\t");
	}
	printf("\n");
	printf("X-xxT :\t");
	for (int i=0;i<n_;i++) {
		for (int j=0;j<n_;j++)
			printf("%8.5f\t",L   [i*n_ + j]);
		printf("\n\t");
	}
	printf("\n");
#endif

	delete [] X;
	delete [] xxT;

#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	int np=n_+1;
	for (int i=0;i<np*np;i++)
		Xtilde[i] = 0;
	Xtilde[0] = 1;
	for (int i=0;i<n_;i++) {

		Xtilde[i+1] = current_sol[i];
		Xtilde[(i+1)*np] = current_sol[i];
		for (int j=i;j<n_;j++) {
			Xtilde[(i+1)*np+(j+1)] = current_sol[indexQ(i,j,n_)];
			Xtilde[(j+1)*np+(i+1)] = current_sol[indexQ(i,j,n_)];
		}
	}
#ifdef QUADRATIC_CUTS_DEBUG
	printf("Xtilde=\n");
	for (int i=0;i<np;i++) {
		for (int j=0;j<np;j++)
			printf("%.5f ",Xtilde[i*np + j]);
		printf("\n");
	}
#endif // QUADRATIC_CUTS_DEBUG
#endif // RECOMPUTE_XTILDE_EV_FROM_SCRATCH
}



void QuadraticCuts::computeEigenvectorsFromCurrentSolution() {
	double *z = new double[n_ * n_];
	double *w = new double[n_];
	int m;
		
	dsyevx_wrapper_only_negative (n_, L, m, w, z,tracer_);
	
	card_ev = 0;
	for (int i=0;i<m;i++) {
		if (w[i] < 0) {
			for (int j=0;j<n_;j++)
				eigenvectors[card_ev][j] = z[i*n_ + j];
			card_ev++;
		} else
			break;
	}
#ifdef QUADRATIC_CUTS_DEBUG
	printf("eigenvectors of X - xx^T:\n");
	for (int i=0;i<card_ev;i++) {
	printf("[%d] ev=%.5f : ",i,w[i]);
	for (int j=0;j<n_;j++)
		printf("%.5f ",eigenvectors[i][j]);
	printf("\n");
	}
	printf("\n");
#endif
	delete [] z;
	delete [] w;

#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	int np = n_+1;
	z = new double[np * np];
	w = new double[np];

	dsyevx_full_wrapper (np, Xtilde, m, w, z, tracer_);

	card_ev_Xtilde = 0;
	for (int i=0;i<m;i++) {
		if (w[i] < 0) {
			for (int j=0;j<np;j++)
				eigenvectors_Xtilde[card_ev_Xtilde][j] = z[i*np + j];
			card_ev_Xtilde++;
		} else
			break;
	}

#ifdef QUADRATIC_CUTS_DEBUG
	printf("eigenvectorsXtilde:\n");
	for (int i=0;i<card_ev_Xtilde;i++) {
	printf("[%d] ev=%.5f : ",i,w[i]);
	for (int j=0;j<np;j++)
		printf("%.5f ",eigenvectors_Xtilde[i][j]);
	printf("\n");
	}
#endif
	delete [] z;
	delete [] w;
#endif
}

void QuadraticCuts::checkQuadraticDiagonalCutsOnCurrentSolution() {
	for (int i=0;i<n_;i++) {
		double value = sol[indexQ(i,i,n_)] - sol[i]*sol[i];
		if (value < - QUADRATIC_CUTS_CHECK_TOLERANCE)
			printf("quadratic cut X_%d,%d - (x_%d)^2 >= 0 violated (lhs=%.5f)\n",i,i,i,value);
	}
}

void QuadraticCuts::checkPreviousQuadraticEVCutsOnCurrentSolution() {
	for(int i=0;i<card_ev;i++) {
		double vTx = 0.0;
		for (int j=0;j<n_;j++) {
			vTx += sol[j] * eigenvectors[i][j];
		}
		double XvvT = 0.0;
		for (int j=0;j<n_;j++) {
			for (int k=j+1;k<n_;k++)
				XvvT += 2*eigenvectors[i][j]*eigenvectors[i][k]*sol[indexQ(j,k,n_)];
			XvvT += eigenvectors[i][j]*eigenvectors[i][j]*sol[indexQ(j,j,n_)];
		}
		double lhs = vTx*vTx - XvvT;

#ifdef QUADRATIC_CUTS_DEBUG
		printf("ev[%d]  v:\t",i);
		for (int k=0;k<n_;k++) {
			printf("%8.5f\t",eigenvectors[i][k]);
		}
		printf("\n");
		printf("ev[%d]vvT:\t",i);
		for (int k=0;k<n_;k++) {
			for (int j=0;j<n_;j++)
				printf("%8.5f\t",eigenvectors[i][k]*eigenvectors[i][j]);
 			printf("\n\t\t");
		}
		printf("\n");
		printf("curr_ev_idx=%d vTx=%.5f vTx^2=%.5f XvvT=%.5f\n",i,vTx,vTx*vTx,XvvT);
		printf("\n");
#endif

		if (lhs > QUADRATIC_CUTS_CHECK_TOLERANCE) {
			printf("quadratic cut from %d-th prev ev (v^Tx)^2 - X.vv^T <=0)violated (lhs=%.5f)\n",i,lhs);
		}
	}


#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	int np = n_+1;
	double **wwT = new double*[np];
	double **myXtilde = new double*[np];
	for (int i=0;i<np;i++) {
		wwT[i] = new double[np];
		myXtilde[i] = new double[np];
	}

	myXtilde[0][0] = 1.0;
	for (int i=0;i<n_;i++) {
		myXtilde[0][i+1] = sol[i];
		myXtilde[i+1][0] = sol[i];
		for (int j=i;j<n_;j++) {
			myXtilde[i+1][j+1] = sol[indexQ(i,j,n_)];
			myXtilde[j+1][i+1] = sol[indexQ(i,j,n_)];
		}
	}

#ifdef QUADRATIC_CUTS_DEBUG
	printf("myXtilde:\n");
	for(int i=0;i<np;i++) {
		for (int j=0;j<np;j++)
			printf("%.2f ",myXtilde[i][j]);
		printf("\n");
	}
#endif

	for(int i=0;i<card_ev_Xtilde;i++) {
		for (int k=0;k<np;k++) {
			for (int j=0;j<np;j++) {
				wwT[k][j] = eigenvectors_Xtilde[i][k] * eigenvectors_Xtilde[i][j];
			}
		}
		double value = 0.0;
		for (int k=0;k<np;k++) {
			for (int j=0;j<np;j++) {
				value+=myXtilde[k][j]*wwT[k][j];
			}
		}
		if (value < -1e-8)
			printf("Original Xtilde ev[%d] cut violated (this should not happen!) = %.8f\n",i,value);
	}

	for (int i=0;i<np;i++) {
		delete [] wwT[i];
		delete [] myXtilde[i];
	}
	delete [] wwT;
#endif // RECOMPUTE_XTILDE_EV_FROM_SCRATCH
}




