/* $Id$
 *
 * Name:    dsyevx_rapper.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <misc_util.hpp>

#include "CoinHelperFunctions.hpp"

#include <tracer.hpp>

#define ABS_TOL_EIG 1e-15



void _dsyevx_value_range_wrapper (int n, double *A, int &m, double * &w, double * &z,double tolerance,double lb_ev, double ub_ev);

void _dsyevx_index_range_wrapper (int n, double *A, int &m, double * &w, double * &z, double tolerance, int firstidx,int lastidx);


#if 0
extern "C" {

/* Lapack routine to compute orthonormal eigenvalues/eigenvectors (in Fortran) */
void dsyevx_ (char   *,
	      char   *,
	      char   *,
	      int    *,
	      double *,
	      int    *,
	      double *,
	      double *,
	      int    *,
	      int    *,
	      double *,
	      int    *,
	      double *,
	      double *,
	      int    *,
	      double *,
	      int    *,
	      int    *,
	      int    *,
	      int    *);
}

#else

/* Lapack routine to compute orthonormal eigenvalues/eigenvectors (in Fortran) */
void dsyevx_ (char   *,
	      char   *,
	      char   *,
	      int    *,
	      double *,
	      int    *,
	      double *,
	      double *,
	      int    *,
	      int    *,
	      double *,
	      int    *,
	      double *,
	      double *,
	      int    *,
	      double *,
	      int    *,
	      int    *,
	      int    *,
	      int    *) {}
#endif



void dsyevx_full_wrapper (int n, double *A, int &m, double * &w, double * &z, Tracer *tracer) {
	tracer->incrementMainTotalEigendecompositions();
	_dsyevx_value_range_wrapper (n,A,m,w,z,ABS_TOL_EIG,-COIN_DBL_MAX,COIN_DBL_MAX);
}

void dsyevx_wrapper_only_positive (int n, double *A, int &m, double * &w, double * &z, Tracer *tracer) {
	tracer->incrementMainTotalEigendecompositions();
	_dsyevx_value_range_wrapper (n,A,m,w,z,ABS_TOL_EIG,0.0,COIN_DBL_MAX);
}

void dsyevx_wrapper_only_negative (int n, double *A, int &m, double * &w, double * &z, Tracer *tracer) {
	tracer->incrementMainTotalEigendecompositions();
	_dsyevx_value_range_wrapper (n,A,m,w,z,ABS_TOL_EIG,-COIN_DBL_MAX,ABS_TOL_EIG);
}

void dsyevx_wrapper_only_most_neg (int n, double *A, int &m, double * &w, double * &z, Tracer *tracer) {
	tracer->incrementMainTotalEigendecompositions();
	_dsyevx_index_range_wrapper (n,A,m,w,z,ABS_TOL_EIG,1,1);
}

void dsyevx_wrapper_first_p (int n, double *A, int &m, double * &w, double * &z, int p, Tracer *tracer) {
	tracer->incrementMainTotalEigendecompositions();
	_dsyevx_index_range_wrapper (n,A,m,w,z,ABS_TOL_EIG,1,p);
}

//########################################################################################


void _dsyevx_value_range_wrapper (int n, double *A, int &m, double * &w, double * &z,double tolerance,double lb_ev, double ub_ev) {
// the lapack call destroys the original matrix A
	if(w == NULL)
		w = new double[n];
	if(z == NULL)
		z = new double[n*n];

	m = n;

	int lwork = 8*n;

	char jobz  = 'V';  // compute both eigenvalues and eigenvectors
	char range = 'V';  // range for selection is on values of eigenvalues
        char uplo  = 'U';  // upper triangular matrix is given

	int il     = 1;   // first  eigenvalue to be returned (not used)
	int iu     = n;   // second                           (not used)
	int info;         // output status
	int lda    = n;   // leading dimension of A
	int ldz    = n;   // leading dimension of z

	int *ifail = new int[n];
	int *iwork = new int[5*n]; 

 
	double abstol = tolerance;	// absolute tolerance
	double vl     = lb_ev;		// minimum eigenvalue wanted
	double vu     = ub_ev;		// maximum

	double *work  = new double[lwork];

	dsyevx_ (&jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

	if (info) {
		printf (":: dsyevx returned status %d\n", info);
#ifdef CHECK
		for(int i=0; i<m; i++) {
			if(ifail[i] > 0) {
				printf("### WARNING: dsyevx_wrapper(): ifail[%d]: %d   curr_ev[%d]=%.18f\n"
					,i, ifail[i],ifail[i],w[ifail[i]]);
			}
		}
#endif
	}

	if ((info == 0) && (m > 0)) { // there is at least one eigenvector
#ifdef TRACE_ALL
		cpp_printvecDBL("eigenvalues", w, m);
#endif
	}
	
	delete [] work;
	delete [] ifail;
	delete [] iwork;
}


//########################################################################################


void _dsyevx_index_range_wrapper (int n, double *A, int &m, double * &w, double * &z, double tolerance, int firstidx,int lastidx) {
// the lapack call destroys the original matrix A
	if(w == NULL)
		w = new double[n];
	if(z == NULL)
		z = new double[n*n];

	m = n;

	int lwork = 8*n;

	char jobz  = 'V';  // compute both eigenvalues and eigenvectors
	char range = 'I';  // range selection
        char uplo  = 'U';  // upper triangular matrix is given

	int il     = firstidx;   // index of the eigenvalue to be returned 1=first
	int iu     = lastidx;   // index of the last eigenvalue to be returuned 1=first
	int info;         // output status
	int lda    = n;   // leading dimension of A
	int ldz    = n;   // leading dimension of z

	int *ifail = new int[n];
	int *iwork = new int[5*n]; 

 
	double abstol = ABS_TOL_EIG;	// absolute tolerance
	double vl     = 0.0;		// not used
	double vu     = 0.0;		// not used

	double *work  = new double[lwork];

	dsyevx_ (&jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

	if (info) {
		printf (":: dsyevx returned status %d\n", info);
#ifdef CHECK
		for(int i=0; i<m; i++) {
			if(ifail[i] > 0) {
				printf("### WARNING: dsyevx_wrapper(): ifail[%d]: %d   curr_ev[%d]=%.18f\n"
					,i, ifail[i],ifail[i],w[ifail[i]]);
			}
		}
#endif
	}

	if ((info == 0) && (m > 0)) { // there is at least one eigenvector
#ifdef TRACE_ALL
		cpp_printvecDBL("eigenvalues", w, m);
#endif
	}
	
	delete [] work;
	delete [] ifail;
	delete [] iwork;
}



