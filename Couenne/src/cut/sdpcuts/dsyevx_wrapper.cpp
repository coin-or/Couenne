/* $Id$
 *
 * Name:    dsyevx_rapper.cpp
 * Authors: Andrea Qualizza
 *          Pietro Belotti
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CoinFinite.hpp"
#include "CouenneConfig.h"

//#define DEBUG

extern "C" {

  /* Lapack routine to compute orthonormal eigenvalues/eigenvectors (in Fortran) */

  void F77_FUNC(dsyevx,DSYEVX) (
				  char   *,
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


int dsyevx_interface (int n, double *A, int &m, 
		      double * &w, 
		      double * &z, // output values 
		      double tolerance,
		      double lb_ev, 
		      double ub_ev,
		      int firstidx,
		      int lastidx) {

#ifdef DEBUG

  printf ("matrix:\n---------------------------------\n");
  for (int   i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      printf ("%g ", A [i*n+j]);
    printf ("\n");
  }
  printf ("---------------------------------\n");
#endif

  if (NULL == w) w = new double [n];
  if (NULL == z) z = new double [n*n];

  m = n;

  int lwork = 8*n;

  char jobz  = 'V';  // compute both eigenvalues and eigenvectors
  char range = 'V';  // range for selection is on values of eigenvalues
  char uplo  = 'U';  // upper triangular matrix is given

  int il  = firstidx; // index of the eigenvalue to be returned 1=first
  int iu  = lastidx;  // index of the last eigenvalue to be returuned 1=first
  int lda = n;        // leading dimension of A
  int ldz = n;        // leading dimension of z

  int info; // output status

  int *ifail = new int [n];
  int *iwork = new int [5*n]; 
 
  double abstol = tolerance;	// absolute tolerance
  double vl     = lb_ev;	// minimum eigenvalue wanted
  double vu     = ub_ev;	// maximum

  double *work  = new double [lwork];

  // Equivalent:
  // Ipopt::IpLapackDsyev (true, n, A, lda, w, info);

  F77_FUNC
    (dsyevx,DSYEVX)
    (&jobz, &range, &uplo, &n, 
     A, &lda, 
     &vl, &vu, &il, &iu,
     &abstol, &m, 
     w, z, &ldz, work, &lwork, iwork, ifail, &info);

  if (info) {
    printf (":: dsyevx returned status %d\n", info);
#ifdef CHECK
    for(int i=0; i<m; i++) {
      if(ifail[i] > 0) {
	printf("### WARNING: dsyevx_wrapper(): ifail[%d]: %d   curr_ev[%d]=%.18f\n"
	       , i, ifail [i], ifail [i], w [ifail [i]]);
      }
    }
#endif
  }

  delete [] work;
  delete [] ifail;
  delete [] iwork;

  return m;
}
