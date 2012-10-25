/* $Id$
 *
 * Name:    dsyevx_wrapper.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef DSYEVX_HPP
#define DSYEVX_HPP

// wrapper for Lapack's Fortran routine to compute all eigenvalues/vectors
void dsyevx_full_wrapper (int n, double *A, int &m, double * &w, double * &z);

// wrapper for Lapack's Fortran routine to compute only the positive eigenvalue/vector
void dsyevx_wrapper_only_positive (int, double *, int &, double * &, double * &);

// wrapper for Lapack's Fortran routine to compute only the negative eigenvalue/vector
void dsyevx_wrapper_only_negative (int, double *, int &, double * &, double * &);

// wrapper for Lapack's Fortran routine to compute only the most negative eigenvalue/vector
void dsyevx_wrapper_only_most_neg (int, double *, int &, double * &, double * &);

// wrapper for Lapack's Fortran routine to compute the first p eigenvalues/vectors
void dsyevx_wrapper_first_p (int, double *, int &, double * &, double * &, int);


#endif
