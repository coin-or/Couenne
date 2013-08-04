/* $Id$
 *
 * Name:    dsyevx_wrapper.hpp
 * Authors: Andrea Qualizza
 *          Pietro Belotti
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef DSYEVX_HPP
#define DSYEVX_HPP

int dsyevx_interface (int n, double *A, int &m, 
		       double * &w, double * &z, // output values 
		       double tolerance,
		       double lb_ev, 
		       double ub_ev,
		       int firstidx,
		       int lastidx);

#endif
