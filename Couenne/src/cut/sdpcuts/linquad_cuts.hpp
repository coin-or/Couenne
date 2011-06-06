/* $Id$
 *
 * Name:    linquad_cuts.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef LINQUAD_CUTS_HPP
#define LINQUAD_CUTS_HPP

#include <CglCutGenerator.hpp>
#include <tracer.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

#define NEWTON_MAX_ITER 10
#define NEWTON_POW_TOLERANCE 1e-12


void linQuadCutGen(const double *sol, OsiCuts &cs);

double f_  (double x);
double fp_ (double x);
double fpp_(double x);
double powNewton(double xc, double yc, double (*f)(double),double (*fp)(double),double (*fpp)(double));
void linQuadCutGen(const double *sol, int n, OsiCuts &cs, Tracer *tracer);
void linQuadCutGenOriginalBounds(const double *xlb, const double *xub, int n, OsiCuts &cs, Tracer *tracer);

void generateTangentDiagonalEntryCut(int n,int i,OsiCuts &cs,double xc,double yc,const double* sol, bool ifViolated);



#endif

