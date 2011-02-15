/* $Id$
 *
 * Name:    linquad_cuts.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <linquad_cuts.hpp>
#include <stdio.h>
#include <math.h>

#include <CglCutGenerator.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>

// finds the value xk. (xk,xk^2) is the point on the parabola y=x^2
// closest to the point (xc,yc)
// 
// xk is the zero of the function 2 x^3 - 2 yc x + x - xc = 0


double f_  (double x) {return x*x;}
double fp_ (double x) {return 2*x;}
double fpp_(double x) {return 2;}

double powNewton(double xc, double yc, double (*f)(double),double (*fp)(double),double (*fpp)(double)) {
	// Find a zero to the function
	//
	// F(x) = x - xc + f'(x) (f(x) - yc)
	//
	// where f(x) is either x^k, exp(x), or log(x).
	// The derivative of F(x) is
	//
	// F'(x) = 1 + f''(x) (f(x) - yc) + (f'(x))^2
	//
	// Apply usual update:
	//
	// x(k+1) = x(k) - F(x(k))/F'(x(k))
	double xk	= xc;
	double fk	= f(xk);
	double fpk	= fp(xk);
	double F	= fpk * fk;
	double Fp	= 1+ fpp(xk) * fk + fpk * fpk;
	

	for (int k = NEWTON_MAX_ITER; k--;) {
		xk -= F / Fp;
		
		fk  = f(xk) - yc;
		fpk = fp(xk);
		F   = xk - xc + fpk * fk;
		
		//    printf ("xk = %g; F = %g, fk = %g, fpk = %g\n", xk, F, fk, fpk);
		
		if (fabs(F) < NEWTON_POW_TOLERANCE) 
			break;
		Fp  = 1 + fpp(xk) * fk + fpk * fpk;
	}
	
	return xk;
}


void linQuadCutGen(const double *sol, int n, OsiCuts &cs,Tracer *tracer) {
	// quadratic diagonal cuts
	//	deepest cut: use tangent to quadratic cut passing throgh the 
	//	closest point to sol on the quadratic function
	//	
	//
	int origcuts = cs.sizeRowCuts();
	Timer linquad_timer;
	linquad_timer.start();
	for(int i=0;i<n;i++) {
		double xc,yc;

		//tangent line of region X_ii >= x_i^2 through the closest point 
		//to (sol[x_i],sol[indexQ(i,i,n)]. (Here Xii is our y)
		xc = powNewton(sol[i],sol[indexQ(i,i,n)],&f_,&fp_,&fpp_);
		yc = f_(xc);
		generateTangentDiagonalEntryCut(n,i,cs,xc,yc,sol,true);

		//tangent line of region X_ii >= x_i^2 through the point (sol[x_i],sol[x_i]*sol[x_i])
		xc = sol[i];
		yc = f_(xc);
		generateTangentDiagonalEntryCut(n,i,cs,xc,yc,sol,true);

		//tangent line of region X_ii >= x_i^2 through the point
		//(sqrt(sol[indexQ(i,i,n)]),sol[x_i]*sol[x_i])
		if(sol[indexQ(i,i,n)] > 0) {
			xc = sqrt(sol[indexQ(i,i,n)]);
			yc = sol[indexQ(i,i,n)];
			generateTangentDiagonalEntryCut(n,i,cs,xc,yc,sol,true);
		}
	}
	tracer->setLinquadTime(linquad_timer.time());
	tracer->setLinquadTotalCuts(cs.sizeRowCuts() - origcuts);
}



void linQuadCutGenOriginalBounds(const double *xlb, const double *xub, int n, OsiCuts &cs,Tracer *tracer) {
#define	LINQUAD_BOUNDS_CUTS_PARTS	3
	int origcuts = cs.sizeRowCuts();
	Timer linquad_timer;
	linquad_timer.start();
	for (int i=0;i<n;i++) {
		double xc,yc;
		double interval = xub[i] - xlb[i];
		double step = interval / LINQUAD_BOUNDS_CUTS_PARTS;
		int cnt=0;
		xc = xub[i];
		while (cnt<=LINQUAD_BOUNDS_CUTS_PARTS) {
			yc = f_(xc);
			generateTangentDiagonalEntryCut(n,i,cs,xc,yc,NULL,false);	
			cnt++;
			xc += step;
		}	
	}
	tracer->setLinquadTime(linquad_timer.time());
	tracer->setLinquadTotalCuts(cs.sizeRowCuts() - origcuts);
}


void generateTangentDiagonalEntryCut(int n,int i,OsiCuts &cs,double xc,double yc,const double* sol, bool ifViolated) {
// (xc,yc) is a point that satisfies X_ii = x_i^2
// now we compute the tangent line of the region at this point which will be
// (X_ii - yc) = 2 xc (x_i - xc) -->  - 2*xc*x_i + X_ii = yc - 2*xc^2
// the cut is simply : - 2 x_i + X_ii >= yc - 2*xc^2
// sol can be NULL if ifViolated=false
		double *coeff;
		int *ind;
		double rhs;
		OsiRowCut *cut;

		cut   = new OsiRowCut;
		coeff = new double [2];
		ind   = new int    [2];
		coeff[0]	= -2.0 *xc;
		coeff[1]	= 1.0 ;
		ind[0]		= i ;
		ind[1]		= indexQ(i,i,n);
		rhs	= yc - 2 * xc * xc;
		cut -> setRow (2, ind, coeff);
		cut -> setLb (rhs);

		if ( !( ifViolated ) || ( cut->violated(sol) <= 0) ) {
			cs.insert (cut);
		}

		delete cut;
		delete [] coeff;
		delete [] ind;
}

