/* $Id$
 *
 * Name:    sdpcuts.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef SDPCUTS_HPP
#define SDPCUTS_HPP

#include <OsiXxxSolverInterface.hpp>

void print_ifdefs(FILE *);

void solver_status(OsiSolverInterface*);

FILE *open_f_res();

FILE *open_short_f_res();

void print_current_sol(int , double , int , int , double ,double , double  );

void print_file_current_sol(FILE* , char* , int , double , int , int , double ,double , double );

void print_file_short_sol(FILE* ,char* , int , double , int , int , double ,double );


#define LP_TOLERANCE 0.0000001
// values smaller than 0.0000001 st_e05 gives a outerapproximation infeasible solution with cplex

#define FEAS_CHECK_NO_VIOLATION 0
#define FEAS_CHECK_BOUNDS_VIOLATION 1
#define FEAS_CHECK_CONSTRAINT_VIOLATION 2
#define FEAS_CHECK_CONSTRAINT_VIOLATION_NO_RECOVER 3

int feasibility_check(const int , const int , const int , const double *, const double ** , const double *, const char *, const double *,const double *,const double *,const double *);

double evaluateSolution(const int , const int ,const double *, const double *, const double * , const double **, const double );


#endif
