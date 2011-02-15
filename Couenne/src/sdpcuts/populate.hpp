/* $Id$
 *
 * Name:    populate.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef POPULATE_HPP
#define POPULATE_HPP

#define COUENNE_EPS 1e-7 //used in createCut(...)

int createCut (OsiCuts &, 
	       double,
	       int,     // sign: -1 is <=, +1 is >=, 0 is =
	       int=-1, double=0,
	       int=-1, double=0,
	       int=-1, double=0,
	       int=-1, double=0, bool=false);

int populateProblem (const char *filename,int *nptr, int *tptr, int *consptr, double **bptr, 
		     double **cptr, double ***Qptr, double *constantptr, double ***origmatptr, 
		     double **origrhsptr, char **origsenseptr, double **xlbptr,double **xubptr,
		     double **ylbptr, double **yubptr,OsiSolverInterface *si);

#endif
