/* $Id$
 *
 * Name:    quadratic_cuts_check.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef QUADRATIC_CUTS_CHECK_HPP
#define QUADRATIC_CUTS_CHECK_HPP

#include <tracer.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))



#define aaRECOMPUTE_XTILDE_EV_FROM_SCRATCH
#define QUADRATIC_CUTS_DEBUG
#define QUADRATIC_CUTS_CHECK_TOLERANCE 1e-8

class QuadraticCuts{
public:
	QuadraticCuts(int n, const double *initial_sol, Tracer *tracer);
	~QuadraticCuts();
	void refresh(const double *current_sol);

private:
	int n_;
	double *L; //L = X - xxT
	double *sol;
	double *previous_sol;
	double **eigenvectors;
	int card_ev;
#ifdef RECOMPUTE_XTILDE_EV_FROM_SCRATCH
	double *Xtilde;
	double **eigenvectors_Xtilde;
	int card_ev_Xtilde;
#endif
	Tracer *tracer_;

	void updateSolution(const double *current_sol);
	void computeEigenvectorsFromCurrentSolution();
	void checkQuadraticDiagonalCutsOnCurrentSolution();
	void checkPreviousQuadraticEVCutsOnCurrentSolution();
};

#endif

