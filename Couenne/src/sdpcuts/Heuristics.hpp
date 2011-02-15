/* $Id$
 *
 * Name:    Heuristics.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef HEURISTICS_HPP
#define HEURISTICS_HPP

#include "CoinPackedVector.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiXxxSolverInterface.hpp"
#include "tracer.hpp"
#include "misc_util.hpp"


class Heuristics {
private:
	const int n_;		// number of variables x_i (involved in quadratic terms)
	const int t_;		// number of variables y_i (NOT involved in quadratic terms)
	int N_;			// number of variables x_i + X_ij  (= n*(n+3)/2)
	const int cons_;	// number of constraints in the original problem
	const double objConst_;	// constant to be added to the objective function value
	const double *b_;	// objective function coefficients of the variables x_i
	const double *c_;	// objective function coefficients of the variables y_i
	const double **Q_;	// objective function coefficients of the variables X_ij
	const double **origMat_;// original problem constraints (linearized)
	const double *origRhs_;	// original problem RHS
	const char *origSense_;	// original problem sense
	const double *xlb_;	// original problem lower bounds for x_i
	const double *xub_;	// original problem upper bounds for x_i
	const double *ylb_;	// original problem lower bounds for y_i
	const double *yub_;	// original problem upper bounds for y_i
	const OsiSolverInterface *si_;// pointer to the main problem

	double currObj_;	// current quadratic objective 
	double bestObj_;	// quadratic objective of bestSol_
	double *bestSol_;	// best solution of original problem; dim. N_ + t_

	bool *heurLbRowAdded_;	// used in heurSi to determine which original 
				// constraints are considered (those with nz at some y_i entries)
				// the other constraints can be discarded

	double *xxTSol_;	// (x,xxT,y) heuristic solution

	// min norm heuristic variables
	double *MNSol_;		// approximate min norm heuristic solution
	double *temp_row_;
	OsiXxxSolverInterface MNLPSi_; 

	OsiXxxSolverInterface heurLPimproveSi_; // heuristic solution improvement LP solver interface
	Tracer *tracer_;	//global tracer

private:
	// Evaluates the current heuristic solution and updates currObj and bestSol if necessary
	int update(double*, double);

	// Tries to improve the current solution by fixing the values for x_i and X_ij 
	// and optimizing over the y_i variables only.
	// Would not produce any improvement if the original problem has no variable y_i (i.e. t=0)
	int heurLP_improveSolution(double*);

	// uses the functions heurLP_improveSolution and feasibility check, and updates current & best solutions
	int processSol(double*, bool, double*, double*);

	double* xxTHeur();
	double* MNHeur();

public:
	// constructor
	Heuristics(const int,const int,const int,const double,const double*,const double*,const double**,const double**,const double*,const char*,const double*,const double*,const double*,const double*,const OsiSolverInterface *si,Tracer *tracer);
	
	// destructor
  	~Heuristics();

	// getters
	double bestObj() {return bestObj_;}
	double* bestSol() {return bestSol_;}
	double currObj() {return currObj_;}

	// run heuristics on the current outer approximation solution
	int run();
};


#endif

