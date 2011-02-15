/* $Id$
 *
 * Name:    tracer.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef TRACER_HPP
#define TRACER_HPP
#include <misc_util.hpp>

#define VERSION_TRACER 1
#define REPORT_DELIMITER "|"
#define REPORT_NO_ENTRY "!"
#define TRACER_GLOBAL_REPORT_FILE "globalreport.txt"
#define TRACER_DETAILED_REPORT_FILE "detailedreport.txt"

#define TRACER_HEURISTICS_COMPARISON_TOLERANCE 1e-8

#define TRACER_INVALID_ENTRY -999999999

#define EXIT_ON_ITER 1000

#define TRACE_ITER_SIZE EXIT_ON_ITER+2


typedef struct bs {
	double time;
	double bound;
	int iter1;
	int iter2;
} bound_struct;

class Tracer {
public:
	Tracer(char *instance_param, int x_vars_param, int X_vars_param, int y_vars_param);
	~Tracer();

	// control functions
	void newIter();
	int iterations() const;
	void detailedReport() const;
	void globalReport() const;

	// functions for main entries
	void setMainBound(double value);
	void setMainIterationTime(double value);
	void setMainLPTime(double value);
	void setMainActiveCuts(int value);
	void setMainAddedCuts(int value);
	void setMainTotalCuts(int value);
	void setMainTotalEigendecompositions(int value);
	void incrementMainTotalEigendecompositions();
	void setMainDeletedCuts(int value);

	// functions for sdp cuts
	void setSDPNumNegativeEV(int value);
	void setSDPMostNegativeEV(double value);
	void setSDPCutsTime(double value);
	void setSDPCutsTotalCuts(int value);


	// functions for sparsify entries
	void setSparsifyTime(double value);
	void setSparsifyTotalCuts(int value);
	void setSparsifyDuplicatedCuts(int value);
	void setSparsifyWiseDecompositions(int value);
	void addSparsifyNz(int value);
	void addSparsifySingleColumnSparsity(int value);
	void addSparsifyColumnPairSparsity(int value);
	void addSparsifyTop20PercCutsViolation(double value);


	// functions for orthocut
	void setOrthocutTime(double value);
	void setOrthocutTotalCuts(int value);


	// functions for linquad cuts
	void setLinquadTime(double value);
	void setLinquadTotalCuts(int value);


	// functions for disjunctive cuts
	void setDisjunctiveCutsTime(double value);
	void setDisjunctiveCutsTotalCuts(int value);


	// functions for heuristics entries
	void setHeuristicsCurrentSolution(double value);
	void setHeuristicsBestSolution(double value);
	void setHeuristicsTime(double value);
	void setHeuristicsxxTSolution(double value);
	void setHeuristicsxxTSolutionLPHeuristicImprovement(double value);
	void setHeuristicsxxTTime(double value);
	void setHeuristicsMNLPSolution(double value);
	void setHeuristicsMNLPSolutionLPHeuristicImprovement(double value);
	void setHeuristicsMNLPTime(double value);
	void setHeuristicsGWSolution(double value);
	void setHeuristicsGWSolutionLPHeuristicImprovement(double value);
	void setHeuristicsGWTime(double value);

	int iter;
	int x_vars;
	int X_vars;
	int y_vars;
	int tot_vars;

	char instance[200];
	
	double	main__iteration_bound[TRACE_ITER_SIZE];
	double	main__time[TRACE_ITER_SIZE];
	double	main__iteration_time[TRACE_ITER_SIZE];
	double	main__lp_time[TRACE_ITER_SIZE];
	Stat	main__lp_time_global_stat;
	int	main__total_cuts[TRACE_ITER_SIZE];
	Stat	main__total_cuts_global_stat;
	int	main__active_cuts[TRACE_ITER_SIZE];
	Stat	main__active_cuts_global_stat;
	int	main__added_cuts[TRACE_ITER_SIZE];
	Stat	main__added_cuts_global_stat;
	int	main__deleted_cuts[TRACE_ITER_SIZE];
	Stat	main__deleted_cuts_global_stat;
	int	main__total_eigendecompositions[TRACE_ITER_SIZE];
	Stat	main__total_eigendecompositions_global_stat; //use only for sum() !!!

	int	sdpcuts__num_negative_ev[TRACE_ITER_SIZE];
	double	sdpcuts__most_negative_ev[TRACE_ITER_SIZE];
	double	sdpcuts__time[TRACE_ITER_SIZE];
	Stat	sdpcuts__time_global_stat;
	int	sdpcuts__total_cuts[TRACE_ITER_SIZE];
	Stat	sdpcuts__total_cuts_global_stat;

	double	sparsify__time[TRACE_ITER_SIZE];
	Stat	sparsify__time_global_stat;
	int	sparsify__total_cuts[TRACE_ITER_SIZE];
	Stat	sparsify__total_cuts_global_stat;
	int	sparsify__duplicated_cuts[TRACE_ITER_SIZE];
	Stat	sparsify__duplicated_cuts_global_stat;
	int	sparsify__wise_decompositions[TRACE_ITER_SIZE];
	Stat	sparsify__wise_decompositions_global_stat;
	Stat	sparsify__nz_global_stat;
	Stat	sparsify__nz_iter_stat[TRACE_ITER_SIZE];
	bool	sparsify__nz_populated;
	Stat	sparsify__single_column_sparsity_global_stat;
	Stat	sparsify__single_column_sparsity_iter_stat[TRACE_ITER_SIZE];
	bool	sparsify__single_column_sparsity_populated;
	Stat	sparsify__column_pairs_sparsity_global_stat;
	Stat	sparsify__column_pairs_sparsity_iter_stat[TRACE_ITER_SIZE];
	bool	sparsify__column_pairs_sparsity_populated;
	Stat	sparsify__top20perc_cuts_violation_global_stat;
	Stat	sparsify__top20perc_cuts_violation_iter_stat[TRACE_ITER_SIZE];
	bool	sparsify__top20perc_cuts_violation_populated;

	int	orthocut__total_cuts[TRACE_ITER_SIZE];
	Stat	orthocut__total_cuts_global_stat;
	double	orthocut__time[TRACE_ITER_SIZE];
	Stat	orthocut__time_global_stat;

	int	linquadcuts__total_cuts[TRACE_ITER_SIZE];
	Stat	linquadcuts__total_cuts_global_stat;
	double	linquadcuts__time[TRACE_ITER_SIZE];
	Stat	linquadcuts__time_global_stat;

	int	disjunctivecuts__total_cuts[TRACE_ITER_SIZE];
	Stat	disjunctivecuts__total_cuts_global_stat;
	double	disjunctivecuts__time[TRACE_ITER_SIZE];
	Stat	disjunctivecuts__time_global_stat;

	double	heuristics__current_solution[TRACE_ITER_SIZE];
	double	heuristics__best_solution[TRACE_ITER_SIZE];
	double	heuristics__time[TRACE_ITER_SIZE];
	Stat	heuristics__time_global_stat;
	double	heuristics__xxT_solution[TRACE_ITER_SIZE];
	double	heuristics__xxT_solution_lp_heuristic_improvement[TRACE_ITER_SIZE];
	Stat	heuristics__xxT_solution_lp_heuristic_improvement_global_stat;
	double	heuristics__xxT_time[TRACE_ITER_SIZE];
	Stat	heuristics__xxT_time_global_stat;
	double	heuristics__MNLP_solution[TRACE_ITER_SIZE];
	double	heuristics__MNLP_solution_lp_heuristic_improvement[TRACE_ITER_SIZE];
	Stat	heuristics__MNLP_solution_lp_heuristic_improvement_global_stat;
	double	heuristics__MNLP_time[TRACE_ITER_SIZE];
	Stat	heuristics__MNLP_time_global_stat;
	double	heuristics__GW_solution[TRACE_ITER_SIZE];
	double	heuristics__GW_solution_lp_heuristic_improvement[TRACE_ITER_SIZE];
	Stat	heuristics__GW_solution_lp_heuristic_improvement_global_stat;
	double	heuristics__GW_time[TRACE_ITER_SIZE];
	Stat	heuristics__GW_time_global_stat;

private:
	bound_struct boundAtTime(double time) const;
	bound_struct boundAtIter(int i) const;
	bound_struct boundAtTimeInterpolated(double time) const;
	void fillVector(double* vector, int length, double value);
	void fillVector(int* vector, int length, int value);
	bool populated(const double* vector)const ;
	bool populated(const int* vector)const ;
};


#endif

