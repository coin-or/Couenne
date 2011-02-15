/* $Id$
 *
 * Name:    tracer.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <sdpcuts.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>

Tracer::Tracer(char *instance_param, int x_vars_param, int X_vars_param, int y_vars_param) {
	iter = 0;
	x_vars = x_vars_param;
	X_vars = X_vars_param;
	y_vars = y_vars_param;
	tot_vars = x_vars + X_vars + y_vars;

	strcpy(instance,instance_param);

	fillVector(main__iteration_bound,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__iteration_time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__lp_time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__active_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__added_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__total_eigendecompositions,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(main__deleted_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(sdpcuts__num_negative_ev,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sdpcuts__most_negative_ev,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sdpcuts__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sdpcuts__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(sparsify__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sparsify__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sparsify__duplicated_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(sparsify__wise_decompositions,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(orthocut__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(orthocut__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(linquadcuts__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(linquadcuts__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(disjunctivecuts__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(disjunctivecuts__total_cuts,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	fillVector(heuristics__current_solution,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__best_solution,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__xxT_solution,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__xxT_solution_lp_heuristic_improvement,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__xxT_time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__MNLP_solution,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__MNLP_solution_lp_heuristic_improvement,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__MNLP_time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__GW_solution,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__GW_solution_lp_heuristic_improvement,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);
	fillVector(heuristics__GW_time,TRACE_ITER_SIZE,TRACER_INVALID_ENTRY);

	sparsify__nz_populated = false;
	sparsify__single_column_sparsity_populated = false;
	sparsify__column_pairs_sparsity_populated = false;
	sparsify__top20perc_cuts_violation_populated = false;
	
}

Tracer::~Tracer() {}


// control functions
void Tracer::newIter(){ 
	iter++;
}
int Tracer::iterations() const{
	return iter+1;
}

void Tracer::detailedReport() const{
	FILE *detailed_report_file = NULL;
	detailed_report_file = fopen(TRACER_DETAILED_REPORT_FILE, "r");
	if(detailed_report_file == NULL) {
		detailed_report_file = fopen(TRACER_DETAILED_REPORT_FILE, "w");
		fprintf(detailed_report_file,"TRACER_VERSION  %d\n",VERSION_TRACER);
		print_ifdefs(detailed_report_file);

		fprintf(detailed_report_file,
		"%-15s %12s %12s %8s %12s %12s %8s %8s %8s %8s %8s %1s %8s %12s %12s %8s %1s %12s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %12s %12s %12s %12s %1s %12s %8s %1s %12s %8s %1s %12s %8s %1s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n"
		,"instance"
		,"tot_time"
		,"iter_time"
		,"iter"
		,"bound"
		,"lp_time"
		,"tot_cuts"
		,"act_cuts"
		,"added_cuts"
		,"del_cuts"
		,"eigendec"
			,REPORT_DELIMITER
		,"num_neg_ev"
		,"most_neg_ev"
		,"sdpcuts_time"
		,"sdpcuts_tot_cuts"
			,REPORT_DELIMITER
		,"sparsify_time"
		,"sparsify_tot_cuts"
		,"sparsify_dup_cuts"
		,"sparsify_wise_eigendec"
		,"sparsify_nz_mean"
		,"sparsify_nz_stddev"
		,"sparsify_nz_min"
		,"sparsify_nz_max"
		,"sparsify_single_column_sparsity_mean"
		,"sparsify_single_column_sparsity_stddev"
		,"sparsify_single_column_sparsity_min"
		,"sparsify_single_column_sparsity_max"
		,"sparsify_column_pairs_sparsity_mean"
		,"sparsify_column_pairs_sparsity_stddev"
		,"sparsify_column_pairs_sparsity_min"
		,"sparsify_column_pairs_sparsity_max"
		,"sparsify_top20perc_cuts_violation_mean"
		,"sparsify_top20perc_cuts_violation_stddev"
		,"sparsify_top20perc_cuts_violation_min"
		,"sparsify_top20perc_cuts_violation_max"
			,REPORT_DELIMITER
		,"orthocut_time"
		,"orthocut_tot_cuts"
			,REPORT_DELIMITER
		,"linquadcuts_time"
		,"linquadcuts_tot_cuts"
			,REPORT_DELIMITER
		,"disjunctivecuts_time"
		,"disjunctivecuts_cuts"
			,REPORT_DELIMITER
		,"heuristics_current_solution"
		,"heuristics_best_solution"
		,"heuristics_time"
		,"heuristics_xxT_solution"
		,"heuristics_xxT_solution_lp_heuristic_improvement"
		,"heuristics_xxT_time"
		,"heuristics_MNLP_solution"
		,"heuristics_MNLP_solution_lp_heuristic_improvement"
		,"heuristics_MNLP_time"
		,"heuristics_GW_solution"
		,"heuristics_GW_solution_lp_heuristic_improvement"
		,"heuristics_GW_time");
	}
	fclose(detailed_report_file);
	detailed_report_file = fopen(TRACER_DETAILED_REPORT_FILE, "a");
	for (int i=0;i<iterations();i++) {
		fprintf(detailed_report_file,"%-15s",instance);
		if (populated(main__time))
			fprintf(detailed_report_file," %12.4f",main__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);

		if (populated(main__iteration_time))
			fprintf(detailed_report_file," %12.4f",main__iteration_time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);

		fprintf(detailed_report_file," %8d",i);

		if (populated(main__iteration_bound))
			fprintf(detailed_report_file," %12.4f",main__iteration_bound[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(main__lp_time))
			fprintf(detailed_report_file," %12.4f",main__lp_time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(main__total_cuts))
			fprintf(detailed_report_file," %8d",main__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(main__active_cuts))
			fprintf(detailed_report_file," %8d",main__active_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(main__added_cuts))
			fprintf(detailed_report_file," %8d",main__added_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(main__deleted_cuts))
			fprintf(detailed_report_file," %8d",main__deleted_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(main__total_eigendecompositions))
			fprintf(detailed_report_file," %8d",main__total_eigendecompositions[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);

		fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if (populated(sdpcuts__num_negative_ev))
			fprintf(detailed_report_file," %8d",sdpcuts__num_negative_ev[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(sdpcuts__most_negative_ev))
			fprintf(detailed_report_file," %12.4f",sdpcuts__most_negative_ev[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if(populated(sdpcuts__time))
			fprintf(detailed_report_file," %12.4f",sdpcuts__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if(populated(sdpcuts__total_cuts))
			fprintf(detailed_report_file," %8d",sdpcuts__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);

		fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if (populated(sparsify__time))
			fprintf(detailed_report_file," %12.4f",sparsify__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(sparsify__total_cuts))
			fprintf(detailed_report_file," %8d",sparsify__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(sparsify__duplicated_cuts))
			fprintf(detailed_report_file," %8d",sparsify__duplicated_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (populated(sparsify__wise_decompositions))
			fprintf(detailed_report_file," %8d",sparsify__wise_decompositions[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
		if (sparsify__nz_populated) {
			fprintf(detailed_report_file
			," %8.2f",sparsify__nz_iter_stat[i].mean());
			fprintf(detailed_report_file
			," %8.2f",sparsify__nz_iter_stat[i].stdDev());
			fprintf(detailed_report_file
			," %8.0f",sparsify__nz_iter_stat[i].min());
			fprintf(detailed_report_file
			," %8.0f",sparsify__nz_iter_stat[i].max());
		}
		else {
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);	
		}
		if (sparsify__single_column_sparsity_populated) {
			fprintf(detailed_report_file
			," %8.2f",sparsify__single_column_sparsity_iter_stat[i].mean());
			fprintf(detailed_report_file
			," %8.2f",sparsify__single_column_sparsity_iter_stat[i].stdDev());
			fprintf(detailed_report_file
			," %8.0f",sparsify__single_column_sparsity_iter_stat[i].min());
			fprintf(detailed_report_file
			," %8.0f",sparsify__single_column_sparsity_iter_stat[i].max());
		}
		else {
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);	
		}
		if (sparsify__column_pairs_sparsity_populated) {
			fprintf(detailed_report_file
			," %8.2f",sparsify__column_pairs_sparsity_iter_stat[i].mean());
			fprintf(detailed_report_file
			," %8.2f",sparsify__column_pairs_sparsity_iter_stat[i].stdDev());
			fprintf(detailed_report_file
			," %8.0f",sparsify__column_pairs_sparsity_iter_stat[i].min());
			fprintf(detailed_report_file
			," %8.0f",sparsify__column_pairs_sparsity_iter_stat[i].max());
		}
		else {
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);	
		}
		if (sparsify__top20perc_cuts_violation_populated) {
			fprintf(detailed_report_file
			," %12.4f",sparsify__top20perc_cuts_violation_iter_stat[i].mean());
			fprintf(detailed_report_file
			," %12.4f",sparsify__top20perc_cuts_violation_iter_stat[i].stdDev());
			fprintf(detailed_report_file
			," %12.4f",sparsify__top20perc_cuts_violation_iter_stat[i].min());
			fprintf(detailed_report_file
			," %12.4f",sparsify__top20perc_cuts_violation_iter_stat[i].max());
		}
		else {
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);	
		}

		fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if(populated(orthocut__time))
			fprintf(detailed_report_file," %12.4f",orthocut__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if(populated(orthocut__total_cuts))
			fprintf(detailed_report_file," %8d",orthocut__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);

		fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if(populated(linquadcuts__time))
			fprintf(detailed_report_file," %12.4f",linquadcuts__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if(populated(linquadcuts__total_cuts))
			fprintf(detailed_report_file," %8d",linquadcuts__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);

			fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if(populated(disjunctivecuts__time))
			fprintf(detailed_report_file," %12.4f",disjunctivecuts__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if(populated(disjunctivecuts__total_cuts))
			fprintf(detailed_report_file," %8d",disjunctivecuts__total_cuts[i]);
		else
			fprintf(detailed_report_file," %8s",REPORT_NO_ENTRY);

		fprintf(detailed_report_file," %1s",REPORT_DELIMITER);

		if (heuristics__current_solution[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file," %12.4f",heuristics__current_solution[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__best_solution[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file," %12.4f",heuristics__best_solution[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(heuristics__time))
			fprintf(detailed_report_file," %12.4f",heuristics__time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__xxT_solution[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__xxT_solution[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__xxT_solution_lp_heuristic_improvement[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__xxT_solution_lp_heuristic_improvement[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(heuristics__xxT_time))
			fprintf(detailed_report_file
			," %12.4f",heuristics__xxT_time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__MNLP_solution[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__MNLP_solution[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__MNLP_solution_lp_heuristic_improvement[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__MNLP_solution_lp_heuristic_improvement[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(heuristics__MNLP_time))
			fprintf(detailed_report_file
			," %12.4f",heuristics__MNLP_time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__GW_solution[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__GW_solution[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (heuristics__GW_solution_lp_heuristic_improvement[i] != TRACER_INVALID_ENTRY)
			fprintf(detailed_report_file
			," %12.4f",heuristics__GW_solution_lp_heuristic_improvement[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		if (populated(heuristics__GW_time))
			fprintf(detailed_report_file
			," %12.4f",heuristics__GW_time[i]);
		else
			fprintf(detailed_report_file," %12s",REPORT_NO_ENTRY);
		
		fprintf(detailed_report_file,"\n");
	}
	fclose(detailed_report_file);
}





void Tracer::globalReport() const{
	int 	bound_at_iter[14]	= {1,2,3,5,10,15,20,30,50,100,200,300,500,1000};
	double	bound_at_time[17]	= {0.5,1,2,3,5,10,15,20,30,60,120,180,300,600,1200,1800,3600};
	double	bound_interp_at_time[17]= {0.5,1,2,3,5,10,15,20,30,60,120,180,300,600,1200,1800,3600};

	int bound_at_iter_entries = 
		sizeof(bound_at_iter)/sizeof(bound_at_iter[0]);
	int bound_at_time_entries = 
		sizeof(bound_at_time)/sizeof(bound_at_time[0]);
	int bound_interp_at_time_entries =
		sizeof(bound_interp_at_time)/sizeof(bound_interp_at_time[0]);

	FILE *global_report_file = NULL;
	global_report_file = fopen(TRACER_GLOBAL_REPORT_FILE, "r");

	if(global_report_file == NULL) {
		global_report_file = fopen(TRACER_GLOBAL_REPORT_FILE, "w");
		fprintf(global_report_file,"TRACER_VERSION  %d\n",VERSION_TRACER);
		print_ifdefs(global_report_file);

		fprintf(global_report_file
		,"%-15s %1s %8s %8s %8s %8s %1s %12s %12s %8s %12s %8s %8s %8s %8s %1s %12s %8s %1s %12s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %12s %12s %12s %12s %1s %12s %8s %1s %12s %8s %1s %12s %8s %1s %12s %12s %8s %8s %12s %8s %8s %12s %8s %8s %12s"
		,"instance"
			,REPORT_DELIMITER
		,"tot_vars"
		,"x_vars"
		,"X_vars"
		,"y_vars"
			,REPORT_DELIMITER
		,"totaltime"
		,"lpsolve_totaltime"
		,"iterations"
		,"bound"
		,"tot_cuts"
		,"generated_cuts"
		,"deleted_cuts"
		,"eigendec"
			,REPORT_DELIMITER
		,"sdpcuts_totaltime"
		,"sdpcuts_tot_cuts"
			,REPORT_DELIMITER
		,"sparsify_totaltime"
		,"sparsify_tot_cuts"
		,"sparsify_dup_cuts"
		,"sparsify_wise_eigendec"
		,"sparsify_nz_mean"
		,"sparsify_nz_stddev"
		,"sparsify_nz_min"
		,"sparsify_nz_max"
		,"sparsify_single_column_sparsity_mean"
		,"sparsify_single_column_sparsity_stddev"
		,"sparsify_single_column_sparsity_min"
		,"sparsify_single_column_sparsity_max"
		,"sparsify_column_pairs_sparsity_mean"
		,"sparsify_column_pairs_sparsity_stddev"
		,"sparsify_column_pairs_sparsity_min"
		,"sparsify_column_pairs_sparsity_max"
		,"sparsify_top20perc_cuts_violation_mean"
		,"sparsify_top20perc_cuts_violation_stddev"
		,"sparsify_top20perc_cuts_violation_min"
		,"sparsify_top20perc_cuts_violation_max"
			,REPORT_DELIMITER
		,"orthocut_totaltime"
		,"orthocut_tot_cuts"
			,REPORT_DELIMITER
		,"linquadcuts_totaltime"
		,"linquadcuts_tot_cuts"
			,REPORT_DELIMITER
		,"disjunctivecuts_totaltime"
		,"disjunctivecuts_tot_cuts"
			,REPORT_DELIMITER
		,"heuristics_best_solution"
		,"heuristics_time"
		,"heuristics_xxT_best"
		,"heuristics_xxT_lp_heuristic_improvements"
		,"heuristics_xxT_time"
		,"heuristics_MNLP_best"
		,"heuristics_MNLP_lp_heuristic_improvements"
		,"heuristics_MNLP_time"
		,"heuristics_GW_best"
		,"heuristics_GW_lp_heuristic_improvements"
		,"heuristics_GW_time"
		);

		for (int i=0;i<bound_at_iter_entries;i++) {
			char string1[20];
			char string2[20];
			sprintf(string1,"bound@iter%d",bound_at_iter[i]);
			sprintf(string2,"time@iter%d",bound_at_iter[i]);
			fprintf(global_report_file," %1s %12s %12s"
			,REPORT_DELIMITER
			,string1
			,string2);
		}
		for (int i=0;i<bound_at_time_entries;i++) {
			char string1[20];
			char string2[20];
			sprintf(string1,"bound@time%.1f",bound_at_time[i]);
			sprintf(string2,"time@time%.1f",bound_at_time[i]);
			fprintf(global_report_file," %1s %12s %12s"
			,REPORT_DELIMITER
			,string1
			,string2);
		}
		for (int i=0;i<bound_interp_at_time_entries;i++) {
			char string1[20];
			sprintf(string1,"bound_interp@time%.1f",bound_interp_at_time[i]);
			fprintf(global_report_file," %1s %12s"
			,REPORT_DELIMITER
			,string1);
		}
		fprintf(global_report_file,"\n");
	}
	
	fclose(global_report_file);
	global_report_file = fopen(TRACER_GLOBAL_REPORT_FILE, "a");

	fprintf(global_report_file,"%-15s",instance);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	fprintf(global_report_file," %8d",tot_vars);
	fprintf(global_report_file," %8d",x_vars);
	fprintf(global_report_file," %8d",X_vars);
	fprintf(global_report_file," %8d",y_vars);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (populated(main__iteration_time))
		fprintf(global_report_file," %12.4f",main__time[iter]);
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(main__lp_time))
		fprintf(global_report_file," %12.4f",main__lp_time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	fprintf(global_report_file," %8d",iter);
	if (populated(main__iteration_bound))
		fprintf(global_report_file," %12.4f",main__iteration_bound[iter]);
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	if (populated(main__total_cuts))
		fprintf(global_report_file," %8d",main__total_cuts[iter]);
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	if (populated(main__added_cuts))
		fprintf(global_report_file," %8.0f",main__added_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	if (populated(main__deleted_cuts))
		fprintf(global_report_file," %8.0f",main__deleted_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	if (populated(main__total_eigendecompositions))
		fprintf(global_report_file," %8.0f",main__total_eigendecompositions_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	
	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (populated(sdpcuts__time))
		fprintf(global_report_file," %12.4f",sdpcuts__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(sdpcuts__total_cuts))
		fprintf(global_report_file," %8.0f",sdpcuts__total_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);
	
	if (populated(sparsify__time))
		fprintf(global_report_file," %12.4f",sparsify__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(sparsify__total_cuts))
		fprintf(global_report_file
		," %8.0f",sparsify__total_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	if (populated(sparsify__duplicated_cuts))
		fprintf(global_report_file
		," %8.0f",sparsify__duplicated_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	if (populated(sparsify__wise_decompositions))
		fprintf(global_report_file
		," %8.0f",sparsify__wise_decompositions_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	if (sparsify__nz_populated) {
		fprintf(global_report_file
		," %8.2f",sparsify__nz_global_stat.mean());
		fprintf(global_report_file
		," %8.2f",sparsify__nz_global_stat.stdDev());
		fprintf(global_report_file
		," %8.0f",sparsify__nz_global_stat.min());
		fprintf(global_report_file
		," %8.0f",sparsify__nz_global_stat.max());
	}
	else {
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	}
	if (sparsify__single_column_sparsity_populated) {
		fprintf(global_report_file
		," %8.2f",sparsify__single_column_sparsity_global_stat.mean());
		fprintf(global_report_file
		," %8.2f",sparsify__single_column_sparsity_global_stat.stdDev());
		fprintf(global_report_file
		," %8.0f",sparsify__single_column_sparsity_global_stat.min());
		fprintf(global_report_file
		," %8.0f",sparsify__single_column_sparsity_global_stat.max());
	}
	else {
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	}
	if (sparsify__column_pairs_sparsity_populated) {
		fprintf(global_report_file
		," %8.2f",sparsify__column_pairs_sparsity_global_stat.mean());
		fprintf(global_report_file
		," %8.2f",sparsify__column_pairs_sparsity_global_stat.stdDev());
		fprintf(global_report_file
		," %8.0f",sparsify__column_pairs_sparsity_global_stat.min());
		fprintf(global_report_file
		," %8.0f",sparsify__column_pairs_sparsity_global_stat.max());
	}
	else {
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);
	}

	if (sparsify__top20perc_cuts_violation_populated) {
		fprintf(global_report_file
		," %12.4f",sparsify__top20perc_cuts_violation_global_stat.mean());
		fprintf(global_report_file
		," %12.4f",sparsify__top20perc_cuts_violation_global_stat.stdDev());
		fprintf(global_report_file
		," %12.4f",sparsify__top20perc_cuts_violation_global_stat.min());
		fprintf(global_report_file
		," %12.4f",sparsify__top20perc_cuts_violation_global_stat.max());
	}
	else {
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	}

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (populated(orthocut__time))
		fprintf(global_report_file," %12.4f",orthocut__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(orthocut__total_cuts))
		fprintf(global_report_file," %8.0f",orthocut__total_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (populated(linquadcuts__time))
		fprintf(global_report_file," %12.4f",linquadcuts__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(linquadcuts__total_cuts))
		fprintf(global_report_file," %8.0f",linquadcuts__total_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (populated(disjunctivecuts__time))
		fprintf(global_report_file," %12.4f",disjunctivecuts__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);
	if (populated(disjunctivecuts__total_cuts))
		fprintf(global_report_file," %8.0f",disjunctivecuts__total_cuts_global_stat.sum());
	else
		fprintf(global_report_file," %8s",REPORT_NO_ENTRY);

	fprintf(global_report_file," %1s",REPORT_DELIMITER);

	if (heuristics__best_solution[iter] != TRACER_INVALID_ENTRY)
		fprintf(global_report_file," %12.4f",heuristics__best_solution[iter]);
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	if (populated(heuristics__time))
		fprintf(global_report_file," %12.4f",heuristics__time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	int heuristics_xxT_best_count = 0;
	for (int i=0;i<iter;i++)
		if  (	(heuristics__current_solution[i] != TRACER_INVALID_ENTRY) && 
			(heuristics__xxT_solution[i] != TRACER_INVALID_ENTRY) &&
			(fabs(heuristics__xxT_solution[i]-heuristics__current_solution[i]) 
				< TRACER_HEURISTICS_COMPARISON_TOLERANCE) )
			heuristics_xxT_best_count++;
	fprintf(global_report_file
	," %8d",heuristics_xxT_best_count);
	fprintf(global_report_file," %8d"
		,heuristics__xxT_solution_lp_heuristic_improvement_global_stat.numNZEntries());

	if (populated(heuristics__xxT_time))
		fprintf(global_report_file
		," %12.4f"
		,heuristics__xxT_time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	int heuristics_MNLP_best_count = 0;
	for (int i=0;i<iterations();i++)
		if  (	(heuristics__current_solution[i] != TRACER_INVALID_ENTRY) && 
			(heuristics__MNLP_solution[i] != TRACER_INVALID_ENTRY) &&
			(fabs(heuristics__MNLP_solution[i]-heuristics__current_solution[i]) 
					< TRACER_HEURISTICS_COMPARISON_TOLERANCE) )
			heuristics_MNLP_best_count++;
	fprintf(global_report_file
	," %8d",heuristics_MNLP_best_count);

	fprintf(global_report_file," %8d"
		,heuristics__MNLP_solution_lp_heuristic_improvement_global_stat.numNZEntries());

	if (populated(heuristics__MNLP_time))
		fprintf(global_report_file
		," %12.4f"
		,heuristics__MNLP_time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	int heuristics_GW_best_count = 0;
	for (int i=0;i<iterations();i++)
		if  (	(heuristics__current_solution[i] != TRACER_INVALID_ENTRY) && 
			(heuristics__GW_solution[i] != TRACER_INVALID_ENTRY) &&
			(fabs(heuristics__GW_solution[i]-heuristics__current_solution[i]) 
				< TRACER_HEURISTICS_COMPARISON_TOLERANCE) )
			heuristics_GW_best_count++;
	fprintf(global_report_file
	," %8d",heuristics_GW_best_count);

	fprintf(global_report_file," %8d"
		,heuristics__GW_solution_lp_heuristic_improvement_global_stat.numNZEntries());

	if (populated(heuristics__GW_time))
		fprintf(global_report_file
		," %12.4f"
		,heuristics__GW_time_global_stat.sum());
	else
		fprintf(global_report_file," %12s",REPORT_NO_ENTRY);

	for (int i=0;i<bound_at_iter_entries;i++) {
		bound_struct bs = boundAtIter(bound_at_iter[i]);
		if ((bs.bound != TRACER_INVALID_ENTRY) && (bs.time != TRACER_INVALID_ENTRY)) {
			fprintf(global_report_file," %1s %12.4f %12.4f"
				,REPORT_DELIMITER
				,bs.bound
				,bs.time);
		} else {
			fprintf(global_report_file," %1s %12s %12s"
				,REPORT_DELIMITER
				,REPORT_NO_ENTRY
				,REPORT_NO_ENTRY);
		}
	}
	for (int i=0;i<bound_at_time_entries;i++) {
		bound_struct bs = boundAtTime(bound_at_time[i]);
		if ((bs.bound != TRACER_INVALID_ENTRY) && (bs.time != TRACER_INVALID_ENTRY)) {
			fprintf(global_report_file," %1s %12.4f %12.4f"
				,REPORT_DELIMITER
				,bs.bound
				,bs.time);
		} else {
			fprintf(global_report_file," %1s %12s %12s"
				,REPORT_DELIMITER
				,REPORT_NO_ENTRY
				,REPORT_NO_ENTRY);
		}
	}
	for (int i=0;i<bound_interp_at_time_entries;i++) {
		bound_struct bs = boundAtTimeInterpolated(bound_interp_at_time[i]);
		if (bs.bound != TRACER_INVALID_ENTRY)
			fprintf(global_report_file," %1s %12.4f"
				,REPORT_DELIMITER
				,bs.bound);
		else
			fprintf(global_report_file," %1s %12s"
				,REPORT_DELIMITER
				,REPORT_NO_ENTRY);
	}

	fprintf(global_report_file,"\n");
	fclose(global_report_file);
}

// functions for main entries
void Tracer::setMainBound(double value) {
	main__iteration_bound[iter] = value;
}
void Tracer::setMainIterationTime(double value) {
	main__iteration_time[iter] = value;
	if (iter == 0)
		main__time[iter] = value;
	else
		main__time[iter] = value + main__time[iter-1];
}
void Tracer::setMainLPTime(double value) {
	main__lp_time[iter] = value;
	main__lp_time_global_stat.addEntry(value);
}
void Tracer::setMainActiveCuts(int value) {
	main__active_cuts[iter] = value;
	main__active_cuts_global_stat.addEntry(value);
}
void Tracer::setMainAddedCuts(int value){
	main__added_cuts[iter] = value;
	main__added_cuts_global_stat.addEntry(value);
}
void Tracer::setMainTotalCuts(int value){
	main__total_cuts[iter] = value;
	main__total_cuts_global_stat.addEntry(value);
}
void Tracer::setMainTotalEigendecompositions(int value){
	main__total_eigendecompositions[iter] = value;
	main__total_eigendecompositions_global_stat.addEntry(value);
}
void Tracer::incrementMainTotalEigendecompositions(){
	if (main__total_eigendecompositions[iter] < 0)
		main__total_eigendecompositions[iter] = 0;
	main__total_eigendecompositions[iter]++;
	main__total_eigendecompositions_global_stat.addEntry(1);
}
void Tracer::setMainDeletedCuts(int value){
	main__deleted_cuts[iter] = value;
	main__deleted_cuts_global_stat.addEntry(value);
}


// functions for sdp cuts entries
void Tracer::setSDPNumNegativeEV(int value){
	sdpcuts__num_negative_ev[iter] = value;
}
void Tracer::setSDPMostNegativeEV(double value){
	sdpcuts__most_negative_ev[iter] = value;
}
void Tracer::setSDPCutsTime(double value) {
	sdpcuts__time_global_stat.addEntry(value);
	sdpcuts__time[iter] = value;
}

void Tracer::setSDPCutsTotalCuts(int value) {
	sdpcuts__total_cuts_global_stat.addEntry(value);
	sdpcuts__total_cuts[iter] = value;
}


// functions for sparsify entries
void Tracer::setSparsifyTime(double value) {
	sparsify__time_global_stat.addEntry(value);
	sparsify__time[iter] = value;
}
void Tracer::setSparsifyTotalCuts(int value){
	sparsify__total_cuts_global_stat.addEntry(value);
	sparsify__total_cuts[iter] = value;
}
void Tracer::setSparsifyDuplicatedCuts(int value){
	sparsify__duplicated_cuts_global_stat.addEntry(value);
	sparsify__duplicated_cuts[iter] = value;
}
void Tracer::setSparsifyWiseDecompositions(int value){
	sparsify__wise_decompositions[iter] = value;
	sparsify__wise_decompositions_global_stat.addEntry(value);
}
void Tracer::addSparsifyNz(int value){
	sparsify__nz_populated = true;
	sparsify__nz_iter_stat[iter].addEntry(value);
	sparsify__nz_global_stat.addEntry(value);
}
void Tracer::addSparsifySingleColumnSparsity(int value){
	sparsify__single_column_sparsity_populated = true;
	sparsify__single_column_sparsity_iter_stat[iter].addEntry(value);
	sparsify__single_column_sparsity_global_stat.addEntry(value);
}
void Tracer::addSparsifyColumnPairSparsity(int value){
	sparsify__column_pairs_sparsity_populated = true;
	sparsify__column_pairs_sparsity_iter_stat[iter].addEntry(value);
	sparsify__column_pairs_sparsity_global_stat.addEntry(value);
}
void Tracer::addSparsifyTop20PercCutsViolation(double value){
	sparsify__top20perc_cuts_violation_populated = true;
	sparsify__top20perc_cuts_violation_iter_stat[iter].addEntry(value);
	sparsify__top20perc_cuts_violation_global_stat.addEntry(value);
}


// functions for orthocut  entries
void Tracer::setOrthocutTime(double value){
	orthocut__time_global_stat.addEntry(value);
	orthocut__time[iter] = value;
}
void Tracer::setOrthocutTotalCuts(int value){
	orthocut__total_cuts_global_stat.addEntry(value);
	orthocut__total_cuts[iter] = value;
}

// functions for linquad cuts entries 
void Tracer::setLinquadTime(double value){
	linquadcuts__time_global_stat.addEntry(value);
	linquadcuts__time[iter] = value;
}
void Tracer::setLinquadTotalCuts(int value){
	linquadcuts__total_cuts_global_stat.addEntry(value);
	linquadcuts__total_cuts[iter] = value;
}

// functions for disjunctive cuts entries
void Tracer::setDisjunctiveCutsTime(double value){
	disjunctivecuts__time_global_stat.addEntry(value);
	disjunctivecuts__time[iter] = value;
}

void Tracer::setDisjunctiveCutsTotalCuts(int value){
	disjunctivecuts__total_cuts_global_stat.addEntry(value);
	disjunctivecuts__total_cuts[iter] = value;
}



// functions for heuristics entries
void Tracer::setHeuristicsCurrentSolution(double value){
	heuristics__current_solution[iter] = value;
}
void Tracer::setHeuristicsBestSolution(double value){
	heuristics__best_solution[iter] = value;
}
void Tracer::setHeuristicsTime(double value){
	heuristics__time[iter] = value;
	heuristics__time_global_stat.addEntry(value);
}
void Tracer::setHeuristicsxxTSolution(double value){
	heuristics__xxT_solution[iter] = value;
}
void Tracer::setHeuristicsxxTSolutionLPHeuristicImprovement(double value){
	heuristics__xxT_solution_lp_heuristic_improvement[iter] = value;
	heuristics__xxT_solution_lp_heuristic_improvement_global_stat.addEntry(value);
}
void Tracer::setHeuristicsMNLPSolution(double value){
	heuristics__MNLP_solution[iter] = value;
}
void Tracer::setHeuristicsMNLPSolutionLPHeuristicImprovement(double value){
	heuristics__MNLP_solution_lp_heuristic_improvement[iter] = value;
	heuristics__MNLP_solution_lp_heuristic_improvement_global_stat.addEntry(value);
}
void Tracer::setHeuristicsGWSolution(double value){
	heuristics__GW_solution[iter] = value;
}
void Tracer::setHeuristicsGWSolutionLPHeuristicImprovement(double value){
	heuristics__GW_solution_lp_heuristic_improvement[iter] = value;
	heuristics__GW_solution_lp_heuristic_improvement_global_stat.addEntry(value);
}
void Tracer::setHeuristicsxxTTime(double value){
	heuristics__xxT_time[iter] = value;
	heuristics__xxT_time_global_stat.addEntry(value);
}
void Tracer::setHeuristicsMNLPTime(double value){
	heuristics__MNLP_time[iter] = value;
	heuristics__MNLP_time_global_stat.addEntry(value);
}
void Tracer::setHeuristicsGWTime(double value){
	heuristics__GW_time[iter] = value;
	heuristics__GW_time_global_stat.addEntry(value);
}



// private functions

bound_struct Tracer::boundAtTime(double time) const {
	bound_struct bs;
	bs.bound = TRACER_INVALID_ENTRY;
	bs.time = TRACER_INVALID_ENTRY;
	bs.iter1 = TRACER_INVALID_ENTRY;
	bs.iter2 = TRACER_INVALID_ENTRY;

	int iter_out = 0;
	if (time > main__time[iter]) {
		return bs;
	}
	
	while (main__time[iter_out]<time) {
		iter_out++;
	}
	bs.bound = main__iteration_bound[iter_out];
	bs.time = main__time[iter_out];
	bs.iter1 = iter_out;
	bs.iter2 = iter_out;
	
	return bs;
}



bound_struct Tracer::boundAtIter(int i) const {
	bound_struct bs;
	bs.bound = TRACER_INVALID_ENTRY;
	bs.time = TRACER_INVALID_ENTRY;
	bs.iter1 = TRACER_INVALID_ENTRY;
	bs.iter2 = TRACER_INVALID_ENTRY;

	int iter_out = 0;
	if ( (i < 0) || (i > iterations()-1) ) {
		return bs;
	}
	
	bs.bound = main__iteration_bound[i];
	bs.time = main__time[i];
	bs.iter1 = i;
	bs.iter2 = i;
	
	return bs;
}

bound_struct Tracer::boundAtTimeInterpolated(double time) const {
	bound_struct bs;
	bs.bound = TRACER_INVALID_ENTRY;
	bs.time = TRACER_INVALID_ENTRY;
	bs.iter1 = TRACER_INVALID_ENTRY;
	bs.iter2 = TRACER_INVALID_ENTRY;

	int iter1 = 0;
	int iter2 = iter;
	if (time > main__time[iter]) {
		return bs;
	}
	if (time < main__time[0]) {
		return bs;
	}
	while (main__time[iter1]<time) {
		iter1++;
	}
	while (main__time[iter2]>time) {
		iter2--;
	}
	//linear interpolation using known points 
	//	(main__iteration_time[iter1], main__iteration_bound[iter1])
	//	(main__iteration_time[iter2], main__iteration_bound[iter2])	
	double interpolated_bound = main__iteration_bound[iter1] + 
			( (time - main__time[iter1]) *
			  (main__iteration_bound[iter2] - main__iteration_bound[iter1]) / 
			  (main__time[iter2] - main__time[iter1]) );
	bs.bound = interpolated_bound;
	bs.time = time;
	bs.iter1 = iter1;
	bs.iter2 = iter2;
	return bs;
}

void Tracer::fillVector(double* vector, int length, double value){
	for (int i=0;i<length;i++)
		vector[i] = value;
}

void Tracer::fillVector(int* vector, int length, int value){
	for (int i=0;i<length;i++)
		vector[i] = value;
}

bool Tracer::populated(const int* vector) const {
	for (int i=0;i<iterations();i++)
		if(vector[i] == TRACER_INVALID_ENTRY)
			return false;
	return true;
}

bool Tracer::populated(const double* vector) const {
	for (int i=0;i<iterations();i++)
		if(vector[i] == TRACER_INVALID_ENTRY)
			return false;
	return true;
}


