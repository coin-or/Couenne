/* $Id$
 *
 * Name:    tracerreport.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tracer.hpp>

#define TIE_PERC 5.0


#define MAX_ENTRIES 50
#define MAX_INSTANCES 1000
#define SEPARATORS "| \t"

#define BOUNDOPTRLTVALUES_DOUBLE_FIELDS 3

#define INVALID_ENTRY -999999999

enum BoundOptRltValuesEntries
{
	BOUNDOPT_BOUND = 0, // SDP+RLT approximate bound (used instead of BOUNDOPT_OPT when the latter is not available)
	BOUNDOPT_OPT = 1,
	BOUNDOPT_RLT = 2
};

class bounds_opt_rlt_values {
private:
	int _bound_opt_rlt_values_num_instances;
	char **_instance_name_lookup;
	double **_instance_data;

public:
	bounds_opt_rlt_values(const char* filename) {
		_instance_name_lookup = (char**) malloc(sizeof(char*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++)
			_instance_name_lookup[i] = (char*) malloc(sizeof(char)*1024);
		_instance_data = (double**) malloc(sizeof(double*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++) {
			_instance_data[i] = (double*) malloc(sizeof(double)*BOUNDOPTRLTVALUES_DOUBLE_FIELDS);
		}
		int status;
		status = readBoundsOptValuesFile(filename);
		if (status)
			exit(1);
	}

	~bounds_opt_rlt_values() {
		for(int i=0;i<MAX_INSTANCES;i++) {
			free(_instance_name_lookup[i]);
			free(_instance_data[i]);
		}
		free(_instance_name_lookup);
		free(_instance_data);
	}
/***********************************************************************/
	int getNumInstances() {return _bound_opt_rlt_values_num_instances;}
/***********************************************************************/
	double getDoubleField(int inst,int field) {return _instance_data[inst][field];}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
	int findInstanceIdx(const char *name) {
		for(int i=0;i<getNumInstances();i++) {
			if(strcmp(name,_instance_name_lookup[i]) == 0)
				return i;
		}
		return -1;
	}
/***********************************************************************/
	char *getInstanceName(int inst) {
		if ((inst >= 0)&& (inst <= getNumInstances()))
			return _instance_name_lookup[inst];
		else
			return NULL;
	}
/*********************************************************************/
	void fprint(FILE *out) {
		for (int i=0;i<getNumInstances();i++) {
			fprintf(out,"%13s\t%10.2f\t%10.2f\t%10.2f\n"
				,getInstanceName(i)
				,getDoubleField(i,BOUNDOPT_BOUND)
				,getDoubleField(i,BOUNDOPT_OPT)
				,getDoubleField(i,BOUNDOPT_RLT));
		}
	}
/***********************************************************************/
private:
	int readBoundsOptValuesFile(const char *filename) {
		FILE *_bound_opt_file = fopen(filename,"r");
		if (_bound_opt_file == NULL) {
			printf("Error opening %s. Exiting.\n",filename);
			return 1;
		}
		char line [ 1024 ];
		int linecnt = 0;
		while ( fgets ( line, sizeof line, _bound_opt_file) != NULL ) {
			if (linecnt > 0) {  //skips the first line in the file (header)
				char *curr;
				curr = strtok ( line, " \t" );
				strcpy(_instance_name_lookup[linecnt-1],curr);
				for (int j=0;j<BOUNDOPTRLTVALUES_DOUBLE_FIELDS;j++) {
					curr = strtok ( NULL, " \t" );
					if (curr == NULL) {
						printf("readBoundsOptValuesFile::missing field %d line [%d]: %s\n"
							,j,linecnt+1,line);
						return 1;
					}

 					_instance_data[linecnt-1][j] = atof(curr);
				}
			}
			linecnt++;
			if (linecnt == MAX_INSTANCES) {
				printf("readBoundsOptRltValuesFile::instances limit exceeded [%d]\n",MAX_INSTANCES);
				return 1;
			}
		}
		_bound_opt_rlt_values_num_instances = linecnt-1;
		fclose(_bound_opt_file);
		return 0;
	}


};


/*********************************************************************/



class GlobalReportInstanceData {
public:
	char 	*instanceName;
	const	char* getInstanceName() {return instanceName;}
	void	setInstanceName(const char *instanceNameParam){ strcpy(instanceName,instanceNameParam);}


	int	boundAtIter_entries;
	double	*boundAtIter_bound;
	double	*boundAtIter_time;
	int	getBoundAtIter_entries() {return boundAtIter_entries;}
	double	getBoundAtIter_bound(int i) {return boundAtIter_bound[i];}
	double	getBoundAtIter_time(int i) {return boundAtIter_time[i];}
	void	setBoundAtIter_bound(int i,double bound) {boundAtIter_bound[i] = bound;}
	void	setBoundAtIter_time(int i,double time) {boundAtIter_time[i] = time;}
	void	setBoundAtIter_entries(int n) {boundAtIter_entries = n;}

	int	boundAtTime_entries;
	double	*boundAtTime_time;
	double	*boundAtTime_bound;
	int	getBoundAtTime_entries() {return boundAtTime_entries;}
	double	getBoundAtTime_time(int i) {return boundAtTime_time[i];}
	double	getBoundAtTime_bound(int i) {return boundAtTime_bound[i];}
 	void	setBoundAtTime_time(int i,double time) {boundAtTime_time[i] = time;}
	void	setBoundAtTime_bound(int i,double bound) {boundAtTime_bound[i] = bound;}
	void	setBoundAtTime_entries(int n) {boundAtTime_entries = n;}

	int	boundInterpAtTime_entries;
	double	*boundInterpAtTime_bound;
	int	getBoundInterpAtTime_entries() {return boundInterpAtTime_entries;}
	double	getBoundInterpAtTime_bound(int i) {return boundInterpAtTime_bound[i];}
	void	setBoundInterpAtTime_bound(int i,double bound) {boundInterpAtTime_bound[i] = bound;}
	void	setBoundInterpAtTime_entries(int n) {boundInterpAtTime_entries = n;}


	GlobalReportInstanceData() {
		instanceName = new char[1024];
		boundAtTime_entries = 0;
		boundAtTime_time = new double[MAX_ENTRIES];
		boundAtTime_bound = new double[MAX_ENTRIES];
		boundAtIter_entries = 0;
		boundAtIter_time = new double[MAX_ENTRIES];
		boundAtIter_bound = new double[MAX_ENTRIES];
		boundInterpAtTime_entries = 0;
		boundInterpAtTime_bound = new double[MAX_ENTRIES];
	}

	~GlobalReportInstanceData() {
		delete [] instanceName;
		delete [] boundAtTime_time;
		delete [] boundAtTime_bound;
		delete [] boundAtIter_time;
		delete [] boundAtIter_bound;
		delete [] boundInterpAtTime_bound;
	}
};



class GlobalReport {
public:
	GlobalReportInstanceData *instances;
	GlobalReportInstanceData *getInstances() {return instances;}

	int	version;
	int	getVersion() {return version;}
	int	instances_card;
	int	getNumInstances() {return instances_card;}
	char	*valid_instances;

	int	boundAtIter_iterheader[MAX_ENTRIES];
	int	boundAtIter_entries;
	int	*getBoundAtIter_iterheader() {return boundAtIter_iterheader;}
	int	getBoundAtIter_entries() {return boundAtIter_entries;}

	double	boundAtTime_timeheader[MAX_ENTRIES];
	int	boundAtTime_entries;
	double	*getBoundAtTime_timeheader() {return boundAtTime_timeheader;}
	int	getBoundAtTime_entries() {return boundAtTime_entries;}

	double	boundInterpAtTime_timeheader[MAX_ENTRIES];
	int	boundInterpAtTime_entries;
	double	*getBoundInterpAtTime_timeheader() {return boundInterpAtTime_timeheader;}
	int	getBoundInterpAtTime_entries() {return boundInterpAtTime_entries;}


	int	boundAtIter_firstEntryIndex;
	int	boundAtTime_firstEntryIndex;
	int	boundInterpAtTime_firstEntryIndex;

	double getDoubleParam(const char *str) const{
		if (strncmp(str,REPORT_NO_ENTRY,1) == 0) {
			return TRACER_INVALID_ENTRY;
		} else {
			return atof(str);
		}
	}

	int findInstanceIdx(const char* instanceName) {
		for (int i=0;i<getNumInstances();i++) {
			if (strcmp(instanceName,instances[i].getInstanceName()) == 0)
				return i;
		}
		return -1;
	}

	GlobalReport(const char* filename, const char* valid_instances) {
		boundAtIter_entries = 0;
		boundAtTime_entries = 0;
		boundInterpAtTime_entries = 0;

		boundAtIter_firstEntryIndex = -1;
		boundAtTime_firstEntryIndex = -1;
		boundInterpAtTime_firstEntryIndex = -1;

		int version = 0;
		char line[64000];

		instances = new GlobalReportInstanceData[MAX_INSTANCES];

		FILE *globalreportfile = fopen(filename,"r");
		int linecnt = 0;
		int instanceidx = 0;
		while ( fgets ( line, sizeof line, globalreportfile) != NULL ) {
			linecnt++;
			if (linecnt == 1) {
				sscanf(line,"TRACER_VERSION %d",&version);
				continue;
			}
			if (linecnt < 10) continue; // skip initial header
			int tokidx = 0;
			char *tok;
			if (linecnt == 10) {
				//column names (we need to parse the bound@iter%d ...)
				tok = strtok(line,SEPARATORS);
				while (tok != NULL) {
					if (strncmp(tok,"bound@iter",10) == 0) {
						int curriter;
						sscanf(tok,"bound@iter%d",&curriter);
						boundAtIter_iterheader[boundAtIter_entries] = curriter;
						if (boundAtIter_firstEntryIndex < 0)
							boundAtIter_firstEntryIndex = tokidx;
					}
					if (strncmp(tok,"time@iter",9) == 0) {
						int curriter;
						sscanf(tok,"time@iter%d",&curriter);
						if (boundAtIter_iterheader[boundAtIter_entries] 
							!= curriter) {
							printf("Error parsing global report header\n");
							exit(1);
						}
						boundAtIter_entries++;
					}
					if (strncmp(tok,"bound@time",10) == 0) {
						double currtime;
						tok += 10;
						currtime = atof(tok);
						boundAtTime_timeheader[boundAtTime_entries] = currtime;
						if (boundAtTime_firstEntryIndex < 0)
							boundAtTime_firstEntryIndex = tokidx;
					}
					if (strncmp(tok,"time@time",9) == 0) {
						double currtime;
						tok += 9;
						currtime = atof(tok);
						if (boundAtTime_timeheader[boundAtTime_entries]
							!= currtime){
							printf("Error parsing global report header\n");
							exit(1);
						}
						boundAtTime_entries++;
					}
					if (strncmp(tok,"bound_interp@time",17) == 0) {
						double currtime;
						tok += 17;
						currtime = atof(tok);
						boundInterpAtTime_timeheader[boundInterpAtTime_entries] = currtime;
						boundInterpAtTime_entries++;
						if (boundInterpAtTime_firstEntryIndex < 0)
							boundInterpAtTime_firstEntryIndex = tokidx;
					}
					tok = strtok(NULL,SEPARATORS);
					tokidx ++;
				}
			} else {
				if (	(boundAtIter_firstEntryIndex<2) || 
					(boundAtTime_firstEntryIndex<2) ||
					(boundInterpAtTime_firstEntryIndex <2) ) {
					printf("Error parsing global report header - missing fields\n");
					exit(1);
				}
				tokidx = 0;
				int paramidx = 0;
				tok = strtok(line,SEPARATORS);
				if (	(valid_instances != NULL) &&
					(strstr(valid_instances,tok) ==NULL) )
					continue; // skip current instance, not in valid_instances
				double *tempparam = new double[1000];
				while (tok != NULL) {
					if (tokidx == 0) {
						instances[instanceidx].setInstanceName(tok);
					}
					if (	(tokidx >= boundAtIter_firstEntryIndex) 
						&& (tokidx < boundInterpAtTime_firstEntryIndex +
								boundInterpAtTime_entries) 
					   ){
						tempparam[paramidx] = getDoubleParam(tok);
						paramidx++;
					}
					tok = strtok(NULL,SEPARATORS);
					tokidx ++;
				}
				paramidx = 0;
				for (int i=0;i<boundAtIter_entries;i++) {
					instances[instanceidx]
						.setBoundAtIter_bound(i,tempparam[paramidx++]);
					instances[instanceidx]
						.setBoundAtIter_time(i,tempparam[paramidx++]);
				}
				for (int i=0;i<boundAtTime_entries;i++) {
					instances[instanceidx]
						.setBoundAtTime_bound(i,tempparam[paramidx++]);
					instances[instanceidx]
						.setBoundAtTime_time(i,tempparam[paramidx++]);
				}
				for (int i=0;i<boundInterpAtTime_entries;i++) {
					instances[instanceidx]
						.setBoundInterpAtTime_bound(i,tempparam[paramidx++]);
				}
				instances[instanceidx].setBoundAtIter_entries
					(boundAtIter_entries);
				instances[instanceidx].setBoundAtTime_entries
					(boundAtTime_entries);
				instances[instanceidx].setBoundInterpAtTime_entries
					(boundInterpAtTime_entries);
				instanceidx++;
//for (int i=0;i<2*boundAtIter_entries+2*boundAtTime_entries+boundInterpAtTime_entries;i++)
//printf("%.0f ",tempparam[i]);
//printf("\n");
				delete [] tempparam; 
			}
		}
		instances_card = instanceidx;
		fclose(globalreportfile);

//		for (int i=0;i<boundAtIter_entries;i++)
//			printf("%d ",boundAtIter_iterheader[i]);
//		printf("\n");
//		for (int i=0;i<boundAtTime_entries;i++)
//			printf("%.1f ",boundAtTime_timeheader[i]);
//		printf("\n");
//		for (int i=0;i<boundInterpAtTime_entries;i++)
//			printf("%.1f ",boundInterpAtTime_timeheader[i]);
//		printf("\n");
//		printf("%d %d %d\n",boundAtIter_firstEntryIndex,boundAtTime_firstEntryIndex,boundInterpAtTime_firstEntryIndex);

/*
	for (int i=0;i<instances_card;i++) {
		printf("instance %s\n",instances[i].getInstanceName());
		for (int j=0;j<boundAtIter_entries;j++) {
			printf("bound@iter%-4d= %12.4f  time@iter%-4d= %12.4f\n",
				boundAtIter_iterheader[j],
				instances[i].getBoundAtIter_bound(j),
				boundAtIter_iterheader[j],
				instances[i].getBoundAtIter_time(j));
		}
		for (int j=0;j<boundAtTime_entries;j++) {
			printf("bound@time%-8.2f= %12.4f  time@time%-8.2f= %12.4f\n",
				boundAtTime_timeheader[j],
				instances[i].getBoundAtTime_bound(j),
				boundAtTime_timeheader[j],
				instances[i].getBoundAtTime_time(j));
		}
		for (int j=0;j<boundAtTime_entries;j++) {
			printf("bound_interp@time%-8.2f= %12.4f\n",
				boundInterpAtTime_timeheader[j],
				instances[i].getBoundInterpAtTime_bound(j));
		}
	}
*/
	}
private:
	
};


#define A_BETTER 1
#define B_BETTER 2
#define AB_TIE   3
#define A_UNC    4
#define B_UNC    5
#define AB_UNC   6
int compare(double a, double b) {
	if ((a == TRACER_INVALID_ENTRY) && (b == TRACER_INVALID_ENTRY))
		return AB_UNC;
	if (a == TRACER_INVALID_ENTRY)
		return A_UNC;
	if (b == TRACER_INVALID_ENTRY)
		return B_UNC;
	if (fabs(a-b)/fabs(a) * 100.0 < TIE_PERC)
		return AB_TIE;
	if (a < b)
		return A_BETTER;
	else
		return B_BETTER;
}

int compareImpr(double a, double b) {
	if ((a == TRACER_INVALID_ENTRY) && (b == TRACER_INVALID_ENTRY))
		return AB_UNC;
	if (a == TRACER_INVALID_ENTRY)
		return A_UNC;
	if (b == TRACER_INVALID_ENTRY)
		return B_UNC;
	if (fabs(a-b) < TIE_PERC)
		return AB_TIE;
	if (a < b)
		return A_BETTER;
	else
		return B_BETTER;
}

int compare_strict(double a ,double b) {
	if ((a == TRACER_INVALID_ENTRY) && (b == TRACER_INVALID_ENTRY))
		return AB_UNC;
	if (a == TRACER_INVALID_ENTRY)
		return A_UNC;
	if (b == TRACER_INVALID_ENTRY)
		return B_UNC;
	if (a == b)
		return AB_TIE;
	if (a < b)
		return A_BETTER;
	else
		return B_BETTER;
}


double computeImprovement(double rlt, double opt, double sdp, double bound) {
	double improvement;
	double optapprox = opt;
	if (optapprox == INVALID_ENTRY){
		optapprox = sdp;
	}
	if ((optapprox != INVALID_ENTRY) && (rlt != INVALID_ENTRY) && (bound != INVALID_ENTRY)) {
		improvement = 100.0*((bound - rlt)/(optapprox - rlt));
	} else
		improvement = INVALID_ENTRY;
	return improvement;
}


int main (int argc, const char **argv) {
	printf("Comparison tolerance %.2f%%\n",(double)TIE_PERC);

	if (argc < 3) {
		printf("Usage: %s globalreport1.txt globalreport2.txt [instances.txt can be nullfile] [bounds+opt+rlt.txt]\n",argv[0]);
		return 1;
	}
	
	char *valid_instances = NULL;
	if (argc > 3) {
		FILE *valid_instances_file;
		valid_instances_file = fopen(argv[3],"r");
		if (valid_instances_file != NULL) {
			char valid_instances_tmp[64000];
			valid_instances = valid_instances_tmp;
			fgets ( valid_instances_tmp, sizeof valid_instances_tmp, valid_instances_file);
			fclose(valid_instances_file);
		}
	}

	bounds_opt_rlt_values *borv = NULL;
	if (argc > 4) {
		borv = new bounds_opt_rlt_values(argv[4]);
	}

	GlobalReport *gr1 = new GlobalReport(argv[1],valid_instances);
	GlobalReport *gr2 = new GlobalReport(argv[2],valid_instances);

	if (gr1->getVersion() != gr2->getVersion()) {
		printf("version mismatch\n");
		exit(1);
	}

	// determine how many instances are not in common
	int notincommon1 = 0;
	int notincommon2 = 0;
	for (int i=0;i<gr1->getNumInstances();i++) {
		int inst_idx1,inst_idx2;
		inst_idx1 = i;
		const char *instance_name = gr1->getInstances()[inst_idx1].getInstanceName();
		inst_idx2 = gr2->findInstanceIdx(instance_name);
		if (inst_idx2 < 0)
			notincommon1++;
	}
	for (int i=0;i<gr2->getNumInstances();i++) {
		int inst_idx1,inst_idx2;
		inst_idx2 = i;
		const char *instance_name = gr2->getInstances()[inst_idx2].getInstanceName();
		inst_idx1 = gr1->findInstanceIdx(instance_name);
		if (inst_idx1 < 0)
			notincommon2++;
	}


	// determine aggregate results for boundAtIter
	int	*boundAtIter_best1 = new int[gr1->getBoundAtIter_entries()];
	int	*boundAtIter_best2 = new int[gr1->getBoundAtIter_entries()];
	int	*boundAtIter_tie   = new int[gr1->getBoundAtIter_entries()];
	int	*boundAtIter_total = new int[gr1->getBoundAtIter_entries()];
	Stat	*improvement1vs2_stat = new Stat[gr1->getBoundAtIter_entries()];	
	for (int j=0;j<gr1->getBoundAtIter_entries();j++) {
		boundAtIter_best1[j] = TRACER_INVALID_ENTRY;
		boundAtIter_best2[j] = TRACER_INVALID_ENTRY;
		boundAtIter_tie[j]   = TRACER_INVALID_ENTRY;
		boundAtIter_total[j] = TRACER_INVALID_ENTRY;
	}

	for (int j=0;j<gr1->boundAtIter_entries;j++) {
		int total_instances = 0;
		int best1 = 0;
		int best2 = 0;
		int ties = 0;
		for (int i=0;i<gr1->getNumInstances();i++) {
			int inst_idx1,inst_idx2;
			inst_idx1 = i;
			const char *instance_name = gr1->getInstances()[inst_idx1].getInstanceName();
			inst_idx2 = gr2->findInstanceIdx(instance_name);
			if (inst_idx2 < 0) 
				continue;
			total_instances++;
			double bound1 = gr1->getInstances()[inst_idx1].getBoundAtIter_bound(j);
			double bound2 = gr2->getInstances()[inst_idx2].getBoundAtIter_bound(j);
double impr1 = INVALID_ENTRY;
double impr2 = INVALID_ENTRY;
if (borv) {
int borv_idx = borv->findInstanceIdx(instance_name);
if (borv_idx < 0) {
	printf("CANNOT FIND INSTANCE\n");
	exit(-1);
}
impr1 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound1);
impr2 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound2);
if ((impr1 != INVALID_ENTRY) && (impr2 != INVALID_ENTRY)) {
	improvement1vs2_stat[j].addEntry(impr2-impr1);

//if ((gr1->getBoundAtIter_iterheader()[j]>= 300) && (impr2-impr1 < 0)) {
//printf("iter=%d impr2VSimpr1=%.2f  %s\n",gr1->getBoundAtIter_iterheader()[j],impr2-impr1,instance_name);
//}
}
}
//			switch(compare(bound1,bound2)) {
			switch(compareImpr(impr2,impr1)) {
				case A_BETTER:
						best1++;
						break;
				case B_BETTER:
						best2++;
						break;
				case AB_TIE:
						ties++;
						break;
				case A_UNC:
						break;
				case B_UNC:
						break;
				case AB_UNC:
						break;
			}
		}
		boundAtIter_best1[j] = best1;
		boundAtIter_best2[j] = best2;
		boundAtIter_tie[j]   = ties;
		boundAtIter_total[j] = total_instances;
	}


//	printf("%13s %6s %6s %6s %6s %6s %6s %4s %4s %4s %4s %4s %6s %8s %8s %8s %8s\n",
//		"b@iter","best1","%best1","best2","%best2","ties","%ties","unc1","tot1","unc2","tot2","best","%gain","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	printf("%13s %6s %6s %6s %6s %6s %6s %6s %6s %4s %8s %8s %8s %8s\n",
		"b@iter","best1","%best1","best2","%best2","ties","%ties","unc","%unc","tot","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	for (int j=0;j<gr1->boundAtIter_entries;j++) {
//		int total_instances_cmp = boundAtIter_best1[j] + boundAtIter_best2[j] + boundAtIter_tie[j];
		int total_instances_cmp = boundAtIter_total[j];
		double best1_perc = 0.0;
		double best2_perc = 0.0;
		double tie_perc = 0.0;
		double unc_perc = 0.0;
		int unc = 0;
		if (total_instances_cmp != 0) {
			best1_perc = boundAtIter_best1[j]*100.0/total_instances_cmp;
			best2_perc = boundAtIter_best2[j]*100.0/total_instances_cmp;
			tie_perc = boundAtIter_tie[j]*100.0/total_instances_cmp;
			unc = boundAtIter_total[j] - boundAtIter_tie[j] -boundAtIter_best2[j] - boundAtIter_best1[j];
			unc_perc = unc*100/total_instances_cmp;
		}
 		printf("%13d %6d %6.2f %6d %6.2f %6d %6.2f %6d %6.2f %4d %8.2f %8.2f %8.2f %8.2f\n",
			gr1->getBoundAtIter_iterheader()[j],
			boundAtIter_best1[j],
			best1_perc,
			boundAtIter_best2[j],
			best2_perc,
			boundAtIter_tie[j],
			tie_perc,
			unc,
			unc_perc,
			total_instances_cmp,
			improvement1vs2_stat[j].mean(),
			improvement1vs2_stat[j].stdDev(),
			improvement1vs2_stat[j].min(),
			improvement1vs2_stat[j].max()
			);
	}

//#########################
printf("\n");
	// determine aggregate results for boundAtTime
	int	*boundAtTime_best1 = new int[gr1->getBoundAtTime_entries()];
	int	*boundAtTime_best2 = new int[gr1->getBoundAtTime_entries()];
	int	*boundAtTime_tie   = new int[gr1->getBoundAtTime_entries()];
	int	*boundAtTime_unc1  = new int[gr1->getBoundAtTime_entries()];
	int	*boundAtTime_unc2  = new int[gr1->getBoundAtTime_entries()];
	int	*boundAtTime_total = new int[gr1->getBoundAtTime_entries()];
	double	*boundAtTime_sumbounds1 = new double[gr1->getBoundAtTime_entries()];
	double	*boundAtTime_sumbounds2 = new double[gr1->getBoundAtTime_entries()];

	for (int j=0;j<gr1->getBoundAtTime_entries();j++) {
		improvement1vs2_stat[j].reset();
		boundAtTime_best1[j] = TRACER_INVALID_ENTRY;
		boundAtTime_best2[j] = TRACER_INVALID_ENTRY;
		boundAtTime_tie[j]   = TRACER_INVALID_ENTRY;
		boundAtTime_unc1[j]  = TRACER_INVALID_ENTRY;
		boundAtTime_unc2[j]  = TRACER_INVALID_ENTRY;
		boundAtTime_total[j] = TRACER_INVALID_ENTRY;
		boundAtTime_sumbounds1[j] = TRACER_INVALID_ENTRY;
		boundAtTime_sumbounds2[j] = TRACER_INVALID_ENTRY;
	}

	for (int j=0;j<gr1->boundAtTime_entries;j++) {
		int total_instances = 0;
		int best1 = 0;
		int best2 = 0;
		int ties = 0;
		for (int i=0;i<gr1->getNumInstances();i++) {
			int inst_idx1,inst_idx2;
			inst_idx1 = i;
			const char *instance_name = gr1->getInstances()[inst_idx1].getInstanceName();
			inst_idx2 = gr2->findInstanceIdx(instance_name);
			if (inst_idx2 < 0) 
				continue;
			total_instances++;
			double bound1 = gr1->getInstances()[inst_idx1].getBoundAtTime_bound(j);
			double bound2 = gr2->getInstances()[inst_idx2].getBoundAtTime_bound(j);
double impr1 = INVALID_ENTRY;
double impr2 = INVALID_ENTRY;
if (borv) {
int borv_idx = borv->findInstanceIdx(instance_name);
if (borv_idx < 0) {
	printf("CANNOT FIND INSTANCE\n");
	exit(-1);
}
impr1 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound1);
impr2 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound2);
if ((impr1 != INVALID_ENTRY) && (impr2 != INVALID_ENTRY)) {
	improvement1vs2_stat[j].addEntry(impr2-impr1);
}
}
//			switch(compare(bound1,bound2)) {
			switch(compare(impr2,impr1)) {
				case A_BETTER:
						best1++;
						break;
				case B_BETTER:
						best2++;		
						break;
				case AB_TIE:
						ties++;
						break;
				case A_UNC:
						break;
				case B_UNC:
						break;
				case AB_UNC:
						break;
			}
		}
		boundAtTime_best1[j] = best1;
		boundAtTime_best2[j] = best2;
		boundAtTime_tie[j]   = ties;
		boundAtTime_total[j] = total_instances;
	}


//	printf("%13s %6s %6s %6s %6s %6s %6s %4s %4s %4s %4s %4s %6s %8s %8s %8s %8s\n",
//		"b@time","best1","%best1","best2","%best2","ties","%ties","unc1","tot1","unc2","tot2","best","%gain","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	printf("%13s %6s %6s %6s %6s %6s %6s %6s %6s %4s %8s %8s %8s %8s\n",
		"b@time","best1","%best1","best2","%best2","ties","%ties","unc","%unc","tot","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	for (int j=0;j<gr1->boundAtTime_entries;j++) {
//		int total_instances_cmp = boundAtTime_best1[j] + boundAtTime_best2[j] + boundAtTime_tie[j];
		int total_instances_cmp = boundAtTime_total[j];
		double best1_perc = 0.0;
		double best2_perc = 0.0;
		double tie_perc = 0.0;
		double unc_perc = 0.0;
		int unc = 0;
		if (total_instances_cmp != 0) {
			best1_perc = boundAtTime_best1[j]*100.0/total_instances_cmp;
			best2_perc = boundAtTime_best2[j]*100.0/total_instances_cmp;
			tie_perc = boundAtTime_tie[j]*100.0/total_instances_cmp;
			unc = boundAtTime_total[j] - boundAtTime_tie[j] - boundAtTime_best2[j] - boundAtTime_best1[j];
			unc_perc = unc*100.0/total_instances_cmp;
		}
 		printf("%13.2f %6d %6.2f %6d %6.2f %6d %6.2f %6d %6.2f %6d %8.2f %8.2f %8.2f %8.2f\n",
			gr1->getBoundAtTime_timeheader()[j],
			boundAtTime_best1[j],
			best1_perc,
			boundAtTime_best2[j],
			best2_perc,
			boundAtTime_tie[j],
			tie_perc,
			unc,
			unc_perc,
			total_instances_cmp,
			improvement1vs2_stat[j].mean(),
			improvement1vs2_stat[j].stdDev(),
			improvement1vs2_stat[j].min(),
			improvement1vs2_stat[j].max());
	}


//#########################
printf("\n");
	// determine aggregate results for boundInterpAtTime
	int	*boundInterpAtTime_best1 = new int[gr1->getBoundInterpAtTime_entries()];
	int	*boundInterpAtTime_best2 = new int[gr1->getBoundInterpAtTime_entries()];
	int	*boundInterpAtTime_tie   = new int[gr1->getBoundInterpAtTime_entries()];
	int	*boundInterpAtTime_unc1  = new int[gr1->getBoundInterpAtTime_entries()];
	int	*boundInterpAtTime_unc2  = new int[gr1->getBoundInterpAtTime_entries()];
	int	*boundInterpAtTime_total = new int[gr1->getBoundInterpAtTime_entries()];
	double	*boundInterpAtTime_sumbounds1 = new double[gr1->getBoundInterpAtTime_entries()];
	double	*boundInterpAtTime_sumbounds2 = new double[gr1->getBoundInterpAtTime_entries()];
	for (int j=0;j<gr1->getBoundInterpAtTime_entries();j++) {
		improvement1vs2_stat[j].reset();
		boundInterpAtTime_best1[j] = TRACER_INVALID_ENTRY;
		boundInterpAtTime_best2[j] = TRACER_INVALID_ENTRY;
		boundInterpAtTime_tie[j]   = TRACER_INVALID_ENTRY;
		boundInterpAtTime_unc1[j]  = TRACER_INVALID_ENTRY;
		boundInterpAtTime_unc2[j]  = TRACER_INVALID_ENTRY;
		boundInterpAtTime_total[j] = TRACER_INVALID_ENTRY;
		boundInterpAtTime_sumbounds1[j] = TRACER_INVALID_ENTRY;
		boundInterpAtTime_sumbounds2[j] = TRACER_INVALID_ENTRY;
	}

	for (int j=0;j<gr1->boundInterpAtTime_entries;j++) {
		int total_instances = 0;
		int best1 = 0;
		int best2 = 0;
		int ties = 0;
		for (int i=0;i<gr1->getNumInstances();i++) {
			int inst_idx1,inst_idx2;
			inst_idx1 = i;
			const char *instance_name = gr1->getInstances()[inst_idx1].getInstanceName();
			inst_idx2 = gr2->findInstanceIdx(instance_name);
			if (inst_idx2 < 0) 
				continue;
			total_instances++;
			double bound1 = gr1->getInstances()[inst_idx1].getBoundInterpAtTime_bound(j);
			double bound2 = gr2->getInstances()[inst_idx2].getBoundInterpAtTime_bound(j);
double impr1 = INVALID_ENTRY;
double impr2 = INVALID_ENTRY;
if (borv) {
int borv_idx = borv->findInstanceIdx(instance_name);
if (borv_idx < 0) {
	printf("CANNOT FIND INSTANCE\n");
	exit(-1);
}
impr1 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound1);
impr2 = computeImprovement(
	borv->getDoubleField(borv_idx,BOUNDOPT_RLT),
	borv->getDoubleField(borv_idx,BOUNDOPT_OPT),
	borv->getDoubleField(borv_idx,BOUNDOPT_BOUND),
	bound2);
if ((impr1 != INVALID_ENTRY) && (impr2 != INVALID_ENTRY)) {
	improvement1vs2_stat[j].addEntry(impr2-impr1);
}
}
//			switch(compare(bound1,bound2)) {
			switch(compare(impr2,impr1)) {
				case A_BETTER:
						best1++;
						break;
				case B_BETTER:
						best2++;		
						break;
				case AB_TIE:
						ties++;
						break;
				case A_UNC:
						break;
				case B_UNC:
						break;
				case AB_UNC:
						break;
			}
		}
		boundInterpAtTime_best1[j] = best1;
		boundInterpAtTime_best2[j] = best2;
		boundInterpAtTime_tie[j]   = ties;
		boundInterpAtTime_total[j] = total_instances;
	}


//	printf("%13s %6s %6s %6s %6s %6s %6s %4s %4s %4s %4s %4s %6s %8s %8s %8s %8s\n",
//		"biterp@time","best1","%best1","best2","%best2","ties","%ties","unc1","tot1","unc2","tot2","best","%gain","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	printf("%13s %6s %6s %6s %6s %6s %6s %6s %6s %4s %8s %8s %8s %8s\n",
		"binterp@time","best1","%best1","best2","%best2","ties","%ties","unc","%unc","tot","2VS1mean","2VS1sdev","2VS1min","2VS1max");
	int totbest1 = 0;
	int totbest2 = 0;
	int totties = 0;
	int totunc1 = 0;
	int totunc2 = 0;
	int totcmp = 0;
	for (int j=0;j<gr1->boundInterpAtTime_entries;j++) {
//		int total_instances_cmp = boundInterpAtTime_best1[j] + boundInterpAtTime_best2[j] + boundInterpAtTime_tie[j];
		int total_instances_cmp = boundInterpAtTime_total[j];
		double best1_perc = 0.0;
		double best2_perc = 0.0;
		double tie_perc = 0.0;
		int unc = 0;
		double unc_perc = 0.0;
		if (total_instances_cmp != 0) {
			best1_perc = boundInterpAtTime_best1[j]*100.0/total_instances_cmp;
			best2_perc = boundInterpAtTime_best2[j]*100.0/total_instances_cmp;
			tie_perc = boundInterpAtTime_tie[j]*100.0/total_instances_cmp;
			unc = boundInterpAtTime_total[j] - boundInterpAtTime_tie[j] - boundInterpAtTime_best2[j] - boundInterpAtTime_best1[j];
			unc_perc = unc*100.0/total_instances_cmp;
		}
// 		printf("%13.2f %6d %6.2f %6d %6.2f %6d %6.2f %4d %4d %4d %4d %4c %6.2f %8.2f %8.2f %8.2f %8.2f\n",
 		printf("%13.2f %6d %6.2f %6d %6.2f %6d %6.2f %6d %6.2f %6d %8.2f %8.2f %8.2f %8.2f\n",
			gr1->getBoundInterpAtTime_timeheader()[j],
			boundInterpAtTime_best1[j],
			best1_perc,
			boundInterpAtTime_best2[j],
			best2_perc,
			boundInterpAtTime_tie[j],
			tie_perc,
			unc,
			unc_perc,
			total_instances_cmp,
			improvement1vs2_stat[j].mean(),
			improvement1vs2_stat[j].stdDev(),
			improvement1vs2_stat[j].min(),
			improvement1vs2_stat[j].max());

		totbest1 += boundInterpAtTime_best1[j];
		totbest2 += boundInterpAtTime_best2[j];
		totties  += boundInterpAtTime_tie[j];
		totunc1  += boundInterpAtTime_unc1[j];
		totunc2  += boundInterpAtTime_unc2[j];
		totcmp   += total_instances_cmp;
	}
	double totbest1_perc = totbest1 *100.0/totcmp;
	double totbest2_perc = totbest2 *100.0/totcmp;
	double totties_perc = totties *100.0/totcmp;
	printf("%13s %6d %6.2f %6d %6.2f %6d %6.2f\n",
		"Aggregate",
		totbest1,
		totbest1_perc,
		totbest2,
		totbest2_perc,
		totties,
		totties_perc
		);




}



