/* $Id$
 *
 * Name:    report2.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <misc_util.hpp>
#include <OsiSolverInterface.hpp> //defines DBL_MAX


//discard time differences below this value
#define MIN_TIME_DIFFERENCE 1.0
#define CMP_TOLERANCE 0.1


/////////////////

#define MAX_INSTANCES 1000
#define SOLUTIONS_DOUBLE_FIELDS 8
#define BOUNDOPTVALUES_DOUBLE_FIELDS 2
#define F_RES_DOUBLE_FIELDS 4
#define F_RES_INT_FIELDS 3
#define F_RES_MAX_ITER 2000
#define F_RES_HEADER_LINES 7
#define TEST_NOT_RUN_STRING "#"
#define F_RES_MULTIPLIER -1.0

#define INVALID_ENTRY -999999999

enum BoundOptValuesEntries
{
	BOUNDOPT_BOUND = 0,
	BOUNDOPT_OPT = 1
};

enum f_resDoubleEntries
{
	F_RES_UBOUND=0,
	F_RES_TIME=1,
	F_RES_CURRHEUR=2,
	F_RES_BESTHEUR=3
};
enum f_resIntEntries
{
	F_RES_ITER=0,
	F_RES_TOTCONS=1,
	F_RES_ITERGENCONS=2
};





/***************************************************************************/

bool doubleEqual(double val1,double val2) {
	double a,b;
	if (val1 > val2) {
		a = val2;
		b = val1;
	} else {
		a = val1;
		b = val2;			
	}
	//a<=b
	if ((b-CMP_TOLERANCE) <= a)
		return true;
	else
		return false;
}

int doubleCmp(double val1,double val2) {
	if (doubleEqual(val1,val2))
		return 0;
	if (val1<val2)
		return 1;
	else
		return 2;
}



class bounds_opt_values {
private:
	int _bound_opt_values_num_instances;
	char **_instance_name_lookup;
	double **_instance_data;

public:
	bounds_opt_values(const char* filename) {
		_instance_name_lookup = (char**) malloc(sizeof(char*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++)
			_instance_name_lookup[i] = (char*) malloc(sizeof(char)*1024);
		_instance_data = (double**) malloc(sizeof(double*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++) {
			_instance_data[i] = (double*) malloc(sizeof(double)*BOUNDOPTVALUES_DOUBLE_FIELDS);
		}
		int status;
		status = readBoundsOptValuesFile(filename);
		if (status)
			exit(1);
	}

	~bounds_opt_values() {
		for(int i=0;i<MAX_INSTANCES;i++) {
			free(_instance_name_lookup[i]);
			free(_instance_data[i]);
		}
		free(_instance_name_lookup);
		free(_instance_data);
	}
/***********************************************************************/
	int getNumInstances() {return _bound_opt_values_num_instances;}
/***********************************************************************/
	double getDoubleField(int inst,int field) {return _instance_data[inst][field];}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
	int instancename2index(char *name) {
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
			fprintf(out,"%13s\t%10.2f\t%10.2f\n"
				,getInstanceName(i)
				,getDoubleField(i,BOUNDOPT_BOUND)
				,getDoubleField(i,BOUNDOPT_OPT));
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
				for (int j=0;j<BOUNDOPTVALUES_DOUBLE_FIELDS;j++) {
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
				printf("readBoundsOptValuesFile::instances limit exceeded [%d]\n",MAX_INSTANCES);
				return 1;
			}
		}
		_bound_opt_values_num_instances = linecnt-1;
		fclose(_bound_opt_file);
		return 0;
	}


};


/*********************************************************************/




class dataset {
private:
	int *_f_res_iter;
	double ***_f_res_double_data;
	int ***_f_res_int_data;
	char *_f_res_name;
	char **_f_res_info; //experiment parameters defined in the f_res file header

	bounds_opt_values *_bov;

	int _num_instances;

public:
	dataset(const char *filename, bounds_opt_values *bovptr) {
		_f_res_name = (char*)malloc(sizeof(char)*256);
		char *name_pos,*last_dot_pos;
		name_pos = strrchr(const_cast <char *> (filename), '/');
		if(name_pos != NULL)
			strcpy(_f_res_name, &(name_pos[1]));
		else
			strcpy(_f_res_name, filename);
		last_dot_pos = strrchr(_f_res_name, '.');
		if(last_dot_pos !=NULL)
			last_dot_pos[0] = '\0';

		_bov = bovptr;
		_f_res_double_data = (double***) malloc(sizeof(double*)*MAX_INSTANCES);
		_f_res_int_data = (int***) malloc(sizeof(int*)*MAX_INSTANCES);
		_f_res_iter = (int*) malloc(sizeof(double*)*MAX_INSTANCES);
		for (int i=0;i<MAX_INSTANCES;i++) {
			_f_res_double_data[i] = (double**) malloc(sizeof(double)*F_RES_MAX_ITER);
			_f_res_int_data[i] = (int**) malloc(sizeof(int)*F_RES_MAX_ITER);
			for(int k=0;k<F_RES_MAX_ITER;k++) {
				_f_res_double_data[i][k] = 
					(double*) malloc(sizeof(double)*F_RES_DOUBLE_FIELDS);
				_f_res_int_data[i][k] = 
					(int*) malloc(sizeof(int)*F_RES_INT_FIELDS);
				for (int j=0;j<F_RES_DOUBLE_FIELDS;j++) {
					_f_res_double_data[i][k][j] = -DBL_MAX;
				}
			}
			_f_res_iter[i] = -1;
		}

		int info_lines = F_RES_HEADER_LINES + 1;
		_f_res_info = (char**) malloc(sizeof(char*)*info_lines);
		for (int i=0;i<info_lines;i++) {
			_f_res_info[i] = (char*) malloc(sizeof(char)*1024);
			_f_res_info[i][0] = '\0';
		}

		int status;
		status = readF_resFile(filename);
		if (status)
			exit(1);
	}
/***********************************************************************/
	~dataset() {
		for (int i=0;i<MAX_INSTANCES;i++) {
			for (int j=0;j<F_RES_MAX_ITER;j++) {
				free(_f_res_double_data[i][j]);
				free(_f_res_int_data[i][j]);
			}
			free(_f_res_double_data[i]);
			free(_f_res_int_data[i]);
		}
		free(_f_res_double_data);
		free(_f_res_int_data);
		free(_f_res_name);
		int info_lines = F_RES_HEADER_LINES + 1;
		for (int i=0;i<info_lines;i++)
			free(_f_res_info[i]);
		free(_f_res_info);
	}

/***********************************************************************/
	char *getName() {return _f_res_name;}
/***********************************************************************/
	int getIter(int inst) {return _f_res_iter[inst];}
/***********************************************************************/
	double getDoubleField(int inst,int iter, int field) {return _f_res_double_data[inst][iter][field];}
/***********************************************************************/
	int getIntField(int inst,int iter, int field) {return _f_res_int_data[inst][iter][field];}
/***********************************************************************/
	void fprintInfo(FILE *out) {
		for(int i=0;i<F_RES_HEADER_LINES;i++)
			fprintf(out,"%s",_f_res_info[i]);
	}
/***********************************************************************/
	void fprint(FILE* out) {
		for (int i=0;i<_bov->getNumInstances();i++) {
			for (int k=0;k<getIter(i);k++) {
				fprintf(out,"%13s\t%10.2f\t%3d\t%7.2f\t%5d\t%5d"
					,_bov->getInstanceName(i)
					,getDoubleField(i,k,F_RES_UBOUND)
					,getIntField(i,k,F_RES_ITER)
					,getDoubleField(i,k,F_RES_TIME)
					,getIntField(i,k,F_RES_TOTCONS)
					,getIntField(i,k,F_RES_ITERGENCONS));
				if (getDoubleField(i,k,F_RES_CURRHEUR) <= -DBL_MAX)
					fprintf(out,"\t%10s","N/A");
				else
					fprintf(out,"\t%10.2f",getDoubleField(i,k,F_RES_CURRHEUR));

				if (getDoubleField(i,k,F_RES_BESTHEUR) <= -DBL_MAX)
					fprintf(out,"\t%10s","N/A");
				else
					fprintf(out,"\t%10.2f",getDoubleField(i,k,F_RES_BESTHEUR));
				fprintf(out,"\n");
			}
		}
	}

/***********************************************************************/
	double timeAtBound(int inst,double value) {
		for (int i=0;i<getIter(inst);i++) {
			if (	(doubleCmp(getDoubleField(inst,i,F_RES_UBOUND),value) == 1) ||
				(doubleCmp(getDoubleField(inst,i,F_RES_UBOUND),value) == 0)) {
//printf("\n f_res%s i=%d bound=%.2f  time=%.2f value=%.2f\n",_f_res_name,i,getDoubleField(inst,i,F_RES_UBOUND),getDoubleField(inst,i,F_RES_TIME),value);
				return getDoubleField(inst,i,F_RES_TIME);
			}
		}
		return -1.0;
	}
/***********************************************************************/
	double timeAtLinearizedBound(int inst,double value) {
		double time1=-DBL_MAX;
		double time2=-DBL_MAX;
		double bound1 = -DBL_MAX;
		double bound2 = -DBL_MAX;
		int iter;
		for (int i=0;i<getIter(inst);i++) {
			if (getDoubleField(inst,i,F_RES_UBOUND) <= value) {
				time2 = getDoubleField(inst,i,F_RES_TIME);
				bound2 = getDoubleField(inst,i,F_RES_UBOUND);
				iter = i;
				break;
			}
		}
		if ((time2 != -DBL_MAX) && (iter >= 1)) {
			time1 = getDoubleField(inst,iter-1,F_RES_TIME);
			bound1 = getDoubleField(inst,iter-1,F_RES_UBOUND);
			double linearized_time;
			linearized_time = ((value - bound1)*(time2 - time1)/(bound2 - bound1)) + time1;


//printf("\ni=%d bound1=%.4f  time1=%.4f  bound2=%.4f time2=%.4f  boundinput=%.4f timelinearized=%.4f\n",
//iter,bound1,time1,bound2,time2,value,linearized_time);
			return linearized_time;
		}
		return -1.0;
	}
/***********************************************************************/
private:
	bool isNumber(char *string) {
		if (string == NULL)
			return false;
		for (int i=0;i< (int)strlen(string);i++) {
			char c=string[i];
			if (!( (c == '.') || (c == '-') || (c == '+') ||(isdigit(c)) ))
				return false;
		}
		return true;
	}
/***********************************************************************/
	int readF_resFile(const char *filename) {
		FILE *f_res_file = fopen(filename,"r");
		if (f_res_file == NULL) {
			printf("Error opening %s. Exiting.\n",filename);
			return 1;
		}
		char line [ 1024 ];
		int linecnt = 0;
		int curr_instance = -1;
		int curr_iter = 0;
		while ( fgets ( line, sizeof line, f_res_file) != NULL ) {
			if (linecnt > F_RES_HEADER_LINES) {
				char *curr;
				char curr_instance_name[256];
				curr = strtok ( line, " \t" );
				if (!(isNumber(curr))) {
					strcpy(curr_instance_name,curr);
					curr_instance = _bov->instancename2index(curr_instance_name);
					if (curr_instance < 0) {
						printf("readF_res_File:: instance %s not in bounds+opt values file\n",curr_instance_name);
						return 1;
					}
					curr_iter = 0;
					if (_f_res_iter[curr_instance] > 0) {
						printf("readF_res_File:: WARNING multiple runs of instance %s (replacing)\n",curr_instance_name);

					}
					_f_res_iter[curr_instance] = 0;
					curr = strtok ( NULL, " \t" );
				}
				if (curr_iter > F_RES_MAX_ITER - 1) {
					printf("readF_res_File::iterations for %s greater than %d\n",curr_instance_name,F_RES_MAX_ITER);
				}
				
				if (curr_instance == -1) {
						printf("readF_res_File::missing instance name on line [%d]: %s\n"
							,linecnt+1,line);
						return 1;
				}

				if (curr == NULL) {
					printf("readF_res_File::missing field 0 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				_f_res_double_data[curr_instance][curr_iter][F_RES_UBOUND] = atof(curr);

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 1 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				_f_res_int_data[curr_instance][curr_iter][F_RES_ITER] = atoi(curr);
				if (curr_iter != _f_res_int_data[curr_instance][curr_iter][F_RES_ITER]) {
					printf("readF_res_File::iteration jump line [%d]:%s\n",linecnt+1,line);
					return 1;
				}

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 2 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				_f_res_double_data[curr_instance][curr_iter][F_RES_TIME] = atof(curr);

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 4 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				_f_res_int_data[curr_instance][curr_iter][F_RES_TOTCONS] = atoi(curr);

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 5 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				_f_res_int_data[curr_instance][curr_iter][F_RES_ITERGENCONS] = atoi(curr);

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 6 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				if (isNumber(curr))
					_f_res_double_data[curr_instance][curr_iter][F_RES_CURRHEUR] =
											atof(curr);
				else
					_f_res_double_data[curr_instance][curr_iter][F_RES_CURRHEUR] =
											-DBL_MAX;

				curr = strtok ( NULL, " \t" );
				if (curr == NULL) {
					printf("readF_res_File::missing field 7 line [%d]: %s\n"
						,linecnt+1,line);
					return 1;
				}
				if (isNumber(curr))
					_f_res_double_data[curr_instance][curr_iter][F_RES_BESTHEUR] =
											atof(curr);
				else
					_f_res_double_data[curr_instance][curr_iter][F_RES_CURRHEUR] =
											-DBL_MAX;

				curr_iter++;
				_f_res_iter[curr_instance]++;
			} else {
				strcpy(_f_res_info[linecnt],line);
			}
			linecnt++;
		}
		fclose(f_res_file);
		return 0;
	}
/***********************************************************************/
	double computeBoundClosedGap(double bound, double rlt, double opt) {
		return ((-bound - rlt) / (opt-rlt))*100;
	}
/***********************************************************************/
};



/***********************************************************************/
/***********************************************************************/
/***********************************************************************/



class Report {

private:
	bounds_opt_values *_bov;
	dataset *_f_res1,*_f_res2;

private:

/*
	void fprintfReportHeader(FILE *out) {
		fprintf(out,". . . | %s . . . . . | %s . . . . | .\n",getF_res(1)->getName(),getF_res(2)->getName());
		fprintf(out,"%-13s %10s %10s | %10s %5s %10s %7s %5s %5s | %10s %5s %10s %7s %5s %5s | %3s\n",
			"instance"
			,"bound"
			,"opt"
			,"bound"
			,"cl.gap"
			,"heursol"
			,"time"
			,"iter"
			,"gencuts"
			,"bound"
			,"cl.gap"
			,"heursol"
			,"time"
			,"iter"
			,"gencuts"
			,"cmp");
	}
*/
/***********************************************************************/
/*	void fprintfRHeader(FILE *out) {
		fprintf(out,"%-15s %3s %4s %4s %8s\n",
			"",
			"alg",
			"prob",
			"fail",
			"test_val");
	}
	void fprintfRLine(FILE *out,int instance, int f_res1_iter, int f_res2_iter) {
		if ((f_res1_iter >=0) && (f_res1_iter < getF_res(1)->getIter(instance))) {
			char tempstring[256];
			sprintf(tempstring,"%s-1",getBoundOptValues()->getInstanceName(instance));
			printf("%-15s %3s %4d %4s %8.2f\n",
				tempstring,
				"A1",
				instance,
				"0",
				getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_TIME));
		}
		if ((f_res2_iter >=0) && (f_res2_iter < getF_res(2)->getIter(instance))) {
			char tempstring[256];
			sprintf(tempstring,"%s-2",getBoundOptValues()->getInstanceName(instance));
			printf("%-15s %3s %4d %4s %8.2f\n",
				tempstring,
				"A2",
				instance,
				"0",
				getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_TIME));
		}
	}
*/
/***********************************************************************/
/*
	void fprintfReportLine(FILE *out,int instance, int f_res1_iter, int f_res2_iter, char *additional_col) {
		//print instance name, rlt,opt,anureet results

		//prints an error line if rlt bounds do no match
		bool rlt_bound_error = false;
		if(getF_res(1)->getIter(instance) >= 0) {
			if (!(doubleEqual(-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND), 
					  getBoundOptValues()->getDoubleField(instance,BOUNDOPT_BOUND)))) {
//				fprintf(out,". . . . . rlt != bound @iter0 = %.2f ",-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND));
				rlt_bound_error = true;
			}
		}
		if(getF_res(2)->getIter(instance) >= 0) {
			if (!(doubleEqual(-getF_res(2)->getDoubleField(instance,0,F_RES_UBOUND), 
					  getBoundOptValues()->getDoubleField(instance,BOUNDOPT_BOUND)))) {
//				fprintf(out,"rlt != bound @iter0 = %.2f",-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND));
				rlt_bound_error = true;
			}
		}
//		if (rlt_bound_error) fprintf(out,"\n");


		fprintf(out,"%-13s %10.2f %10.2f | "
				,getBoundOptValues()->getInstanceName(instance)
				,getBoundOptValues()->getDoubleField(instance,BOUNDOPT_BOUND)
				,getBoundOptValues()->getDoubleField(instance,BOUNDOPT_OPT));

		//print line from f_res1
		if ((f_res1_iter >= 0) && (f_res1_iter < getF_res(1)->getIter(instance))) {
			fprintf(out,"%10.2f %5.2f"
				,getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_UBOUND)*F_RES_MULTIPLIER
				,getF_res(1)->boundClosedGap(instance,f_res1_iter));

			if (getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_BESTHEUR) <= -DBL_MAX)
				fprintf(out,"\t%10s","N/A ");
			else
				fprintf(out,"\t%10.2f "
					,getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_BESTHEUR)*F_RES_MULTIPLIER);

			fprintf(out,"%7.2f %5d %5d"
				,getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_TIME)
				,getF_res(1)->getIntField(instance,f_res1_iter,F_RES_ITER)
				,getF_res(1)->getIntField(instance,f_res1_iter,F_RES_TOTCONS));
		} else {
			fprintf(out,"%10s %5s %10s %7s %5s %5s"
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING);
		}

		fprintf(out," | ");

		//print line from f_res2
		if ((f_res2_iter >= 0) && (f_res2_iter < getF_res(2)->getIter(instance))) {
			fprintf(out,"%10.2f %5.2f"
				,getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_UBOUND)*F_RES_MULTIPLIER
				,getF_res(2)->boundClosedGap(instance,f_res2_iter));

			if (getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_BESTHEUR) <= -DBL_MAX)
				fprintf(out,"\t%10s","N/A ");
			else
				fprintf(out,"\t%10.2f "
					,getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_BESTHEUR)*F_RES_MULTIPLIER);

			fprintf(out,"%7.2f %5d %5d"
				,getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_TIME)
				,getF_res(2)->getIntField(instance,f_res2_iter,F_RES_ITER)
				,getF_res(2)->getIntField(instance,f_res2_iter,F_RES_TOTCONS));
		} else {
			fprintf(out,"%10s %5s %10s %7s %5s"
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING
				,TEST_NOT_RUN_STRING);
		}

		//print comparison column
		fprintf(out," | %3s\n",additional_col);
	}
*/
/***********************************************************************/
public: 
	Report(const char *bound_opt_values_filename, const char *f_res1_filename, const char *f_res2_filename) {
		_bov = new bounds_opt_values(bound_opt_values_filename);
		_f_res1 = new dataset(f_res1_filename,_bov);
		_f_res2 = new dataset(f_res2_filename,_bov);
	}
/***********************************************************************/

/***********************************************************************/
	dataset* getF_res(int i) {
		if (i == 1)
			return _f_res1;
		if (i == 2)
			return _f_res2;
		return NULL;
	}
/***********************************************************************/
	bounds_opt_values *getBoundOptValues() {return _bov;}
/***********************************************************************/
	void fprint(FILE *out) {
		getBoundOptValues()->fprint(out);
		getF_res(1)->fprint(out);
		getF_res(2)->fprint(out);
	}


	void report(double perc) {
		Stat improvement_stat;
		printf("time improvement when bound get below %.2f x SDP+RLT_BOUND\n",perc);
		printf("%14s %13s %13s %13s %12s %12s %10s\n",
			"instance","rlt_bound","sdp_bound","bound_threshold","time1","time2","improvement");
		for(int i=0;i<getBoundOptValues()->getNumInstances();i++) {
			char* instancename = getBoundOptValues()->getInstanceName(i);
			double rlt_bound = getF_res(1)->getDoubleField(i,0,F_RES_UBOUND);
			// just to check :)
			if (rlt_bound != getF_res(2)->getDoubleField(i,0,F_RES_UBOUND)) {
				printf("RLT bounds mismatch for %s!!!\n",instancename);
				exit(0);
			}
			double sdp_bound = getBoundOptValues()->getDoubleField(i,BOUNDOPT_BOUND);
			printf("%14s ",instancename);
			if (rlt_bound == -DBL_MAX) {
				printf("%13s ","N/A");
			}
			else {
				double bound_threshold = sdp_bound + ((perc)*(rlt_bound - sdp_bound));
//				double time1 = getF_res(1)->timeAtBound(i,bound_threshold);
//				double time2 = getF_res(2)->timeAtBound(i,bound_threshold);	
				double time1 = getF_res(1)->timeAtLinearizedBound(i,bound_threshold);
				double time2 = getF_res(2)->timeAtLinearizedBound(i,bound_threshold);

				printf("%13.4f ",rlt_bound);
				printf("%13.4f ",sdp_bound);
				printf("%13.4f ",bound_threshold);
				if (time1 >= 0)	printf("%12.2f ",time1);
				else		printf("%12s ","N/A");
				if (time2 >= 0)	printf("%12.2f ",time2);
				else		printf("%12s ","N/A");
				if ((time1 >= 0) && (time2 >= 0)) {
					double maxtime = time1;
					if (time1 < time2)
						maxtime = time2;
					if (maxtime == 0)
						maxtime = 0.000000001;
					double improvement = ((time1-time2)*100.0/maxtime);
					if (fabs(time1-time2) < MIN_TIME_DIFFERENCE)
						improvement = 0.0;
					printf("%10.2f %% ",improvement);
					improvement_stat.addEntry(improvement);
				}
			}
			printf("\n");
		}
		
		printf("\n Mean=%.2f %%  Stddev= %2f  Min= %.2f %% Max= %.2f %%\n",improvement_stat.mean(),improvement_stat.stdDev(),improvement_stat.min(),improvement_stat.max());
	}


	void reportMaxDiff() {
		Stat improvement_stat;
		printf("%14s %13s %13s %13s %12s %12s %10s %10s\n",
			"instance","rlt_bound","sdp_bound","bound","time1","time2","maximprtimetime","maximprtimebnd");
		for(int i=0;i<getBoundOptValues()->getNumInstances();i++) {
			char* instancename = getBoundOptValues()->getInstanceName(i);
			double rlt_bound = getF_res(1)->getDoubleField(i,0,F_RES_UBOUND);
			// just to check :)
			if (rlt_bound != getF_res(2)->getDoubleField(i,0,F_RES_UBOUND)) {
				printf("RLT bounds mismatch for %s!!!\n",instancename);
				exit(0);
			}
			double sdp_bound = getBoundOptValues()->getDoubleField(i,BOUNDOPT_BOUND);
			printf("%14s ",instancename);
			if (rlt_bound == -DBL_MAX) {
				printf("%12s ","N/A");
			}
			else {
				double maximprtime = INVALID_ENTRY;
				double maximprbound = INVALID_ENTRY;
				double besttime1 = -INVALID_ENTRY;
				double besttime2 = -INVALID_ENTRY;
				double bestbound = -INVALID_ENTRY;
				for (int j=0;j<getF_res(1)->getIter(i);j++) {
					double time1 = getF_res(1)->getDoubleField(i,j,F_RES_TIME);
					double bound1 = getF_res(1)->getDoubleField(i,j,F_RES_UBOUND);
					double time2 = getF_res(2)->timeAtLinearizedBound(i,bound1);
					if ((time1 >= 0) && (time2 >= 0)) {
						double maxtime = time1;
						if (time1 < time2)
							maxtime = time2;
						if (maxtime == 0)
							maxtime = 0.000000001;
						double improvement = ((time1-time2)*100.0/maxtime);
						if (fabs(time1-time2) < MIN_TIME_DIFFERENCE)
							improvement = 0.0;
						if (improvement > maximprtime) {
							maximprtime = improvement;
							besttime1 = time1;
							besttime2 = time2;
							bestbound = bound1;
							double currbound = (bestbound - sdp_bound)*100.0/(rlt_bound - sdp_bound);
//							if (fabs(currbound) < 0.01)
//								currbound = fabs(currbound);
							maximprbound = currbound;
							
						}
					}
				}
				printf("%13.4f ",rlt_bound);
				printf("%13.4f ",sdp_bound);
				if ((besttime1 >= 0) && (besttime2 >= 0) && (maximprtime != INVALID_ENTRY)){
					printf("%13.4f ",bestbound);
					printf("%12.2f ",besttime1);
					printf("%12.2f ",besttime2);
					printf("%10.2f %% ",maximprtime);
					printf("%10.2f %% ",maximprbound);
					improvement_stat.addEntry(maximprtime);

				} else {
					printf("%13s ","N/A");
					printf("%12s ","N/A");
					printf("%12s ","N/A");
					printf("%10s ","N/A");
					printf("%10s ","N/A");
				}
			}
			printf("\n");
		}
	}

};




/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
#if 0
int main (int argc, const char **argv) {
	if (argc < 5) {
		printf("Usage: %s [boundoptfile] [f_res1.xxx] [f_res2.xxx] [perc]\n",argv[0]);
		return 1;
	}
	Report *report = new Report(argv[1],argv[2],argv[3]);

//	report->fprint(stdout);

	report->report(atof(argv[4])/100.0);
//	report->reportMaxDiff();


}

#endif
