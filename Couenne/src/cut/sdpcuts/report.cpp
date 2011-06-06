/* $Id$
 *
 * Name:    report.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>

#include "CouenneConfig.h"
#include "CoinFinite.hpp" //defines COIN_DBL_MAX

#define MAX_INSTANCES 1000
#define SOLUTIONS_DOUBLE_FIELDS 8
#define BOUNDOPTVALUES_DOUBLE_FIELDS 2
#define F_RES_DOUBLE_FIELDS 4
#define F_RES_INT_FIELDS 3
#define F_RES_MAX_ITER 1000
#define F_RES_HEADER_LINES 7
#define TEST_NOT_RUN_STRING "#"
#define F_RES_MULTIPLIER -1.0

#define CMP_TOLERANCE 0.1

#define INVALID_ENTRY -999999999

enum SolutionEntries
{
	SOLUTIONS_RLT=0,
	SOLUTIONS_OPT=1,
	SOLUTIONS_V1GAP=2,
	SOLUTIONS_V2GAP=3,
	SOLUTIONS_V3GAP=4,
	SOLUTIONS_V1TIME=5,
	SOLUTIONS_V2TIME=6,
	SOLUTIONS_V3TIME=7
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


class solutions {
private:
	int _solutions_num_instances;
	char **_instance_name_lookup;
	double **_instance_data;

public:
	solutions(const char* filename) {
		_instance_name_lookup = (char**) malloc(sizeof(char*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++)
			_instance_name_lookup[i] = (char*) malloc(sizeof(char)*1024);
		_instance_data = (double**) malloc(sizeof(double*)*MAX_INSTANCES);
		for(int i=0;i<MAX_INSTANCES;i++) {
			_instance_data[i] = (double*) malloc(sizeof(double)*SOLUTIONS_DOUBLE_FIELDS);
		}
		int status;
		status = readSolutionFile(filename);
		if (status)
			exit(1);
	}
/***********************************************************************/
	~solutions() {
		for(int i=0;i<MAX_INSTANCES;i++) {
			free(_instance_name_lookup[i]);
			free(_instance_data[i]);
		}
		free(_instance_name_lookup);
		free(_instance_data);
	}
/***********************************************************************/
	int getNumInstances() {return _solutions_num_instances;}
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
/***********************************************************************/
	void fprint(FILE *out) {
		for (int i=0;i<getNumInstances();i++) {
			fprintf(out,"%13s\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f\t%7.2f\t%7.2f\t%7.2f\n"
				,getInstanceName(i)
				,getDoubleField(i,SOLUTIONS_RLT)
				,getDoubleField(i,SOLUTIONS_OPT)
				,getDoubleField(i,SOLUTIONS_V1GAP)
				,getDoubleField(i,SOLUTIONS_V2GAP)
				,getDoubleField(i,SOLUTIONS_V3GAP)
				,getDoubleField(i,SOLUTIONS_V1TIME)
				,getDoubleField(i,SOLUTIONS_V2TIME)
				,getDoubleField(i,SOLUTIONS_V3TIME));
		}
	}
/***********************************************************************/
private:
	int readSolutionFile(const char *filename) {
		FILE *_solutions_file = fopen(filename,"r");
		if (_solutions_file == NULL) {
			printf("Error opening %s. Exiting.\n",filename);
			return 1;
		}
		char line [ 1024 ];
		int linecnt = 0;
		while ( fgets ( line, sizeof line, _solutions_file) != NULL ) {
			if (linecnt > 0) {  //skips the first line in the file (header)
				char *curr;
				curr = strtok ( line, " \t" );
				strcpy(_instance_name_lookup[linecnt-1],curr);
				for (int j=0;j<SOLUTIONS_DOUBLE_FIELDS;j++) {
					curr = strtok ( NULL, " \t" );
					if (curr == NULL) {
						printf("readSolutionFile::missing field %d line [%d]: %s\n"
							,j,linecnt+1,line);
						return 1;
					}

 					_instance_data[linecnt-1][j] = atof(curr);
				}
			}
			linecnt++;
			if (linecnt == MAX_INSTANCES) {
				printf("readSolutionFile::instances limit exceeded [%d]\n",MAX_INSTANCES);
				return 1;
			}
		}
		_solutions_num_instances = linecnt-1;
		fclose(_solutions_file);
		return 0;
	}
};



/***************************************************************************/




class dataset {
private:
	int *_f_res_iter;
	double ***_f_res_double_data;
	int ***_f_res_int_data;
	char *_f_res_name;
	char **_f_res_info; //experiment parameters defined in the f_res file header

	solutions *_sol;

	int _num_instances;

public:
	dataset(const char *filename, solutions *solptr) {
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

		_sol = solptr;
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
					_f_res_double_data[i][k][j] = -COIN_DBL_MAX;
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
		for (int i=0;i<_sol->getNumInstances();i++) {
			for (int k=0;k<getIter(i);k++) {
				fprintf(out,"%13s\t%10.2f\t%3d\t%7.2f\t%5d\t%5d"
					,_sol->getInstanceName(i)
					,getDoubleField(i,k,F_RES_UBOUND)
					,getIntField(i,k,F_RES_ITER)
					,getDoubleField(i,k,F_RES_TIME)
					,getIntField(i,k,F_RES_TOTCONS)
					,getIntField(i,k,F_RES_ITERGENCONS));
				if (getDoubleField(i,k,F_RES_CURRHEUR) <= -COIN_DBL_MAX)
					fprintf(out,"\t%10s","N/A");
				else
					fprintf(out,"\t%10.2f",getDoubleField(i,k,F_RES_CURRHEUR));

				if (getDoubleField(i,k,F_RES_BESTHEUR) <= -COIN_DBL_MAX)
					fprintf(out,"\t%10s","N/A");
				else
					fprintf(out,"\t%10.2f",getDoubleField(i,k,F_RES_BESTHEUR));
				fprintf(out,"\n");
			}
		}
	}
/***********************************************************************/
	double boundClosedGap(int inst, int iter) {
 		if ((inst < _sol->getNumInstances()) && (iter < getIter(inst)) )
			return 
			fabs(computeBoundClosedGap(getDoubleField(inst,iter,F_RES_UBOUND)
						,_sol->getDoubleField(inst,SOLUTIONS_RLT)
					,_sol->getDoubleField(inst,SOLUTIONS_OPT)));
		else {
			printf("boundClosedGapF_res1:: ERROR\n");
			return -COIN_DBL_MAX;
		}
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
					curr_instance = _sol->instancename2index(curr_instance_name);
					if (curr_instance < 0) {
						printf("readF_res_File:: instance %s not in solution file\n",curr_instance_name);
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
											-COIN_DBL_MAX;

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
											-COIN_DBL_MAX;

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
	solutions *_sol;
	dataset *_f_res1,*_f_res2;


private:
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
	void fprintfReportHeader(FILE *out) {
		fprintf(out,". . . | %s . | %s . . . . . | %s . . . . | .\n","asx",getF_res(1)->getName(),getF_res(2)->getName());
		fprintf(out,"%-13s %10s %10s | %5s %7s | %10s %5s %10s %7s %5s %5s | %10s %5s %10s %7s %5s %5s | %3s\n",
			"instance"
			,"rlt"
			,"opt"
			,"cl.gap"
			,"time"
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
/***********************************************************************/
	void fprintfRHeader(FILE *out) {
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
			sprintf(tempstring,"%s-1",getSolutions()->getInstanceName(instance));
			printf("%-15s %3s %4d %4s %8.2f\n",
				tempstring,
				"A1",
				instance,
				"0",
				getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_TIME));
		}
		if ((f_res2_iter >=0) && (f_res2_iter < getF_res(2)->getIter(instance))) {
			char tempstring[256];
			sprintf(tempstring,"%s-2",getSolutions()->getInstanceName(instance));
			printf("%-15s %3s %4d %4s %8.2f\n",
				tempstring,
				"A2",
				instance,
				"0",
				getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_TIME));
		}
	}

/***********************************************************************/
	void fprintfReportLine(FILE *out,int instance, int f_res1_iter, int f_res2_iter, char *additional_col) {
		//print instance name, rlt,opt,anureet results

		//prints an error line if rlt bounds do no match
		bool rlt_bound_error = false;
		if(getF_res(1)->getIter(instance) >= 0) {
			if (!(doubleEqual(-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND), 
					  getSolutions()->getDoubleField(instance,SOLUTIONS_RLT)))) {
//				fprintf(out,". . . . . rlt != bound @iter0 = %.2f ",-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND));
				rlt_bound_error = true;
			}
		}
		if(getF_res(2)->getIter(instance) >= 0) {
			if (!(doubleEqual(-getF_res(2)->getDoubleField(instance,0,F_RES_UBOUND), 
					  getSolutions()->getDoubleField(instance,SOLUTIONS_RLT)))) {
//				fprintf(out,"rlt != bound @iter0 = %.2f",-getF_res(1)->getDoubleField(instance,0,F_RES_UBOUND));
				rlt_bound_error = true;
			}
		}
//		if (rlt_bound_error) fprintf(out,"\n");


		fprintf(out,"%-13s %10.2f %10.2f | %5.2f %7.2f | "
				,getSolutions()->getInstanceName(instance)
				,getSolutions()->getDoubleField(instance,SOLUTIONS_RLT)
				,getSolutions()->getDoubleField(instance,SOLUTIONS_OPT)
				,getSolutions()->getDoubleField(instance,SOLUTIONS_V1GAP)
				,getSolutions()->getDoubleField(instance,SOLUTIONS_V1TIME));

		//print line from f_res1
		if ((f_res1_iter >= 0) && (f_res1_iter < getF_res(1)->getIter(instance))) {
			fprintf(out,"%10.2f %5.2f"
				,getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_UBOUND)*F_RES_MULTIPLIER
				,getF_res(1)->boundClosedGap(instance,f_res1_iter));

			if (getF_res(1)->getDoubleField(instance,f_res1_iter,F_RES_BESTHEUR) <= -COIN_DBL_MAX)
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

			if (getF_res(2)->getDoubleField(instance,f_res2_iter,F_RES_BESTHEUR) <= -COIN_DBL_MAX)
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
/***********************************************************************/
public: 
	Report(const char *solutions_filename, const char *f_res1_filename, const char *f_res2_filename) {
		_sol = new solutions(solutions_filename);
		_f_res1 = new dataset(f_res1_filename,_sol);
		_f_res2 = new dataset(f_res2_filename,_sol);
	}
/***********************************************************************/
	void compareAtLastIter(FILE *out) {
		fprintfReportHeader(out);
		
		int bound_tot = 0;
		int bound_tot_1 = 0;
		int bound_tot_2 = 0;
		int bound_tot_tail = 0;

		for(int i=0;i<getSolutions()->getNumInstances();i++) {
			char string[256] = "#";
			int f_res1_iter = getF_res(1)->getIter(i)-1;
			int f_res2_iter = getF_res(2)->getIter(i)-1;
			if ((f_res1_iter>=0)&&(f_res2_iter>=0)) {
				bound_tot++;
				switch (doubleCmp(- getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_UBOUND),
						  - getF_res(2)->getDoubleField(i,f_res2_iter,F_RES_UBOUND))) {
					case 0:
						bound_tot_tail++;
						strcpy(string,"=");
						break;
					case 1:
						bound_tot_2++;
						strcpy(string,"2");
						break;
					case 2:
						bound_tot_1++;
						strcpy(string,"1");
						break;
					default: 
						break;
				}
			}
			fprintfReportLine(out,i,f_res1_iter,f_res2_iter,string);
		}
		double bound_tot_1_perc,bound_tot_2_perc,bound_tot_tail_perc;
		bound_tot_1_perc = 100 *((double) bound_tot_1 ) / ((double) bound_tot);
		bound_tot_2_perc = 100 *((double) bound_tot_2 ) / ((double) bound_tot);
		bound_tot_tail_perc = 100 *((double) bound_tot_tail ) / ((double) bound_tot);
		fprintf(out,"\n\n");
		fprintf(out,"Results for best bound on %d significant instances.\n",bound_tot);
		fprintf(out,"  Best bound:  algorithm_1 = %d [%.2f%%]   algorithm_2 = %d [%.2f%%]  ties = %d [%.2f%%]\n",bound_tot_1,bound_tot_1_perc,bound_tot_2,bound_tot_2_perc,bound_tot_tail,bound_tot_tail_perc);
	}
/***********************************************************************/
	void compareAtIterAlgo1(FILE *out,int iter) {
		fprintfReportHeader(out);
//	compare bound at iteration
//	given an iteration number, it determines the bound in f_res1 at that iteration (or at its last iter if the case) and then finds the iteration in f_res2 with the same or better bound. Prints those results.
//	does not count zero gap things

		int time_cmp_tot = 0;
		int time_cmp_tot_1 = 0;
		int time_cmp_tot_2 = 0;
		int time_cmp_tot_tail = 0;

		for(int i=0;i<getSolutions()->getNumInstances();i++) {
			if (iter < 0)
				return;

			char string[256] = "#";
			int f_res1_iter,f_res2_iter;
			f_res1_iter = iter;

			double bound1_at_iter;
//clean out this a little

			if ((f_res1_iter >= getF_res(1)->getIter(i)-1) && (getF_res(1)->getIter(i) >= 0)) {
				f_res1_iter = getF_res(1)->getIter(i)-1;
				bound1_at_iter = getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_UBOUND);
			} else	if (getF_res(1)->getIter(i) < 0) {
					bound1_at_iter = -COIN_DBL_MAX;
					f_res1_iter = -1;
			} else {
				bound1_at_iter = getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_UBOUND);
			}

			//finds a bound in f_res2 which is equal or better than bound1_at_iter
			f_res2_iter = -1;
			for (int j=0;j<getF_res(2)->getIter(i)-1;j++) {
//				printf(":%d  alg1=%.2f alg2=%.2f\n",j,- bound1_at_iter,- getF_res(2)->getDoubleField(i,j,F_RES_UBOUND));
				if (doubleCmp(- bound1_at_iter , 
						- getF_res(2)->getDoubleField(i,j,F_RES_UBOUND)) != 2) {
					//either equal or bound2 better than bound1
					f_res2_iter = j;
					break;
				}
			}

			if ((f_res1_iter == -1) && (f_res2_iter == -1)) {
				strcpy(string,"#");
			}
			else if (f_res2_iter == -1) {
				strcpy(string,"1");
				time_cmp_tot_1++;
				time_cmp_tot++;
			} else {
				time_cmp_tot++;
				switch (doubleCmp(getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_TIME),
				                  getF_res(2)->getDoubleField(i,f_res2_iter,F_RES_TIME)) ) {
					case 0:
						time_cmp_tot_tail++;
						strcpy(string,"=");
						break;
					case 1:
						time_cmp_tot_1++;
						strcpy(string,"1");
						break;
					case 2:
						time_cmp_tot_2++;
						strcpy(string,"2");
						break;
					default: 
						break;
				}
			}

			fprintfReportLine(out,i,f_res1_iter,f_res2_iter,string);
		}

		double time_cmp_tot_1_perc,time_cmp_tot_2_perc,time_cmp_tot_tail_perc;
		time_cmp_tot_1_perc = 100 *((double) time_cmp_tot_1 ) / ((double) time_cmp_tot);
		time_cmp_tot_2_perc = 100 *((double) time_cmp_tot_2 ) / ((double) time_cmp_tot);
		time_cmp_tot_tail_perc = 100 *((double) time_cmp_tot_tail ) / ((double) time_cmp_tot);
		fprintf(out,"\n\n");
		fprintf(out,"Results for best time@bound@iter_alg1=%d on %d significant instances.\n",iter,time_cmp_tot);
		fprintf(out,"  Best time:  algorithm_1 = %d [%.2f%%]   algorithm_2 = %d [%.2f%%]  ties = %d [%.2f%%]\n",time_cmp_tot_1,time_cmp_tot_1_perc,time_cmp_tot_2,time_cmp_tot_2_perc,time_cmp_tot_tail,time_cmp_tot_tail_perc);

	}
/***********************************************************************/
	void compareAtIterAlgo1R(FILE *out,int iter) {
//	compare bound at iteration
//	given an iteration number, it determines the bound in f_res1 at that iteration (or at its last iter if the case) and then finds the iteration in f_res2 with the same or better bound. Prints those results.
//	does not count zero gap things
		
		fprintfRHeader(out);

		for(int i=0;i<getSolutions()->getNumInstances();i++) {
			if (iter < 0)
				return;

			int f_res1_iter,f_res2_iter;
			f_res1_iter = iter;

			double bound1_at_iter;
//clean out this a little

			if ((f_res1_iter >= getF_res(1)->getIter(i)-1) && (getF_res(1)->getIter(i) >= 0)) {
				f_res1_iter = getF_res(1)->getIter(i)-1;
				bound1_at_iter = getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_UBOUND);
			} else	if (getF_res(1)->getIter(i) < 0) {
					bound1_at_iter = -COIN_DBL_MAX;
					f_res1_iter = -1;
			} else {
				bound1_at_iter = getF_res(1)->getDoubleField(i,f_res1_iter,F_RES_UBOUND);
			}

			//finds a bound in f_res2 which is equal or better than bound1_at_iter
			f_res2_iter = -1;
			for (int j=0;j<getF_res(2)->getIter(i)-1;j++) {
//				printf(":%d  alg1=%.2f alg2=%.2f\n",j,- bound1_at_iter,- getF_res(2)->getDoubleField(i,j,F_RES_UBOUND));
				if (doubleCmp(- bound1_at_iter , 
						- getF_res(2)->getDoubleField(i,j,F_RES_UBOUND)) != 2) {
					//either equal or bound2 better than bound1
					f_res2_iter = j;
					break;
				}
			}

			fprintfRLine(out,i,f_res1_iter,f_res2_iter);

		}
	}
/***********************************************************************/
	dataset* getF_res(int i) {
		if (i == 1)
			return _f_res1;
		if (i == 2)
			return _f_res2;
		return NULL;
	}
/***********************************************************************/
	solutions *getSolutions() {return _sol;}
/***********************************************************************/
};




/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

#if 0
int main (int argc, const char **argv) {
	if (argc < 4) {
		printf("Usage: %s [solutionfile] [f_res1.xxx] [f_res2.xxx]\n",argv[0]);
		return 1;
	}
	Report *report = new Report(argv[1],argv[2],argv[3]);

//	report->getSolutions()->fprint(stdout);

//	report->getF_res(1)->fprint(stdout);

//	report->compareAtLastIter(stdout);

	report->compareAtIterAlgo1(stdout,10);

//	report->compareAtIterAlgo1R(stdout,5);
	
	


//	compare bound at iteration
//	given an iteration number, it determines the bound in f_res1 at that iteration (or at its last iter if the case) and then finds the iteration in f_res2 with the same or better bound. Prints those results.
//	does not count zero gap things

// comparison at time
// comparison at iteration
// comparison at x percent of sdp bound given by V1

//consider only those with gap > 0

// comparison at x percent from optimal solution

// we want results of the form:
// in x instances we found a solution h% from the optimal and so on
// show if current instance comparison is ok or not (out of parameters)
	return 0;
}



#endif
