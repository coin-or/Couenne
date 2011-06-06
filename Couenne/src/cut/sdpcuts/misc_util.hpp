/* $Id$
 *
 * Name:    misc_util.hpp
 * Author:  Andrea Qualizza
 * Purpose: utilities for sdpcuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef MISCUTIL_HPP
#define MISCUTIL_HPP

#include <stdio.h>
#include <stdlib.h>
#include "OsiXxxSolverInterface.hpp"
#include "OsiSolverInterface.hpp"

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

class Stat {
	public: 

		Stat();
		Stat(const double *vector,int n);
		Stat(const int *vector,int n);
		void addEntry(double value);
		int numEntries() const;
		int numNZEntries() const;
		double mean() const;
		double stdDev() const;
		double sum() const;
		double min() const;
		double max() const;
		int minIndex() const;
		int maxIndex() const;
		void reset();
	private:
		int _n;
		int _nz_entries;
		double _sum;
		double _sum_squares;
		double _min;
		double _max;
		int _minIndex;
		int _maxIndex;
};

class Timer {
	public:
		Timer();
		~Timer();
		void start();
		void pause();
		void restore();
		double time();
	private:
		double _starttime;
		Timer *_pausetimer;
		bool _pause;
		double starttime();
};

extern void cpp_fprintvecINT(FILE *, char const *, const int *, int, 
			     const int);
extern void cpp_fprintvecINT(FILE *, char const *, const int *, int);
extern void cpp_printvecINT(char const *, const int *, int);

extern void cpp_printmatINT(char const *, const int **, int, int);
extern void cpp_printmatINT(char const *, int **, int, int);

extern void cpp_fprintmatINT(FILE *, char const *, const int **, int, int);
extern void cpp_fprintmatINT(FILE *, char const *, int **, int, int);

extern void cpp_printvecCHAR(char const *, const char *, int);
extern void cpp_printvecCHAR(char const *, char *, int);
extern void cpp_fprintvecCHAR(FILE *, char const *, const char *, int);
extern void cpp_fprintvecCHAR(FILE *, char const *, char *, int);

extern void cpp_fprintvecDBL(FILE *, char const *, const double *, 
			     const int, const int, char const *);
extern void cpp_fprintvecDBL(FILE *, char const *, const double *, 
			     const int, const int, const int, const int);
extern void cpp_fprintvecDBLg(FILE *, char const *, const double *, 
			      const int, const int, const int, const int);
extern void cpp_fprintvecDBL(FILE *, char const *, const double *, const int);
extern void cpp_printvecDBL(char const *, const double *, const int);
extern void cpp_printmatDBL(char const *, const double **, int, int);
extern void cpp_printmatDBL(char const *, double **, int, int);

extern void cpp_fprintvecDBL(FILE *, char const *, const double *, 
			     char const *);
extern void cpp_fprintvecDBL(FILE *, char const *, double *, int);
extern void cpp_fprintvecDBLg(FILE *, char const *, const double *, 
			      const int, const int, const int);
extern void cpp_fprintmatDBL(FILE *, char const *, const double **, int, int);
extern void cpp_fprintmatDBL(FILE *, char const *, double **, int, int);

extern void cpp_myfree_and_null(char **);

extern double cpp_genalea(int *);

extern void cpp_quicksort_dec(int, int, int *, double *);
extern void cpp_quicksortINT_dec(int, int, int *, int *);
extern void cpp_quicksort_inc(int, int, int *, double *);
extern void cpp_quicksortINT_inc(int, int, int *, int *);


extern void print_barQ(FILE *f, char *header,
		double **Q, const double *b, const int n);

extern void get_vec_from_matbar(double **mat, const int dim, double *v,
			 const int include_entry_00);

extern void get_mat_from_vec(const double *v, const int n, const double entry_00,
		      double **mat);

extern void print_mat_from_vec(FILE *f, char *header,
			const double *v, const int n, const double entry_00);

extern void get_LPsol_vec_from_vvT(const double *v, const int n, double *vec,
			    const int include_entry_00);

extern void get_mat_from_vvT(const double *v, const int n, double **mat);

extern void print_mat_from_vvT(FILE *f, char *header,
			const double *v, const int n);

extern void fprintvecmat(FILE *f, char *header, const double *v, const int n,
		  const double entry_00);

extern void check_prod_row_barQ(double **Q, const double *b, const int n,
			 const OsiSolverInterface *si,
			 const int from);

// heuristic solution from eigenvextors of Q for positive eigenvalue
extern void heur_sol_eigenv_Q(const int n, double *b, double **Q, 
		       double *heur_sol, double *heur_sol_val);

#endif
