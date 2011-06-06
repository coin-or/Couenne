/* $Id$
 *
 * Name:    misc_util.cpp
 * Author:  Andrea Qualizza
 * Purpose: utilities for sdpcuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <cstdio>
#include <math.h>
#include "CoinTime.hpp"

#include "misc_util.hpp"

Stat::Stat() { reset(); }
Stat::Stat(const double *vector,int n) {
	reset();
	for(int i=0;i<n;i++)
		addEntry(vector[i]);
}
Stat::Stat(const int *vector,int n) {
	reset();
	for(int i=0;i<n;i++)
		addEntry((double)vector[i]);
}
void Stat::reset() { _n = 0; _nz_entries = 0; _sum = 0; _sum_squares = 0;}
void Stat::addEntry(double value) {
	if (fabs(value) > 1e-8)
		_nz_entries++;
	_sum         += value;
	_sum_squares += value*value;
	_n++;
	if (_n == 1) {
		_min = value;
		_max = value;
		_minIndex = 1;
		_maxIndex = 1;
	} else {
		if (value < _min) {
			_min = value;
 			_minIndex = _n;
		}
		if (value > _max) {
			_max = value;
			_maxIndex = _n;
		}
	}
}
int Stat::numEntries() const {return _n;}
int Stat::numNZEntries() const {return _nz_entries;}
double Stat::mean() const { 
	if (_n>0)
		return _sum/_n;
	else
		return 0.0;
}
double Stat::stdDev() const { 
	if (_n>0)
		return sqrt( fabs((_sum_squares / _n) - (mean()*mean())) ); //fabs is to ensure that when the difference is close to 0.0 the value is actually NOT negative, otherwise we can have nan errors
	else
		return 0.0;
}
double Stat::sum() const {return _sum;}
double Stat::min() const { 
	if (_n>0)
		return _min;
	else
		return 0.0;
	}
double Stat::max() const { 
	if (_n>0)
		return _max;
	else
		return 0.0;
	}
int Stat::minIndex() const { 
	if (_n>0)
		return _minIndex;
	else
		return -1;
	}
int Stat::maxIndex() const { 
	if (_n>0)
		return _maxIndex;
	else
		return -1;
	}
/**********************************************************/
Timer::Timer() { _pause = false; _starttime = -99999999;}
Timer::~Timer() { if (_pause) delete _pausetimer; }
void Timer::start() {
  _starttime = CoinCpuTime ();
}

double Timer::time() {
	if (_starttime == -99999999)
		return 0.0;
	if (_pause) {
		return fabs( _pausetimer->starttime() - _starttime );
	}

	return (fabs (CoinCpuTime () - _starttime));
}

void Timer::pause() {
	if (_pause) {
		fprintf(stderr,"timer already in pause state\n");
		fprintf(stdout,"timer already in pause state\n");
	} else {
		_pause = true;
		_pausetimer = new Timer();
		_pausetimer->start();
	}
}
void Timer::restore() {
	if (_pause) {
		_pause = false;
		_starttime += _pausetimer->time();
		delete _pausetimer;
	} else {
		fprintf(stderr,"timer not in pause state\n");
		fprintf(stdout,"timer not in pause state\n");
	}
}

double Timer::starttime() {
	return _starttime;
}
/**********************************************************/
void cpp_fprintvecINT(FILE *file, 
		      char const *vecstr, const int *x, const int n,
		      const int numberAcross)
{
  int num, fromto, upto, j, i;

  num = (n/numberAcross) + 1;
  if(vecstr != NULL) {
    fprintf(file, "%s : \n", vecstr);
  }
  for(j=0; j<num; j++){
    fromto = numberAcross * j;
    upto = numberAcross * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      fprintf(file, " %4d", x[i]);
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintvecINT */

/**********************************************************/
void cpp_fprintvecINT(FILE *file, 
		      char const *vecstr, const int *x, const int n)
{
  cpp_fprintvecINT(file, vecstr, x, n, 10);
} /* cpp_fprintvecINT */

/**********************************************************/
void cpp_printvecINT(char const *vecstr, const int *x, const int n)
{
  cpp_fprintvecINT(stdout, vecstr, x, n, 10);
} /* cpp_printvecINT */

/**********************************************************/
void cpp_fprintmatINT(FILE *file, char const *vecstr, const int **x, 
		      const int m, 
		      const int n)
{
  int i, j;

  fprintf(file, "%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(file, " %4d", x[i][j]);
    }
      fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintmatINT */

/**********************************************************/
void cpp_fprintmatINT(FILE *file, char const *vecstr, int **x, const int m, 
		  const int n)
{
  int i, j;

  fprintf(file, "%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(file, " %4d", x[i][j]);
    }
      fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintmatINT */

/**********************************************************/
void cpp_printmatINT(char const *vecstr, const int **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* cpp_printmatINT */

/**********************************************************/
void cpp_printmatINT(char const *vecstr, int **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* cpp_printmatINT */

/**********************************************************/
void cpp_printvecCHAR(char const *vecstr, const char *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10 * j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i = fromto; i<upto; i++)
      printf(" %c", x[i]);
    printf("\n");
  }
  printf("\n");
} /* cpp_printvecCHAR */

/**********************************************************/
void cpp_printvecCHAR(char const *vecstr, char *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10 * j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i = fromto; i<upto; i++)
      printf(" %c", x[i]);
    printf("\n");
  }
  printf("\n");
} /* cpp_printvecCHAR */

/**********************************************************/
void cpp_fprintvecCHAR(FILE *file, char const *vecstr, const char *x, const int n)
{
  int num, fromto, upto, j, i;

  num = n/10 + 1;
  fprintf(file, "%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10 * j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      fprintf(file, " %c", x[i]);
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintvecCHAR */

/**********************************************************/
void cpp_fprintvecCHAR(FILE *file, char const *vecstr, char *x, const int n)
{
  int num, fromto, upto, j, i;

  num = n/10 + 1;
  fprintf(file, "%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10 * j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      fprintf(file, " %c", x[i]);
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintvecCHAR */

/**********************************************************/
void cpp_fprintvecDBL(FILE *file, char const *vecstr, const double *x, 
		      const int n, const int numberAcross, 
		      char *form)
{
  int num, fromto, upto, j, i;

  num = (n/numberAcross) + 1;
  if(vecstr != NULL) {
    fprintf(file,"%s:\n",vecstr);
  }
  for(j=0;j<num; j++){
    fromto = numberAcross * j;
    upto = numberAcross * (j + 1);
    if(n <= upto) upto = n;
    for(i = fromto; i < upto; i++)
      fprintf(file, form, x[i]);
    fprintf(file,"\n");
  }
  fprintf(file,"\n");
} /* cpp_fprintvecDBL */

/**********************************************************/
void cpp_fprintvecDBL(FILE *file, char const *vecstr, const double *x, 
		      const int n, const int numberAcross, 
		      const int print_pos, const int decimals)
{
  char form[15];
  sprintf(form, " %%%d.%df", print_pos, decimals);
  cpp_fprintvecDBL(file, vecstr, x, n, numberAcross, form);
} /* cpp_fprintvecDBL */

/**********************************************************/
void cpp_fprintvecDBLg(FILE *file, char const *vecstr, const double *x, 
		       const int n, const int numberAcross, 
		       const int print_pos, const int print_digits)
{
  char form[15];
  sprintf(form, " %%%d.%dg", print_pos, print_digits);
  cpp_fprintvecDBL(file, vecstr, x, n, numberAcross, form);
} /* cpp_fprintvecDBL */

/**********************************************************/
void cpp_fprintvecDBL(FILE *file, char const *vecstr, const double *x, const int n)
{
  cpp_fprintvecDBL(file, vecstr, x, n, 10, 7, 3);
} /* cpp_fprintvecDBL */

/**********************************************************/
void cpp_printvecDBL(char const *vecstr, const double *x, const int n)
{
  cpp_fprintvecDBL(stdout, vecstr, x, n, 10, 7, 3);
} /* cpp_printvecDBL */

/**********************************************************/
void cpp_printmatDBL(char const *vecstr, const double **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4.3f", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* cpp_printmatDBL */

/**********************************************************/
void cpp_printmatDBL(char const *vecstr, double **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4.3f", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* cpp_printmatDBL */

/**********************************************************/
void cpp_fprintmatDBL(FILE *file, char const *vecstr, const double **x, const int m, 
		  const int n)
{
  int i, j;

  fprintf(file, "%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(file, " %4.3f", x[i][j]);
    }
      fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintmatDBL */

/**********************************************************/
void cpp_fprintmatDBL(FILE *file, char const *vecstr, double **x, const int m, 
		  const int n)
{
  int i, j;

  fprintf(file, "%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(file, " %4.3f", x[i][j]);
    }
      fprintf(file, "\n");
  }
  fprintf(file, "\n");
} /* cpp_fprintmatDBL */

/**********************************************************/
/***
void cpp_myfree_and_null (char **ptr)
{
  if (*ptr != NULL) {
    free(*ptr);
    *ptr = NULL;
  }
} 
***/
/* cpp_free_and_null */
/**********************************************************/
double cpp_genalea (int *x0)
{
  int m = 2147483647;
  int a = 16807 ;
  int b = 127773 ;
  int c = 2836 ;
  int x1, k;

  k = (int) ((*x0)/b) ;
  x1 = a*(*x0 - k*b) - k*c ;
  if(x1 < 0) x1 = x1 + m;
  *x0 = x1;

  return((double)x1/(double)m);
/*
  if(((double)x1/(double)m > 0.0001) && ((double)x1/(double)m < 0.99999))
    return((double)x1/(double)m);
  else return(genalea(x0));
*/
} /* cpp_genalea */


/**********************************************************/
int cpp_partition_dec(const int q, const int r, int *key, double *value)
{
  double v, t;
  int i, j, index;

  v = value[r];
  i = q-1;
  j = r;
  do{     
    i++;j--;
    while ( value[i]> v) i++;
    while ( (value[j]< v)&&(j>q) ) j--;
    t = value[i];
    value[i] = value[j];
    value[j] = t;
    index = key[i];
    key[i] = key[j];
    key[j] = index;
  }while (j>i);

  value[j] = value[i];
  value[i] = value[r];
  value[r] = t;
  key[j] = key[i];
  key[i] = key[r];
  key[r] = index;
  return(i);
}

/**********************************************************/
void cpp_quicksort_dec(const int k, const int l, int *key, double *value)
                         /* Quicksort in decreasing order  */
{
  int p, i, q, r;
  int stack[1000];

  if(l == 1) {
    return;
  }

  stack[0] = stack[1] = 0;
  q = k;
  r = l-1;
  p = 2;
  do{
    if (r>q){
      i = cpp_partition_dec(q,r,key,value);
      if (i-q>r-i){
        stack[p] = q;
        stack[p+1] = i - 1;
        q = i + 1;
      }
      else{
        stack[p] = i + 1;
        stack[p+1] = r;
        r = i - 1;
      }
      p += 2;
    }
    else{
      p -= 2;
      q = stack[p];
      r = stack[p+1];
    }
  }while (p);
}

/**********************************************************/
int cpp_partitionINT_dec(const int q, const int r, int *key, int *value)
{
  int v, t;
  int i, j, index;

  v = value[r];
  i = q-1;
  j = r;
  do{     
    i++;j--;
    while ( value[i]> v) i++;
    while ( (value[j]< v)&&(j>q) ) j--;
    t = value[i];
    value[i] = value[j];
    value[j] = t;
    index = key[i];
    key[i] = key[j];
    key[j] = index;
  }while (j>i);

  value[j] = value[i];
  value[i] = value[r];
  value[r] = t;
  key[j] = key[i];
  key[i] = key[r];
  key[r] = index;
  return(i);
}

/**********************************************************/
void cpp_quicksortINT_dec(const int k, const int l, int *key, int *value)
                         /* Quicksort in decreasing order  */
{
  int p, i, q, r;
  int stack[1000];

  if(l == 1) {
    return;
  }

  stack[0] = stack[1] = 0;
  q = k;
  r = l-1;
  p = 2;
  do{
    if (r>q){
      i = cpp_partitionINT_dec(q,r,key,value);
      if (i-q>r-i){
        stack[p] = q;
        stack[p+1] = i - 1;
        q = i + 1;
      }
      else{
        stack[p] = i + 1;
        stack[p+1] = r;
        r = i - 1;
      }
      p += 2;
    }
    else{
      p -= 2;
      q = stack[p];
      r = stack[p+1];
    }
  }while (p);
}

/**********************************************************/
int cpp_partition_inc(const int q, const int r, int *key, double *value)
{
  double v, t;
  int i, j, index;

  v = value[r];
  i = q-1;
  j = r;
  do{     
    i++;j--;
    while ( value[i]< v) i++;
    while ( (value[j]> v)&&(j>q) ) j--;
    t = value[i];
    value[i] = value[j];
    value[j] = t;
    index = key[i];
    key[i] = key[j];
    key[j] = index;
  }while (j>i);

  value[j] = value[i];
  value[i] = value[r];
  value[r] = t;
  key[j] = key[i];
  key[i] = key[r];
  key[r] = index;
  return(i);
} /* cpp_partition_inc */

/**********************************************************/
void cpp_quicksort_inc(const int k, const int l, int *key, double *value)
                         /* Quicksort in increasing order  */
{
  int p, i, q, r;
  int stack[1000];


  if(l == 1) {
    return;
  }

  stack[0] = stack[1] = 0;
  q = k;
  r = l-1;
  p = 2;
  do{
    if (r>q){
      i = cpp_partition_inc(q,r,key,value);
      if (i-q>r-i){
        stack[p] = q;
        stack[p+1] = i - 1;
        q = i + 1;
      }
      else{
        stack[p] = i + 1;
        stack[p+1] = r;
        r = i - 1;
      }
      p += 2;
    }
    else{
      p -= 2;
      q = stack[p];
      r = stack[p+1];
    }
  }while (p);
} /* cpp_quicksort_inc */

/**********************************************************/
int cpp_partitionINT_inc(const int q, const int r, int *key, int *value)
{
  int v, t;
  int i, j, index;

  v = value[r];
  i = q-1;
  j = r;
  do{     
    i++;j--;
    while ( value[i]< v) i++;
    while ( (value[j]> v)&&(j>q) ) j--;
    t = value[i];
    value[i] = value[j];
    value[j] = t;
    index = key[i];
    key[i] = key[j];
    key[j] = index;
  }while (j>i);

  value[j] = value[i];
  value[i] = value[r];
  value[r] = t;
  key[j] = key[i];
  key[i] = key[r];
  key[r] = index;
  return(i);
} /* cpp_partitionINT_inc */

/**********************************************************/
void cpp_quicksortINT_inc(const int k, const int l, int *key, int *value)
                         /* Quicksort in increasing order for integers */
{
  int p, i, q, r;
  int stack[1000];


  if(l == 1) {
    return;
  }

  stack[0] = stack[1] = 0;
  q = k;
  r = l-1;
  p = 2;
  do{
    if (r>q){
      i = cpp_partitionINT_inc(q,r,key,value);
      if (i-q>r-i){
        stack[p] = q;
        stack[p+1] = i - 1;
        q = i + 1;
      }
      else{
        stack[p] = i + 1;
        stack[p+1] = r;
        r = i - 1;
      }
      p += 2;
    }
    else{
      p -= 2;
      q = stack[p];
      r = stack[p+1];
    }
  }while (p);
} /* cpp_quicksortINT_inc */




/***********************************************************************/
void get_barQ(double **Q, const double *b, const int n, double **mat) {
  
  int np = n+1;
  mat[0][0] = 0;

  for(int i=1; i<np; i++) {
    mat[0][i] = b[i-1];
    mat[i][0] = b[i-1];
  }
  
  for(int i=1; i<np; i++) {
    for(int j=i; j<np; j++) {
      mat[i][j] = Q[i-1][j-1];
      mat[j][i] = Q[i-1][j-1];
    }
  }
} /* get_barQ */

/***********************************************************************/
void print_barQ(FILE *f, char *header,
		double **Q, const double *b, const int n) {
  
  int i, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
  }

  get_barQ(Q, b, n, mat);
  cpp_fprintmatDBL(f, header, mat, np, np);

  for (i=0; i<np; i++) {
    free (mat[i]);
  }
  free (mat);


} /* print_barQ */

/***********************************************************************/
void get_vec_from_matbar(double **mat, const int dim, double *v,
			 const int include_entry_00) {

  int ind = 0;

  if(include_entry_00) {
    v[ind] = mat[0][0];
    ind++;
  }

  for(int j=1; j<dim; j++) {
    v[ind] = mat[0][j];
    ind++;
  }
  for(int i=1; i<dim; i++) {
    for(int j=i; j<dim; j++) {
      v[ind] = mat[i][j];
      ind++;
    }
  }
} /* get_vec_from_matbar */

/***********************************************************************/
void get_mat_from_vec(const double *v, const int n, const double entry_00,
		      double **mat) {
  
  mat[0][0] = entry_00;

  for(int i=1; i<n+1; i++) {
    mat[0][i] = v[i-1];
    mat[i][0] = v[i-1];
  }
  
  for(int i=1; i<n+1; i++) {
    for(int j=i; j<n+1; j++) {
      int ind = indexQ(i-1, j-1, n);
      mat[i][j] = v[ind];
      mat[j][i] = v[ind];
    }
  }
} /* get_mat_from_vec */

/***********************************************************************/
void print_mat_from_vec(FILE *f, char *header,
			const double *v, const int n, const double entry_00) {
  
  int i, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
  }

  get_mat_from_vec(v, n, entry_00, mat);
  cpp_fprintmatDBL(f, header, mat, np, np);

  for (i=0; i<np; i++) {
    free (mat[i]);
  }
  free (mat);
} /* print_mat_from_vec */

/***********************************************************************/
void get_LPsol_vec_from_vvT(const double *v, const int n, double *vec,
			    const int include_entry_00) {

  int i, j, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
  }
  
  mat[0][0] = 1;
  
  for(i=1; i<np; i++) {
    mat[0][i] = v[i-1];
    mat[i][0] = v[i-1];
  }
  
  for(i=1; i<np; i++) {
    for(j=1; j<np; j++) {
      mat[i][j] = v[i-1] * v[j-1];
    }
  }
  get_vec_from_matbar(mat, np, vec, include_entry_00);
  
  for (i=0; i<np; i++) {
    free(mat [i]);
  }
  free(mat);
} /* get_LPsol_vec_from_vvT */

/***********************************************************************/
void get_mat_from_vvT(const double *v, const int n, double **mat) {
  
  for(int i=0; i<n+1; i++) {
    for(int j=0; j<n+1; j++) {
      mat[i][j] = v[i] * v[j];
    }
  }
} /* get_mat_from_vvT */

/***********************************************************************/
void print_mat_from_vvT(FILE *f, char *header,
			const double *v, const int n) {
  
  int i, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
  }

  get_mat_from_vvT(v, n, mat);
  cpp_fprintmatDBL(f, header, mat, np, np);

  for (i=0; i<np; i++) {
    free (mat[i]);
  }
  free (mat);
} /* print_mat_from_vvT */

/***********************************************************************/
void fprintvecmat(FILE *f, char *header, const double *v, const int n,
		  const double entry_00) {

  int i, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
  }

  get_mat_from_vec(v, n, entry_00, mat);
  cpp_fprintmatDBL(f, header, mat, np, np);

  for (i=0; i<np; i++) {
    free (mat[i]);
  }
  free (mat);

} /* printvecmat */


/***********************************************************************/
void check_prod_row_barQ(double **Q, const double *b, const int n,
			 const OsiSolverInterface *si,
			 const int from) {
  int i, j, k, np = n+1;
  double **mat = (double **) malloc (np * sizeof (double *));
  double **mat_row = (double **) malloc (np * sizeof (double *));

  for (i=0; i<np; i++) {
    mat[i] = (double *) malloc (np * sizeof (double));
    mat_row[i] = (double *) malloc (np * sizeof (double));
  }

  get_barQ(Q, b, n, mat);

  int curr_nrows = si->getNumRows();
  const double *y = si->getRowPrice();
  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  const double *rhs = si->getRightHandSide();
  
  int row_size = np * np;

  double *row_v = new double[row_size];

  for(int i = from; i<curr_nrows; i++) {
    if(fabs(y[i]) > 1e-5) {
      
      const int *rowind = byRow->getVector(i).getIndices();
      int card_row = byRow->getVector(i).getNumElements();
      const double *rowelem = byRow->getVector(i).getElements();

      for(j=0; j<row_size; j++) {
	row_v[j] = 0;
      }
      
      for(j=0; j<card_row; j++) {
	row_v[rowind[j]] = rowelem[j];
      }
      
      get_mat_from_vec(row_v, n, -rhs[i], mat_row);

      double ckmuQ = 0;
      
      for (j=0; j<np; j++) {
	for (k=0; k<np; k++) {
	  ckmuQ += mat[j][j] * mat_row[j][k]; 
	}
      }      

      if(ckmuQ > 0) {
	printf("### WARNING: SdpCutGen::check_prod_row_barQ(): row %d with dual: %8.4f  muQ: %8.4f\n", i, y[i], ckmuQ);

#ifdef TRACE_ALL
	cpp_printmatDBL("mat_row", mat_row, np, np);
	cpp_printmatDBL("mat_barQ", mat, np, np);
#endif

	exit(1);
      }
      if(ckmuQ > 1e30) {
	exit(1);
      }
    }
  }

  for (i=0; i<np; i++) {
    free (mat[i]);
    free (mat_row[i]);
  }
  free (mat);
  free (mat_row);
  delete[] row_v;

} /* check_prod_row_barQ */


