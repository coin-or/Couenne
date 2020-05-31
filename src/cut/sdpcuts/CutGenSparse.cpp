/*
 *
 * Name:    CutGenSparse.cpp
 * Authors: Andrea Qualizza
 *          Pietro Belotti
 * Purpose: Generation of sparse cuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "operators/CouenneExprMul.hpp"
#include "CouenneSdpCuts.hpp"

#include "dsyevx_wrapper.hpp"

#ifdef _WIN32
#define drand48() ((double) (rand () * (RAND_MAX + 1) + rand ()) / (RAND_MAX + 1) * (RAND_MAX + 1))
#endif

//#define DEBUG

const bool WISE_SPARSIFY = true;

#define SPARSIFY_MAX_CARD  10000
#define WISE_SPARSIFY_GAP  0.0001

#define SPARSIFY_OLD_DELTA 0.50
#define SPARSIFY_NEW_DELTA 0.50

#define SPARSIFY_OLD_NZ_THRESHOLD 0.50
#define SPARSIFY_NEW_NZ_THRESHOLD 0.70

#define EV_TOL 1e-13

using namespace Couenne;

/************************************************************************/
void CouenneSdpCuts::sparsify2 (const int n,
				const double *A, double **sparse_v_mat,
				int *card_v_mat, int min_nz, int *evdec_num) const {

  bool *del_idx = NULL;

  int 
    running_n = n,
    best_idx,
    rnsq     = (running_n - 1) * (running_n - 1),
    card_ev_best = running_n - 1;

  double

    *matrix = CoinCopyOfArray (A, n*n),

    best_val,

    *matrixCopy = CoinCopyOfArray (matrix, running_n * running_n),

    *T          = new double [rnsq],
    *Tcopy      = new double [rnsq],
    *Tbest      = new double [rnsq],

    *wbest      = new double [running_n - 1],
    *zbest      = new double [rnsq],

    *w          = NULL,//new double [running_n - 1],
    *z          = NULL;//new double [rnsq];

  // remove one column/row at a time to get smaller and smaller minor

  while (running_n > 1) {

    rnsq = (running_n - 1) * (running_n - 1);

    best_val = 0.;
    best_idx = -1;

    for (int k=0; k < running_n; ++k) {

      /// construct two matrices T and Tcopy that cut k-th row and k-th column

      for (int i=0, ii=0; i<running_n; i++) {

	if (i==k) continue;

	for (int j=0, jj=0; j<running_n; j++) {

	  if (j==k) continue;

	  int
	    idx1 = (running_n - 1) * ii + jj,
	    idx2 = (running_n - 1) * jj + ii;

	  double val2 = matrixCopy [running_n*i + j];

	  T     [idx1] = 
	  T     [idx2] = 
	  Tcopy [idx1] = 
	  Tcopy [idx2] = val2;

	  ++jj;
	}

	++ii;
      }

      int card_ev;

      (*evdec_num)++;

      //--------------------------------------------------------------------------------------------------------------------------------
      dsyevx_interface (running_n - 1, T, card_ev, w, z, EV_TOL, -COIN_DBL_MAX, 0., 1, (running_n - 1 == min_nz) ? (running_n - 1) : 1);
      //--------------------------------------------------------------------------------------------------------------------------------

      double val = w [0]; // minimum eigenvalue

      if (val < best_val) {

	best_val = val;
	best_idx = k;

	std::memcpy (Tbest, Tcopy, rnsq            * sizeof (double));
	std::memcpy (zbest, z,     rnsq            * sizeof (double));
	std::memcpy (wbest, w,     (running_n - 1) * sizeof (double));

	card_ev_best = card_ev;
      }

      delete [] w;
      delete [] z;

      w = z = NULL;
    }

    // For this value of k, now we have in Tbest the best minor of size
    // running_n

    if (best_idx >= 0) {

      if (del_idx == NULL) {
	del_idx = new bool[n];
	CoinFillN (del_idx, n, false);
      }

      int cnt_idx_orig = 0;
      int cnt_idx_minor = 0;

      while (cnt_idx_minor < running_n) {
	if (del_idx [cnt_idx_orig] == false) {
	  if (cnt_idx_minor == best_idx) {
	    del_idx [cnt_idx_orig] = true;
	    break;
	  }

	  cnt_idx_minor++;
	}

	cnt_idx_orig++;
      }

      if (running_n - 1 == min_nz) {

	for (int i=0; i < card_ev_best && wbest [i] < 0; i++) {

	  CoinFillN (sparse_v_mat [i], n, 0.);

	  double *curr_ev = zbest + i * (running_n - 1);

	  for (int idx_orig = 0, idx_minor = 0; idx_orig < n; ++idx_orig)

	    if (!(del_idx [idx_orig]))
	      sparse_v_mat [i] [idx_orig] = curr_ev [idx_minor++];

	  ++ *card_v_mat;
	}

	break; // done when reached min_nz dimension
      }
    }

    CoinCopyN (Tbest, (n-1) * (n-1), matrixCopy);

    --running_n;

  } // end while

  delete [] del_idx;

  delete [] z;
  delete [] w;

  delete [] T;
  delete [] Tcopy;
  delete [] matrixCopy;

  delete [] Tbest;
  delete [] zbest;
  delete [] wbest;

  delete [] matrix;
} // sparsify2 ()


// Adds SDP cuts using negative eigenvectors where small (smaller than
// 1 / (10 sqrt n)) components are fixed to zero

/************************************************************************/
void CouenneSdpCuts::additionalSDPcuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					CouenneExprMatrix *minor, 
					int n, 
					const double *A, 
					const double *vector,
					int **indA) const {
  int
    *indices = new int [n],
    cnt = 0;

  double threshold = 1 / (10 * sqrt ((double) n));

  for (int i=0; i < n; i++)
    indices [i] = ((fabs (vector [i]) > threshold) ? (cnt++) : -1);

  double *subA = new double [cnt*cnt];

  for (register int i=0, k=0; i<n; i++)

    if (indices [i] >= 0) {
      
      for (register int j=0, k2 = 0; j<n; j++)

	if (indices [j] >= 0) {
	  subA [cnt * k  + k2] = 
	  subA [cnt * k2 + k ] = A [n*i + j];
	  ++k2;
	}

      ++k;
    }

  double *w = NULL, *z = NULL;
  int m;

  //printf ("calling dsyevx ADDITIONAL SDP CUTS\n");

  //////////////////////////////////////////////////////
  dsyevx_interface (cnt, subA, m, w, z, EV_TOL, -COIN_DBL_MAX, onlyNegEV_ ? 0. : COIN_DBL_MAX, 1, cnt);
  //////////////////////////////////////////////////////

  double
    *v    = new double [n],
    *newv = new double [n];

  for (int k=0; k<m; k++) {

    if (onlyNegEV_ && (w [k] >= 0.))
      break;

    double *zbase = z + k * cnt;

    for (int j=0; j<cnt; j++)
      v [j] = *zbase++;

    for(int j=0; j<n; j++)
      newv [j] = (indices [j] >= 0) ? v [indices [j]] : 0.;

    genSDPcut (si, cs, minor, newv, newv, indA);
  }

  delete [] v;
  delete [] newv;

  delete [] w;
  delete [] z;

  delete [] subA;
  delete [] indices;
} // additionalSDPcuts


/************************************************************************/
void CouenneSdpCuts::update_sparsify_structures (const int n, double *v,
						 double* margin, double* A, double *lhs, 
						 const int *zeroed, int evidx, bool decompose, 
						 int *evdec_num) const {

  int minor_n = n;

  if (zeroed != NULL) {
    for (int i=0; i<n; i++)
      if (zeroed [i] == 0)
	--minor_n;
  }

  if (decompose && (minor_n > 2)) {

    /*
      if (minor_n < n) {
      add_v_cut(n, loc_selected, loc_card_selected, locv, 
      init_card_selected, &has_init_vect,
      selected, &card_selected, &new_selected, 
      trace_bin, trace_bin_size,
      sparse_v_mat, card_v_mat);
      }
    */

    (*evdec_num)++;

    double *minor_A = new double [n*n];
    double *minor_w = new double [n];
    double *minor_z = new double [n*n];

    // prepare active submatrix (induced by zeroed vector)

    int
      ii = 0,
      jj = 0;

    for (int i=0;i<n;i++) {
      if (zeroed[i] == 0)
	continue;

      jj = 0;

      for (int j=0;j<n;j++) {
	if (zeroed [j] == 0)
	  continue;
	minor_A [minor_n*ii + jj] = A [n * i + j];
	jj++;
      }

      ii++;
    }

    int m;

    //printf ("calling dsyevx UPDATE_SPARSIFY_STRUCTURES\n");

    //----------------------------------------------------------------------------------------
    dsyevx_interface (minor_n, minor_A, m, minor_w, minor_z, EV_TOL, -COIN_DBL_MAX, 0., 1, 1);
    //----------------------------------------------------------------------------------------

    // update v (reindex the evidx-th eigenvector entries)
    ii = 0;
    for (int i=0;i<n;i++) {
      v[i] = 0;
      if (zeroed[i] == 0)
	continue;
      v[i] = minor_z[ii];
      ii++;
    }

    delete [] minor_A;
    delete [] minor_w;
    delete [] minor_z;
  }

  for   (int i=0; i<n; ++i)
    for (int j=0; j<n; ++j) {
      A [i*n + j] = //*= v[i] * v[j];
      A [j*n + i] *= v[i] * v[j];
      if ((zeroed != NULL) && (zeroed [j] == 0))
	A [i*n + j] = A [j*n + i] = 0;
    }

  *lhs = 0.;

  for (int i=0; i<n; i++) {

    margin[i] = 0;

    for(int j=0; j<n; j++)
      margin[i] += A [i*n + j];

    *lhs += margin[i];
  }
}


/************************************************************************/
void CouenneSdpCuts::zero_comp (const int ind_i, 
				const double delta,
				const int n, 
				const int *selected,
				int *loc_selected, 
				int *ploc_card_selected, 
				int *ploc_card_new_selected, 
				double *ploc_lhs, 
				double *locmargin, 
				double *locmat, 
				double *locv, 
				const int evidx, 
				bool wise, 
				int *evdec_num, 
				double *recomp_gap, 
				double *threshold) const {

  //double  curr_lhs = (*ploc_lhs);
  static int zerocount;
  bool local_wise = false;
  if (wise && (*ploc_lhs - delta > *threshold)) {
    (*threshold) = (*ploc_lhs)-delta + (*recomp_gap);
    local_wise = true;
  }

  zerocount++;

  loc_selected[ind_i] = 0;
  (*ploc_card_selected)--;

  if (selected [ind_i] != 1)
    (*ploc_card_new_selected)--;

  (*ploc_lhs) -= delta;

  update_sparsify_structures (n,locv, locmargin, locmat, ploc_lhs, loc_selected, evidx, local_wise, evdec_num);

}


/************************************************************************/
void CouenneSdpCuts::zero_unified     (enum zero_type type,
				       const int n, 
				       const int *order,
				       const int *selected,
				       const int min_card_new_selected,
				       const double min_delta, 
				       const int start_point, 
				       const int curr_i, 
				       int *loc_selected, 
				       int *ploc_card_selected, 
				       int *ploc_card_new_selected, 
				       double *ploc_lhs, 
				       double *locmargin, 
				       double *locmat, 
				       int *pnchanged, 
				       double *locv, 
				       const int evidx, 
				       bool wise, 
				       double *recomp_gap, 
				       double *threshold,
				       int *evdec_num) const {

  int curr_ind = curr_i;

  *pnchanged = 0;
  for (int i=0; i<n; i++) {
    
    if (++curr_ind == n)
      curr_ind = 0;
    
    int ind_i = order[curr_ind];

    if (((type ==  POS_DELTA || type == VALID_DELTA) && 
	 ((((selected [ind_i] == 0) && (min_card_new_selected >= *ploc_card_new_selected)) || 
	   (curr_ind == start_point) || 
	   (loc_selected[ind_i] == 0)))) 
	||
	((type == SELECTED) && 
	 ((selected[ind_i] == 0) || (loc_selected[ind_i] == 0))))
      continue;
    
    double delta = 2 * locmargin [ind_i] - locmat [ind_i * n + ind_i];
    if (((type == VALID_DELTA || type == SELECTED) && (*ploc_lhs - delta < min_delta)) ||
	((type == POS_DELTA)                       && (delta > 0))) {

      zero_comp (ind_i, delta, n, selected, loc_selected, 
		 ploc_card_selected, ploc_card_new_selected, 
		 ploc_lhs, locmargin, locmat, locv, evidx, wise, evdec_num, recomp_gap, threshold);
      (*pnchanged)++;
    }
  }
}


/************************************************************************/
void CouenneSdpCuts::add_v_cut(const int n,
			       const int *loc_selected, 
			       const int loc_card_selected,
			       const double *locv,
			       const int init_card_selected, int *has_init_vect,
			       int *selected, int *pcard_selected,
			       int *pnew_selected,
			       double **sparse_v_mat,
			       int *pcard_v_mat) const {

  *pnew_selected = 0;

  for (int i=0; i<n; i++) {
    if (loc_selected[i]) {
      sparse_v_mat [*pcard_v_mat] [i] = locv [i];
      if(selected[i] == 0) {
	selected[i] = 1;
	(*pcard_selected)++;
	(*pnew_selected)++;
      }
    }
    else
      sparse_v_mat [*pcard_v_mat][i] = 0;
  }

#ifdef NORMALIZE_SPARSE_CUTS
  //normalization (setting vector norm to 1)
  double curr_norm = 0.0;
  for (int i=0;i<n;i++) {
    curr_norm += fabs(sparse_v_mat[*pcard_v_mat][i]);
  }
  for (int i=0;i<n;i++) {
    if (sparse_v_mat[*pcard_v_mat][i] != 0.0)
      sparse_v_mat[*pcard_v_mat][i] /= curr_norm;
  }
#endif

  if (loc_card_selected + init_card_selected == n) {
    if (*has_init_vect == 1) return;
    else
       (*has_init_vect) = 1;
  }

  (*pcard_v_mat)++;
}


/************************************************************************/
void CouenneSdpCuts::sparsify (bool use_new_sparsify,
			       const int evidx, const double eigen_val, 
			       const double *v, const int n,
			       const double *A, double **sparse_v_mat,
			       int *card_v_mat, int *evdec_num) const {

  int nchanged = 0,
    min_number_new_per_cut = 1,
    loc_card_new_selected  = 0,
    card_selected          = 0,
    loc_card_selected      = 0,

    *selected     = new int [n],
    *loc_selected = new int [n], 
    *order        = new int [n];
	
  double
    min_delta,
    is_zero = 1 / (10 * sqrt ((double) n)),
    lhs     = 0., 
    loc_lhs = 0.,

    *margin    = new double  [n],
    *locv      = new double  [n],
    *locv_orig = new double  [n],
    *locmargin = new double  [n],
    *locmat    = new double  [n*n],
    *mat       = CoinCopyOfArray (A, n*n);

  *card_v_mat = 0;

  for (int i=0; i<n; i++) {

    order [i] = i;

    // zero small components in v
    if (fabs (v[i]) < is_zero) {

      locv_orig [i] = 0;
      selected  [i] = -1; // -1: ind will be set to 0 in loc_selected
      card_selected++;

    } else {
      selected  [i] = 0;
      locv_orig [i] = v[i];
    }
  }

  /// Knuth shuffle for creating a random order
  for (int i=0; i<n; ++i) {

    int 
      newpos = i + (int) floor (((double) (n - i) - 1.e-3) * drand48 ()),
      tmp    = order [newpos];

    order [newpos] = order [i];
    order [i] = tmp;
  }

  // printf ("matrix:\n=================================\n");
  // for (int   i=0; i<n; ++i) {
  //   for (int j=0; j<n; ++j)
  //     printf ("%g ", mat [i*n+j]);
  //   printf ("\n");
  // }
  // printf ("=================================\n");

  update_sparsify_structures (n, locv_orig, margin, mat, &lhs, NULL, evidx, false, evdec_num);

  int
    init_card_selected = card_selected, // to recognize if cut from original
    has_init_vect = 0,                  // vector is generated
    start_point = -1;                   // order [start_point]: index that should not be removed

  min_delta = lhs * (use_new_sparsify ? SPARSIFY_NEW_DELTA : SPARSIFY_OLD_DELTA); // do not weaken the cut too much

  while (card_selected < n) {

    for (int i=0; i<n; i++)

      if (selected [order [i]] == 0) {
	start_point = i;
	break;
      }
    
    loc_card_selected = n;
    loc_card_new_selected = n;
    loc_lhs = lhs;

    double
      recomp_gap = fabs (lhs * WISE_SPARSIFY_GAP),
      threshold  = lhs + recomp_gap;

    // restore locv (might have been changed by WISE_SPARSIFY during the generation of the last sparse cut)
    CoinCopyN (locv_orig, n, locv);
    CoinCopyN (mat, n*n, locmat);

    // printf ("LOCMAT :\n:::::::::::::::::::::::::::::::::\n");
    // for (int   i=0; i<n; ++i) {
    //   for (int j=0; j<n; ++j)
    // 	printf ("%g ", locmat [i*n+j]);
    //   printf ("\n");
    // }
    // printf (":::::::::::::::::::::::::::::::::\n");

    for(int i=0; i<n; i++) {
      if(selected[i] == -1) {
	loc_selected[i] = 0;
	loc_card_selected--;
	loc_card_new_selected--;
      } else {
	loc_selected[i] = 1;
				
	if (selected[i] == 1)
	  loc_card_new_selected--;
      }

      locmargin[i] = margin[i];
    }

    if (loc_lhs >= min_delta) { // CASE 1 ////////////////////

      // use vector as is

      card_selected = n;

      if (onlyNegEV_) {

	int new_selected = 0;
			
	add_v_cut (n, loc_selected, loc_card_selected, locv, 
		   init_card_selected, &has_init_vect,
		   selected, &card_selected, &new_selected, 
		   sparse_v_mat, card_v_mat);
      }

    } else { // CASE 2 ///////////////////////////////////////

      int changed = 1;
			
      while (changed) {

	int curr_i = start_point;

	changed = 0;

	int curr_nchanged = -1; 

	while (curr_nchanged) {

	  zero_unified (SELECTED,
			n, order, selected, min_number_new_per_cut,
			min_delta, start_point,
			curr_i, loc_selected, 
			&loc_card_selected, &loc_card_new_selected, 
			&loc_lhs, locmargin, locmat, 
			&curr_nchanged,locv,evidx, use_new_sparsify, 
			&recomp_gap,
			&threshold,
			evdec_num);
					
	  if (curr_nchanged)
	    nchanged += curr_nchanged;
	}

	curr_nchanged = -1;
	
	while (curr_nchanged) {

	  zero_unified (POS_DELTA,
			n, order, selected, min_number_new_per_cut,
			min_delta, // unused
			start_point, start_point, loc_selected, 
			&loc_card_selected, &loc_card_new_selected, 
			&loc_lhs, locmargin, locmat, 
			&curr_nchanged,locv,evidx,use_new_sparsify, &recomp_gap,&threshold,
			evdec_num);

	  if (curr_nchanged) {
	    nchanged += curr_nchanged;
	    changed = 1;
	  }
	} /* while(pos_nchanged != 0) */

	if (changed)
	  continue;

	curr_i = start_point;

	curr_nchanged = -1;

	if (curr_nchanged) {

	  zero_unified (VALID_DELTA,
			n, order, selected, min_number_new_per_cut,
			min_delta, start_point,
			curr_i, loc_selected, 
			&loc_card_selected, &loc_card_new_selected, 
			&loc_lhs, locmargin, locmat, 
			&curr_nchanged,locv,evidx, use_new_sparsify, &recomp_gap,&threshold,
			evdec_num);

	  if (curr_nchanged) {
	    nchanged += curr_nchanged;
	    changed = 1;
	  }
	}
      } /* while(changed) */

      if ((loc_card_selected < n * (use_new_sparsify ? SPARSIFY_NEW_NZ_THRESHOLD : SPARSIFY_OLD_NZ_THRESHOLD)) || (*card_v_mat == 0)) {

	int new_selected = 0;

	add_v_cut (n, loc_selected, loc_card_selected, locv, 
		   init_card_selected, &has_init_vect,
		   selected, &card_selected, &new_selected, 
		   sparse_v_mat, card_v_mat);
      } else {
	selected [order [start_point]] = 1;
	card_selected++;
      }
    }
  } /* while (card_selected < n) */

  delete[] order;
	
  delete [] mat;
  delete [] locmat;

  delete[] locv;
  delete[] locv_orig;
  delete[] margin;
  delete[] locmargin;
	
  delete[] selected;
  delete[] loc_selected;
}
