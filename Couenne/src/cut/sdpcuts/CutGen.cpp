/* $Id$
 *
 * Name:    CutGen.cpp
 * Author:  Andrea Qualizza
 * Purpose: Generation of all cust
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "CoinTime.hpp"

#include "CutGen.hpp"
#include "misc_util.hpp"
#include "dsyevx_wrapper.hpp"
#include "orthocut.hpp"
#include "disjunctive_cuts.hpp"
#include "linquad_cuts.hpp"
#include "rlt_cuts.hpp"

static int decomposition_counter;

/************************************************************************/
// sdpcut separator
void CutGen::generateCuts (const OsiSolverInterface &si, OsiCuts &cs, 
			      const CglTreeInfo info) const {
#ifdef TRACE_CUT_TIME
	Timer cut_timer;
	cut_timer.start();
#endif

	int np = n_ + 1;
	int ncuts = cs.sizeRowCuts ();
	const double *sol  = si.getColSolution ();
	double *A = new double[np*np];

	int origsdpcuts_card = 0;
	int duplicate_cuts = 0;

	Timer origsdpcuts_timer;
	Timer sparsify_timer;

	origsdpcuts_timer.start();
	/*
	   A =  /1 x^T\
	        \x  X /
	*/
	A[0] = 1;
	for (int i=0;i<n_;i++) {
		A[np*(i+1)] = sol[i];
		//A[1+i] = sol[i];
		// we need only the upper triangular part, so we can comment out the previous line
	}
	for (int i=0;i<n_;i++) {
		for (int j=i;j<n_;j++) {
			A[(j+1)*np+(i+1)] = sol[indexQ (i, j, n_)];
			//A[(i+1)*np+(j+1)] = sol[indexQ (i, j, n_)]; 
			// we need only the upper triangular part, so we can comment out the previous line
		}
	}

	double *Acopy = new double[np*np];
	for (int i=0;i<np*np;i++)
		Acopy[i] = A[i];

	double *w = NULL, *z = NULL;
	int m;
	dsyevx_wrapper_only_negative (np, A, m, w, z,tracer_);
	//remember that A gets destroyed

	double **work_ev = new double*[m];
	for (int i=0;i<m;i++) {
		work_ev[i] = new double[np];
	
		double *curr_ev ;
		curr_ev = z + (i*np);
#ifdef SCALE_EIGENV
		double scaling_factor = sqrt(np);
		for (int j=0;j<np;j++)
			work_ev[i][j] = curr_ev[j] * scaling_factor;
#else
		for (int j=0;j<np;j++)
			work_ev[i][j] = curr_ev[j];
#endif // SCALE_EIGENV
	}

	for (int i=0;i<m;i++) {
#ifdef INCL_ORIG
#ifdef ONLY_NEG_EIGENV
		if(w [i] > 0)
			break;
#endif

#ifdef ONLY_MOST_NEG
		if(i > 0)
			break;
#endif
		genSDPcut(si,cs,work_ev[i],work_ev[i],removeduplicates_,&duplicate_cuts);
		origsdpcuts_card++;
#endif //INCL_ORIG
	}
	origsdpcuts_timer.pause();

	tracer_->setSDPNumNegativeEV(m);
	if (m >= 1) 
		tracer_->setSDPMostNegativeEV(w[0]);

#if (defined TRACE_ALL)||(defined TRACE_EIGENVALUES)
	printf("Eigenvalues: ");
	for(int i=0;i<m;i++)
		printf("%.5f ",w[i]);
	printf("\n");
#endif

#ifdef RLT_SEPARATION
	rltCutsGen(sol, n_, cs, xlb_, xub_, m, tracer_);
#endif

#ifdef ORTHOCUT
	orthoCutGen(sol, n_, cs, z, w, m, tracer_);
#endif

#ifdef DISJUNCTIVE_CUTS
	disjunctiveCutGen(si, cs, sol, n_, tracer_);
#endif

#ifdef LINQUAD_CUTS
	linQuadCutGen(si.getColSolution(),n_,cs,tracer_);
#endif

//FILE *fsparsifycmp = fopen("sparsifycmp.txt","a");
//compareSparsify(si,n_,m,sol,z,w,fsparsifycmp); 
//fclose(fsparsifycmp);

	int ncuts_beforesparsify = cs.sizeRowCuts ();
	int wise_evdec_num = 0;

	sparsify_timer.start();
#if (defined SPARSIFY) || (defined SPARSIFY2) || (defined WISE_SPARSIFY)
	int card_sparse_v_mat = 0;
	double **sparse_v_mat = new double*[SPARSIFY_MAX_CARD];
	for (int i=0; i<SPARSIFY_MAX_CARD; i++)
		sparse_v_mat[i] = new double[np];
#endif

#ifdef SPARSIFY2
	int min_nz;
	
	min_nz = ceil(np*0.70);
	card_sparse_v_mat = 0;

	sparsify2(n_,sol,sparse_v_mat,&card_sparse_v_mat,min_nz,&wise_evdec_num);
	for(int k=0; k<card_sparse_v_mat; k++) {
		genSDPcut (si, cs, sparse_v_mat[k], sparse_v_mat[k],removeduplicates_,&duplicate_cuts);
	}
#endif // SPARSIFY2

	double *v;
	for (int i=0;i<m;i++) {
		v = work_ev[i];
#ifdef ONLY_NEG_EIGENV
		if(w [i] > 0)
			break;
#endif

#ifdef ONLY_MOST_NEG
		if(i > 0)
			break;
#endif

#if (defined TRACE_ALL) || (defined TRACE_USED_EIGENVECTORS)
		printf("Used Eigenvector idx:%d - (Corresp. to Eigenvalue = %.5f)\n",i,w[i]);
		for(int j=0;j<np;j++)
			printf("%.5f ",v[j]);
		printf("\n");
#endif

#if  (defined SPARSIFY) || (defined WISE_SPARSIFY)
		card_sparse_v_mat = 0;
		double *work = new double[np];
#if (defined WISE_SPARSIFY)
		sparsify_new(i,w[i], v, n_, sol, sparse_v_mat, &card_sparse_v_mat,work,
			true,&wise_evdec_num);
#else
		sparsify(i,w[i], v, n_, sol, sparse_v_mat, &card_sparse_v_mat,work,
			false,&wise_evdec_num);
#endif // (defined WISE_SPARSIFY)
		delete [] work;

		for(int k=0; k<card_sparse_v_mat; k++) {

#ifdef  ADD_RAND_PERC_SPARSIFY_CUTS
			if ( cpp_genalea(seed_) <= ADD_RAND_PERC_SPARSIFY_CUTS)
#endif
			genSDPcut (si, cs, sparse_v_mat[k], sparse_v_mat[k],removeduplicates_,
					&duplicate_cuts);


#ifdef SPARSIFY_MINOR_SDP_CUTS
			additionalSDPcuts(si,cs, np, Acopy, sparse_v_mat[k],
					&duplicate_cuts);
#endif

#ifdef TRACE_ALL
			print_mat_from_vvT(stdout, "generated cut", sparse_v_mat[k], n);
#endif
		}
#endif // (defined SPARSIFY) || (defined WISE_SPARSIFY)

#ifdef TRACE_SPARSIFY_CUT_DEPTH
		sparsify_timer.pause();
		for(int k=0; k<card_sparse_v_mat; k++) {
			double violation = sparse_v_mat[k][0]*sparse_v_mat[k][0];
			for (int p=0;p<n_;p++) {
				violation += 2*sparse_v_mat[k][0]*sparse_v_mat[k][p+1]*sol[p];
				violation += sparse_v_mat[k][p+1]*sparse_v_mat[k][p+1]*sol[indexQ(p,p,n_)];
				for (int q=p+1;q<n_;q++)
					violation += 2*sparse_v_mat[k][p+1]*sparse_v_mat[k][q+1]*sol[indexQ(p,q,n_)];
			}
			printf("violation = %.8f\n",violation);
		}
		sparsify_timer.restore();
#endif
	}


#if (defined SPARSIFY) || (defined SPARSIFY2) || (defined WISE_SPARSIFY)
	for(int i=0;i<SPARSIFY_MAX_CARD;i++)
		delete [] sparse_v_mat[i];
	delete [] sparse_v_mat;
#endif


#if (defined ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS)  ||  (defined TRACE_SPARSIFY_CUT_DEPTH)
	int card_cuts = cs.sizeRowCuts();
	double *cuts_violations = new double[card_cuts - ncuts_beforesparsify];
	int *cuts_indices = new int[card_cuts - ncuts_beforesparsify];
	for (int i=ncuts_beforesparsify;i<card_cuts;i++)
		cuts_indices[i-ncuts_beforesparsify] = i;
	for (int i=ncuts_beforesparsify;i<card_cuts;i++) {
		const OsiRowCut *curr_cut = cs.rowCutPtr(i);
		const double *curr_cut_elem = curr_cut->row().getElements();
		const int *curr_cut_ind = curr_cut->row().getIndices();
		
		const double rhs = curr_cut->rhs();
		double lhs = 0.0;
		for(int j=0;j<curr_cut->row().getNumElements();j++) {
			lhs += curr_cut_elem[j] * sol[curr_cut_ind[j]];
		}
		switch(curr_cut->sense()) {
			case 'L': 
				cuts_violations[i-ncuts_beforesparsify] = lhs - rhs;
				break;
			case 'G': 
				cuts_violations[i-ncuts_beforesparsify] = - lhs + rhs;
				break;
			default: //'E' 
				cuts_violations[i-ncuts_beforesparsify] = max((lhs - rhs),(- lhs + rhs));
		}
	}
	cpp_quicksort_dec(0,card_cuts - ncuts_beforesparsify, cuts_indices, cuts_violations);

//printf("card_cuts=%d ",card_cuts);
//printf(" ncuts_beforesparsify=%d ",ncuts_beforesparsify);
//printf(" card_cuts-ncuts_beforesparsify=%d\n",card_cuts-ncuts_beforesparsify);

#ifdef  ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS
	int card_to_be_removed = (int) floor((card_cuts - ncuts_beforesparsify) *
						(1-ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS));
	if (card_to_be_removed > 0 ) {
		int *cuts_indices_to_be_removed = new int[card_to_be_removed];
		for (int i=0;i<card_to_be_removed;i++) {
			cuts_indices_to_be_removed[i] = cuts_indices[i];
		}
		int *unuseful_vector = new int[card_to_be_removed];
		cpp_quicksortINT_dec(0,card_to_be_removed,unuseful_vector,cuts_indices_to_be_removed);
		//remove cuts from the last one !
		for (int i=0;i<card_to_be_removed;i++)
			cs.eraseRowCut(cuts_indices_to_be_removed[i]);
		delete [] cuts_indices_to_be_removed;
		delete [] unuseful_vector;
	}
#endif // ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS

#ifdef TRACE_SPARSIFY_CUT_DEPTH
	for (int i=0;i<card_cuts- ncuts_beforesparsify;i++)
		printf("sparse_cut[%4d]_violation=%.8f\n",i,cuts_violations[i]);
#endif
	delete [] cuts_violations;
	delete [] cuts_indices;
#endif // (defined TRACE_SPARSIFY_CUT_DEPTH) || (defined ADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS)


	delete [] z;
	delete [] w;
	delete [] A;
	delete [] Acopy;
	for (int i=0;i<m;i++)
		delete [] work_ev[i];
	delete [] work_ev;

#ifdef TRACE_CUT_TIME
	double time = cut_timer.time();
	double timepercut = time / (cs.sizeRowCuts() - ncuts);
	printf("tot_cutgen_time=%.6f time_per_cut=%.6f\n",time,timepercut);
#endif

	tracer_->setSDPCutsTime(origsdpcuts_timer.time());
	tracer_->setSDPCutsTotalCuts(origsdpcuts_card);

#if (defined SPARSIFY) || (defined SPARSIFY2) || (defined WISE_SPARSIFY)
	globaltimer_->pause();
	tracer_->setSparsifyTime(sparsify_timer.time());
	tracer_->setSparsifyTotalCuts(cs.sizeRowCuts() - ncuts_beforesparsify);
	tracer_->setSparsifyDuplicatedCuts(duplicate_cuts);
	tracer_->setSparsifyWiseDecompositions(wise_evdec_num);

	double *violation = new double[cs.sizeRowCuts()-ncuts_beforesparsify];
	int *notused = new int[cs.sizeRowCuts()-ncuts_beforesparsify];
	int *cols_sparsity = new int[np];
	for (int i=0;i<np;i++)
		cols_sparsity[i] = 0;
	int **pair_sparsity = new int*[np];
	for (int i=0;i<np;i++) {
		pair_sparsity[i] = new int[np];
		for (int j=0;j<np;j++)
			pair_sparsity[i][j] = 0;
	}
	double *sparse_vector = new double[np];
	for (int i=0;i<cs.sizeRowCuts()-ncuts_beforesparsify;i++) {
		const OsiRowCut *cut = cs.rowCutPtr(i+ncuts_beforesparsify);
		tracer_->addSparsifyNz(cut->row().getNumElements());
		violation[i] = cut->violated(sol);
		const double *elements = cut->row().getElements();
		const int *indices = cut->row().getIndices();
		int numelements = cut->row().getNumElements();
		for (int k=0;k<np;k++)
			sparse_vector[k] = 0.0;
		sparse_vector[0] = sqrt(cut->rhs());
		for (int k=0;k<n_;k++) {
			int idx = indexQ(k,k,n_);
			for (int j=0;j<numelements;j++)
				if ((indices[j] == idx) && (elements[j] != 0.0))
					sparse_vector[k+1] = sqrt(elements[j]);
		}
		for (int k=0;k<np;k++) {
			if (sparse_vector[k] != 0.0) {
				cols_sparsity[k]++;
				for (int p=k+1;p<np;p++)
					if (sparse_vector[p] != 0.0)
						pair_sparsity[k][p]++;
			}
		}
	}
	cpp_quicksort_dec(0,cs.sizeRowCuts()-ncuts_beforesparsify, notused,violation);
	for (int i=0;i<ceil((cs.sizeRowCuts()-ncuts_beforesparsify)*0.20);i++)
		tracer_->addSparsifyTop20PercCutsViolation(violation[i]);
	for (int j=0;j<np;j++) {
		tracer_->addSparsifySingleColumnSparsity(cols_sparsity[j]);
		for (int p=j+1;p<np;p++)
			tracer_->addSparsifyColumnPairSparsity(pair_sparsity[j][p]);
	}

	delete [] cols_sparsity;
	for (int i=0;i<np;i++)
		delete [] pair_sparsity[i];
	delete [] pair_sparsity;
	delete [] violation;
	delete [] notused;
	delete [] sparse_vector;
	globaltimer_->restore();
#endif

}




double violation_from_v(int n, double *v1, double *v2, const double *sol) {
	double violation = (v1[0]*v2[0]);
	for (int i=0;i<n;i++) {
		violation += v1[0] * v2[i+1] * sol[i];
		violation += v2[0] * v1[i+1] * sol[i];
		violation += v1[i+1] * v2[i+1] * sol[indexQ(i,i,n)];
		for (int j=i+1;j<n;j++) {
			violation += v1[i+1] * v2[j+1] * sol[indexQ(i,j,n)];
			violation += v2[i+1] * v1[j+1] * sol[indexQ(i,j,n)];
		}
	}
	return violation;
}

void CutGen::compareSparsify(const OsiSolverInterface &si,int n, int m, const double *sol, double *z, double *w,FILE *out) const {


	static int iter;
	iter++;
	
	int duplicate_cuts1 = 0;
	int duplicate_cuts2 = 0;
	
	int np = n+1;
	int card_sparse_v_mat1 = 0;
	double **sparse_v_mat1 = new double*[SPARSIFY_MAX_CARD];
	for (int i=0; i<SPARSIFY_MAX_CARD; i++)
		sparse_v_mat1[i] = new double[np];
	int card_sparse_v_mat2 = 0;
	double **sparse_v_mat2 = new double*[SPARSIFY_MAX_CARD];
	for (int i=0; i<SPARSIFY_MAX_CARD; i++)
		sparse_v_mat2[i] = new double[np];

	OsiSolverInterface *sirlt = si.clone(true);
	sirlt->unmarkHotStart();
	Timer timer_rlt;
	timer_rlt.start();	
	sirlt->initialSolve();
	double time_rlt = timer_rlt.time();


	OsiCuts cs_origev;
	for(int i=0;i<m;i++) {
		double *curr_ev ;
		curr_ev = z + (i*np);
		genSDPcut (si, cs_origev, curr_ev, curr_ev,false,&duplicate_cuts1);
	}




//first sparsify
	decomposition_counter = 0;
	int evdec_num1 = 0;
	int *sparsify_cols1 = new int[np];
	for (int i=0;i<np;i++)
		sparsify_cols1[i] = 0;
	int **pair_sparsity1 = new int*[np];
	for (int i=0;i<np;i++) {
		pair_sparsity1[i] = new int[np];
		for (int j=0;j<np;j++)
			pair_sparsity1[i][j] = 0;
	}
	double time1 = 0,starttime1 = 0;

	starttime1 = CoinCpuTime ();


	int gencuts1 = 0;
	OsiCuts cs1;
	double *workevs1 = new double[m];
	for(int i=0;i<m;i++) {
		double *curr_ev ;
		curr_ev = z + (i*np);
		double *work_ev = new double[np];
		sparsify(i,w[i], curr_ev, n, sol, sparse_v_mat1, &card_sparse_v_mat1,work_ev,false,&evdec_num1);
		workevs1[i] = violation_from_v(n_,work_ev,work_ev,sol);
		gencuts1 += card_sparse_v_mat1;
		for (int k=0;k<card_sparse_v_mat1;k++) {
			genSDPcut (si, cs1, sparse_v_mat1[k], sparse_v_mat1[k],true,&duplicate_cuts1);
			for (int j=0;j<np;j++)
				if (sparse_v_mat1[k][j] != 0.0) {
					sparsify_cols1[j]++;
					for (int p=j+1;p<np;p++) {
						if (sparse_v_mat1[k][p] != 0.0)
							pair_sparsity1[j][p]++;
					}
				}
		}
	}

	int decomposition_counter1 = decomposition_counter;

	time1 = CoinCpuTime () - starttime1;
	double time_per_cut1 = time1 / gencuts1;
	int dupcuts1 = gencuts1 - cs1.sizeRowCuts();
	double dupcuts_perc1 = (100.0 * ((double) dupcuts1)) / ((double)gencuts1);

	OsiSolverInterface *si1_lp = si.clone(true);
	double starttime1_lp = CoinCpuTime ();
	si1_lp->applyCuts (cs1);
	si1_lp->resolve();
	double time1_lp = CoinCpuTime () - starttime1_lp;
	double obj1 = si1_lp->getObjValue() + objConst_;

	OsiSolverInterface *si1_lporigcuts = si.clone(true);
	double starttime1_lporigcuts = CoinCpuTime ();
	si1_lporigcuts->applyCuts(cs_origev);
	si1_lporigcuts->applyCuts(cs1);
	si1_lporigcuts->resolve();
	double time1_lporigcuts = CoinCpuTime () - starttime1_lporigcuts;
	double obj1_withorigev = si1_lporigcuts->getObjValue() + objConst_;


	double nz_mean1 = 0.0;
	int nz_min1 = 10000000;
	int nz_max1 = 0;
	for (int i=0;i<cs1.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs1.rowCutPtr(i);
		nz_mean1 += (double) cut->row().getNumElements();
		if (cut->row().getNumElements() > nz_max1)
			nz_max1 = cut->row().getNumElements();
		if (cut->row().getNumElements() < nz_min1)
			nz_min1 = cut->row().getNumElements();
	}
	nz_mean1 /= cs1.sizeRowCuts();
	double nz_std_dev1 = 0.0;
	for (int i=0;i<cs1.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs1.rowCutPtr(i);
		nz_std_dev1 += (((double) cut->row().getNumElements()) - nz_mean1)
				*(((double) cut->row().getNumElements()) - nz_mean1);
	}
	nz_std_dev1 /= cs1.sizeRowCuts();
	nz_std_dev1 = sqrt(nz_std_dev1);

//second sparsify procedure
	decomposition_counter = 0;
	int evdec_num2 = 0;
	int *sparsify_cols2 = new int[np];
	for (int i=0;i<np;i++)
		sparsify_cols2[i] = 0;
	int **pair_sparsity2 = new int*[np];
	for (int i=0;i<np;i++) {
		pair_sparsity2[i] = new int[np];
		for (int j=0;j<np;j++)
			pair_sparsity2[i][j] = 0;
	}
	double time2 = 0,starttime2 = 0;

	starttime2 = CoinCpuTime ();
	int gencuts2 = 0;
	OsiCuts cs2;
	double *workevs2 = new double[m];
	for(int i=0;i<m;i++) {
		double *curr_ev ;
		curr_ev = z + (i*np);
		double *work_ev = new double[np];
		sparsify_new(i,w[i], curr_ev, n, sol, sparse_v_mat2, &card_sparse_v_mat2,work_ev,true,&evdec_num2);
		workevs2[i] = violation_from_v(n_,work_ev,work_ev,sol);
		gencuts2 += card_sparse_v_mat2;
		for (int k=0;k<card_sparse_v_mat2;k++) {
			genSDPcut (si, cs2, sparse_v_mat2[k], sparse_v_mat2[k],true,&duplicate_cuts2);
			for (int j=0;j<np;j++)
				if (sparse_v_mat2[k][j] != 0.0) {
					sparsify_cols2[j]++;
					for (int p=j+1;p<np;p++) {
						if (sparse_v_mat2[k][p] != 0.0)
							pair_sparsity2[j][p]++;
					}
				}
		}
	}

	int decomposition_counter2 = decomposition_counter;

	time2 = CoinCpuTime () - starttime2;
	double time_per_cut2 = time2 / gencuts2;
	int dupcuts2 = gencuts2 - cs2.sizeRowCuts();
	double dupcuts_perc2 = (100.0 * ((double) dupcuts2)) / ((double)gencuts2);
	OsiSolverInterface *si2_lp = si.clone(true);
	double starttime2_lp = CoinCpuTime ();
	si2_lp->applyCuts (cs2);
	si2_lp->resolve();
	double time2_lp = CoinCpuTime () - starttime2_lp;
	double obj2 = si2_lp->getObjValue() + objConst_;
	OsiSolverInterface *si2_lporigcuts = si.clone(true);
	double starttime2_lporigcuts = CoinCpuTime ();
	si2_lporigcuts->applyCuts(cs_origev);
	si2_lporigcuts->applyCuts (cs2);
	si2_lporigcuts->resolve();
	double time2_lporigcuts = CoinCpuTime () - starttime2_lporigcuts;
	double obj2_withorigev = si2_lporigcuts->getObjValue() + objConst_;

	double nz_mean2 = 0.0;
	int nz_min2 = 10000000;
	int nz_max2 = 0;
	for (int i=0;i<cs2.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs2.rowCutPtr(i);
		nz_mean2 += (double) cut->row().getNumElements();
		if (cut->row().getNumElements() > nz_max2)
			nz_max2 = cut->row().getNumElements();
		if (cut->row().getNumElements() < nz_min2)
			nz_min2 = cut->row().getNumElements();
	}
	nz_mean2 /= cs2.sizeRowCuts();
	double nz_std_dev2 = 0.0;
	for (int i=0;i<cs2.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs2.rowCutPtr(i);
		nz_std_dev2 += (((double) cut->row().getNumElements()) - nz_mean2)
				*(((double) cut->row().getNumElements()) - nz_mean2);
	}
	nz_std_dev2 /= cs2.sizeRowCuts();
	nz_std_dev2 = sqrt(nz_std_dev2);


	fprintf(out,"     Cuts Dup Dup%%  TotTime    TimePerCut NZmean NZSD Bound     Bound(+origev)\n");
	fprintf(out,"%2d S1 %3d %3d %2.2f%% %.8f %.8f %5.2f %5.2f %.4f %.4f\n",
		iter,gencuts1,dupcuts1,dupcuts_perc1,time1,time_per_cut1,nz_mean1,nz_std_dev1,obj1,obj1_withorigev);
	fprintf(out,"%2d S2 %3d %3d %2.2f%% %.8f %.8f %5.2f %5.2f %.4f %.4f\n",
		iter,gencuts2,dupcuts2,dupcuts_perc2,time2,time_per_cut2,nz_mean2,nz_std_dev2,obj2,obj2_withorigev);
	fprintf(out,"\n");
	fprintf(out,"      Bound(onlysp)  Time            Bound(+origev) Time       Decomp\n");
	fprintf(out,"%2d S1 %.4f      %.8f      %.4f      %.8f %d\n",
		iter,obj1,time1_lp,obj1_withorigev,time1_lporigcuts,decomposition_counter1);
	fprintf(out,"%2d S2 %.4f      %.8f      %.4f      %.8f %d\n",
		iter,obj2,time2_lp,obj2_withorigev,time2_lporigcuts,decomposition_counter2);

	fprintf(out,"\n");


#ifdef DETAILED_SPARSIFY_COMPARISON
	fprintf(out,"   orig_ev     S1_work     S2_work\n");
	for (int i=0;i<m;i++) {
		fprintf(out,"%2d %2.8f %2.8f %2.8f\n",iter,w[i],workevs1[i],workevs2[i]);
	}
	fprintf(out,"\n");
#endif
	fprintf(out,"sparsity statistics\n");
#ifdef DETAILED_SPARSIFY_COMPARISON
	fprintf(out,"idx: ");
	for (int i=0;i<np;i++) {
		fprintf(out,"%2d ",i);
	}
	fprintf(out,"\n");
	fprintf(out,"S1 : ");
	for (int i=0;i<np;i++) {
		fprintf(out,"%2d ",sparsify_cols1[i]);
	}
	fprintf(out,"\n");
	fprintf(out,"S2 : ");
	for (int i=0;i<np;i++) {
		fprintf(out,"%2d ",sparsify_cols2[i]);
	}
	fprintf(out,"\n");
#endif
	double mean_cols1 = 0.0;
	double mean_cols2 = 0.0;
	int max_cols1 = 0;
	int max_cols2 = 0;
	int min_cols1 = 1000000000;
	int min_cols2 = 1000000000;
	for (int i=0;i<np;i++) {
		mean_cols1 += sparsify_cols1[i];
		mean_cols2 += sparsify_cols2[i];
		if (sparsify_cols1[i] > max_cols1)
			max_cols1 = sparsify_cols1[i];
		if (sparsify_cols2[i] > max_cols2)
			max_cols2 = sparsify_cols2[i];
		if (sparsify_cols1[i] < min_cols1)
			min_cols1 = sparsify_cols1[i];
		if (sparsify_cols2[i] < min_cols2)
			min_cols2 = sparsify_cols2[i];
	}
	mean_cols1 /= np;
	mean_cols2 /= np;
	double stddev_cols1 = 0.0;
	double stddev_cols2 = 0.0;
	for (int i=0;i<np;i++) {
		stddev_cols1 += ((sparsify_cols1[i]-mean_cols1)*(sparsify_cols1[i]-mean_cols1));
		stddev_cols2 += ((sparsify_cols2[i]-mean_cols2)*(sparsify_cols2[i]-mean_cols2));
	}
	stddev_cols1 /= np;
	stddev_cols1 = sqrt(stddev_cols1);
	stddev_cols2 /= np;
	stddev_cols2 = sqrt(stddev_cols2);

	fprintf(out,"single column sparsity:\n");
#ifdef DETAILED_SPARSIFY_COMPARISON
	fprintf(out,"     idx  S1  S2\n");
	for (int i=0;i<np;i++) {
		fprintf(out,"%2d   %3d %3d %3d\n",iter,i,sparsify_cols1[i],sparsify_cols2[i]);
	}
#endif
	fprintf(out,"S1 mean = %.2f  stddev = %.2f  max = %d  min = %d\n"
		,mean_cols1,stddev_cols1,max_cols1,min_cols1);
	fprintf(out,"S2 mean = %.2f  stddev = %.2f  max = %d  min = %d\n"
		,mean_cols2,stddev_cols2,max_cols2,min_cols2);
	//column pair sparsity
	double mean_pair_sparsity1 = 0.0;
	double mean_pair_sparsity2 = 0.0;
	int max_pair_sparsity1 = 0;
	int max_pair_sparsity2 = 0;
	int min_pair_sparsity1 = 1000000000;
	int min_pair_sparsity2 = 1000000000;
	for (int i=0;i<np;i++)
		for (int j=i+1;j<np;j++) {
			mean_pair_sparsity1 += pair_sparsity1[i][j];
			mean_pair_sparsity2 += pair_sparsity2[i][j];
			if (pair_sparsity1[i][j] > max_pair_sparsity1)
				max_pair_sparsity1 = pair_sparsity1[i][j];
			if (pair_sparsity2[i][j] > max_pair_sparsity2)
				max_pair_sparsity2 = pair_sparsity2[i][j];
			if (pair_sparsity1[i][j] < min_pair_sparsity1)
				min_pair_sparsity1 = pair_sparsity1[i][j];
			if (pair_sparsity2[i][j] < min_pair_sparsity2)
				min_pair_sparsity2 = pair_sparsity2[i][j];
		}
	mean_pair_sparsity1 /= (np*(np-1))/2;
	mean_pair_sparsity2 /= (np*(np-1))/2;
	double stddev_pair_sparsity1 = 0.0;
	double stddev_pair_sparsity2 = 0.0;
	for (int i=0;i<np;i++)
		for (int j=i+1;j<np;j++) {
			stddev_pair_sparsity1 += (pair_sparsity1[i][j]-mean_pair_sparsity1) *
						 (pair_sparsity1[i][j]-mean_pair_sparsity1);
			stddev_pair_sparsity2 += (pair_sparsity2[i][j]-mean_pair_sparsity2) *
						 (pair_sparsity2[i][j]-mean_pair_sparsity2);
		}
	stddev_pair_sparsity1 /= (np*(np-1))/2;
	stddev_pair_sparsity2 /= (np*(np-1))/2;
	stddev_pair_sparsity1 = sqrt(stddev_pair_sparsity1);
	stddev_pair_sparsity2 = sqrt(stddev_pair_sparsity2);
	fprintf(out,"column pair sparsity:\n");
#ifdef DETAILED_SPARSIFY_COMPARISON
	fprintf(out,"column pair sparsity:\n");
	fprintf(out,"     idx1 idx2  S1  S2\n");
	for (int i=0;i<np;i++)
		for (int j=i+1;j<np;j++) {
			fprintf(out,"%2d   %3d  %3d  %3d  %3d\n"
				,iter,i,j,pair_sparsity1[i][j],pair_sparsity2[i][j]);
		}
#endif
	fprintf(out,"S1 mean = %.2f  stddev = %.2f  max = %d min = %d\n",mean_pair_sparsity1,stddev_pair_sparsity1,max_pair_sparsity1,min_pair_sparsity1);
	fprintf(out,"S2 mean = %.2f  stddev = %.2f  max = %d  min = %d\n",mean_pair_sparsity2,stddev_pair_sparsity2,max_pair_sparsity2,min_pair_sparsity1);
	fprintf(out,"\n");



	double *violations1 = new double[cs1.sizeRowCuts()];
	int *nz1 = new int[cs1.sizeRowCuts()];
	for (int i=0;i<cs1.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs1.rowCutPtr(i);
		violations1[i] = cut->violated(sol);
		nz1[i] = cut->row().getNumElements();
	}
	cpp_quicksort_dec(0,cs1.sizeRowCuts(), nz1,violations1);

	double *violations2 = new double[cs2.sizeRowCuts()];
	int *nz2 = new int[cs2.sizeRowCuts()];
	for (int i=0;i<cs2.sizeRowCuts();i++) {
		const OsiRowCut *cut = cs2.rowCutPtr(i);
		violations2[i] = cut->violated(sol);
		nz2[i] = cut->row().getNumElements();
	}
	cpp_quicksort_dec(0,cs2.sizeRowCuts(), nz2,violations2);
	int max_card = CoinMax(cs1.sizeRowCuts(),cs2.sizeRowCuts());

	int min_card = CoinMin(cs1.sizeRowCuts(),cs2.sizeRowCuts());
	int top_card = (int) (ceil(min_card * 0.20));
	fprintf(out,"Top 20%% cuts statistics (cardinality=%d)\n",top_card);
	double mean_viol1 = 0.0;
	double mean_viol2 = 0.0;
	for(int i=0;i<top_card;i++) {
		mean_viol1 += violations1[i];
		mean_viol2 += violations2[i];
	}
printf("--->%.18f %.18f %d %.18f %.18f\n",mean_viol1,mean_viol2,top_card,mean_viol1/top_card,mean_viol2/top_card);
	mean_viol1 /= top_card;
	mean_viol2 /= top_card;
	double stddev_viol1 = 0.0;
	double stddev_viol2 = 0.0;
	for(int i=0;i<top_card;i++) {
		stddev_viol1 += (violations1[i] - mean_viol1)*(violations1[i] - mean_viol1);
		stddev_viol2 += (violations2[i] - mean_viol2)*(violations2[i] - mean_viol2);
	}
	stddev_viol1 /= top_card;
	stddev_viol2 /= top_card;
	stddev_viol1 = sqrt(stddev_viol1);
	stddev_viol2 = sqrt(stddev_viol2);
	fprintf(out,"S1 mean violation = %.4f  stddev = %.4f\n",mean_viol1,stddev_viol1);
	fprintf(out,"S2 mean violation = %.4f  stddev = %.4f\n",mean_viol2,stddev_viol2);
#ifdef DETAILED_SPARSIFY_COMPARISON
	fprintf(out,"   S1                  S2               [sparsified cuts only]\n");
	fprintf(out,"   violation    cutNZ        violation    cutNZ  \n");
	for (int i=0;i<max_card;i++) {
		fprintf(out,"%2d ",iter);
		if (i+1<= cs1.sizeRowCuts())
			fprintf(out,"%.8f    %4d",violations1[i],nz1[i]);
		else
			fprintf(out,"    N/A           ");
		fprintf(out,"       ");
		if (i+1<= cs2.sizeRowCuts())
			fprintf(out," %.8f    %4d",violations2[i],nz2[i]);
		else
			fprintf(out,"         N/A      ");
		fprintf(out,"\n");
	}
#endif

	fprintf(out,"\n-------------------------------------------\n\n");



	spartrace.generated_cuts1[iter-1] = gencuts1;
	spartrace.generated_cuts2[iter-1] = gencuts2;
	spartrace.duplicate1[iter-1] = dupcuts1;
	spartrace.duplicate2[iter-1] = dupcuts2;
	spartrace.sparsifytime1[iter-1] = time1;
	spartrace.sparsifytime2[iter-1] = time2;
	spartrace.boundtime1[iter-1] = time1_lporigcuts;
	spartrace.boundtime2[iter-1] = time2_lporigcuts;

	spartrace.nzmean1[iter-1] = nz_mean1;
	spartrace.nzmean2[iter-1] = nz_mean2;
	spartrace.nzmin1[iter-1] = nz_min1;
	spartrace.nzmin2[iter-1] = nz_min2;
	spartrace.nzmax1[iter-1] = nz_max1;
	spartrace.nzmax2[iter-1] = nz_max2;
	spartrace.decomp1[iter-1] = decomposition_counter1;
	spartrace.decomp2[iter-1] = decomposition_counter2;
	spartrace.single_column_sparsity_mean1[iter-1] = mean_cols1;
	spartrace.single_column_sparsity_mean2[iter-1] = mean_cols2;
	spartrace.single_column_sparsity_max1[iter-1] = max_cols1;
	spartrace.single_column_sparsity_max2[iter-1] = max_cols2;
	spartrace.single_column_sparsity_min1[iter-1] = min_cols1;
	spartrace.single_column_sparsity_min2[iter-1] = min_cols2;
	spartrace.column_pair_sparsity_mean1[iter-1] = mean_pair_sparsity1;
	spartrace.column_pair_sparsity_mean2[iter-1] = mean_pair_sparsity2;
	spartrace.column_pair_sparsity_max1[iter-1] = max_pair_sparsity1;
	spartrace.column_pair_sparsity_max2[iter-1] = max_pair_sparsity2;
	spartrace.column_pair_sparsity_min1[iter-1] = min_pair_sparsity1;
	spartrace.column_pair_sparsity_min2[iter-1] = min_pair_sparsity2;
	spartrace.top_cuts_mean_violation1[iter-1] = mean_viol1;
	spartrace.top_cuts_mean_violation2[iter-1] = mean_viol2;
	spartrace.bounds1[iter-1] = obj1_withorigev;
	spartrace.bounds2[iter-1] = obj2_withorigev;
	if (iter == 1)
		spartrace.times1[iter-1] = time_rlt + time1_lporigcuts;
	else
		spartrace.times1[iter-1] = spartrace.times1[iter-2] + time1_lporigcuts;
	if (iter == 1)
		spartrace.times2[iter-1] = time_rlt + time2_lporigcuts;
	else
		spartrace.times2[iter-1] = spartrace.times2[iter-2] + time2_lporigcuts;
	*spartrace.iterations = iter;

	for(int i=0;i<SPARSIFY_MAX_CARD;i++)
		delete [] sparse_v_mat1[i];
	delete [] sparse_v_mat1;
	for(int i=0;i<SPARSIFY_MAX_CARD;i++)
		delete [] sparse_v_mat2[i];
	delete [] sparse_v_mat2;

	delete [] violations1;
	delete [] violations2;
	delete [] nz1;
	delete [] nz2;
	delete [] sparsify_cols1; 
	delete [] sparsify_cols2;
} //compareSparsify()

/************************************************************************/
void CutGen::genSDPcut (const OsiSolverInterface &si,
			   OsiCuts &cs, double *v1, double *v2, bool checkduplicates, int *duplicate_cuts) const {
	
	int nterms = 0;
	int np     = n_+1;
	
	OsiRowCut *cut   = new OsiRowCut;
	double    *coeff = new double [N_];
	int       *ind   = new int    [N_];
	
	// coefficients for X_ij
	for (int i=1; i<np; i++)
		for (int j=i; j<np; j++) {
			double coeff0 = v1 [i] * v2 [j] + v1 [j] * v2 [i];
			if (coeff0 != 0.0) {
				coeff [nterms] = (i==j) ? (0.5 * coeff0) : (coeff0);
				ind   [nterms++] = indexQ (i-1, j-1, n_);
			}
		}

	// coefficients for x_i
	for (int i=1; i<np; i++) {
		double coeff0 = v1 [i] * v2 [0] + v1 [0] * v2 [i];
		if (coeff0 != 0.0) {
			coeff [nterms]   = coeff0;
			ind   [nterms++] = i-1;
		}
	}
	
	cut -> setRow (nterms, ind, coeff);
	cut -> setLb (- *v1 * *v2);



	if(nterms > 0) {
#ifdef RAND_CUT_ADD
		if ( cpp_genalea(seed_) <= RAND_CUT_ADD)
#endif
		{

		if (!(checkduplicates)) {
			globaltimer_->pause();
		}
		CoinAbsFltEq treatAsSame = CoinAbsFltEq(1.0e-8);
		int initial = cs.sizeRowCuts();
		cs.insertIfNotDuplicate (*cut, treatAsSame);
		int final = cs.sizeRowCuts();
		if (initial == final) {
			(*duplicate_cuts) ++;
			// if flag was false, we still add the duplicate cut
			if (!(checkduplicates))
				cs.insert (cut);
		}
		if (!(checkduplicates)) {
			globaltimer_->restore();
		}

		}
	}



	delete cut;
	delete [] ind;
	delete [] coeff;
}


/************************************************************************/
void CutGen::updateSol() {
	heuristics_->run();
}

/************************************************************************/
// constructor
CutGen::CutGen (const int n, 
			const int t,
			const int cons,
			const double objConst,
			const double *b,
			const double *c,
			const double **Q,
			const double **origMat,
			const double *origRhs,
			const char *origSense,
			const double *xlb,
			const double *xub,
			const double *ylb,
			const double *yub,
			OsiSolverInterface *si,
			Timer *globaltimer,
			Tracer *tracer
			):

	n_ (n),
	t_ (t),
	cons_ (cons),
	objConst_ (objConst),
	si_ (si),
	globaltimer_ (globaltimer),
	tracer_ (tracer) {

	N_ = n*(n+3)/2;

	heuristics_ = new
			Heuristics(n,t,cons,objConst,b,c,Q,origMat,origRhs,origSense,xlb,xub,ylb,yub,si,tracer);

	b_ = new double[n];
	c_ = new double[t];
	Q_ = new double*[n];
 	for (int i=0; i<n;i++)
		Q_[i] = new double[n];

 	for (int i=0; i<n;i++) {
		b_ [i] = b [i];
		for (int j=0; j<n;j++) 
			Q_ [i] [j] = Q [i] [j];
	}

	for (int i=0; i<t ;i++) {
		c_[i] = c[i];
	}

	origRhs_ = new double[cons];
	origSense_ = new char[cons];
	origMat_ = new double*[cons];
	for (int i=0; i<cons;i++)
		origMat_[i] = new double[N_+t];
	for (int i=0; i<cons;i++) {
		for(int j=0;j<N_+t;j++)
			origMat_[i][j] = origMat[i][j];
		origRhs_[i] = origRhs[i];
		origSense_[i] = origSense[i];
	}

	xlb_ = new double[n];
	xub_ = new double[n];
	ylb_ = new double[t];
	yub_ = new double[t];

	for (int i=0; i<n;i++) {
		xlb_[i] = xlb[i];
		xub_[i] = xub[i];
	}
	for (int i=0; i<t;i++) {
		ylb_[i] = ylb[i];
		yub_[i] = yub[i];
	}
	
	seed_ = new int[0];
	*seed_ = time(0);

	max_nb_cuts = 100000;

#ifdef SPARSIFY_REMOVE_DUPLICATES
	removeduplicates_ = true;
#else
	removeduplicates_ = false;
#endif



#if 0
	spartrace.generated_cuts1 = new int[EXIT_ON_ITER];
	spartrace.generated_cuts2 = new int[EXIT_ON_ITER];
	spartrace.duplicate1 = new int[EXIT_ON_ITER];
	spartrace.duplicate2 = new int[EXIT_ON_ITER];
	spartrace.sparsifytime1 = new double[EXIT_ON_ITER];
	spartrace.sparsifytime2 = new double[EXIT_ON_ITER];
	spartrace.boundtime1 = new double[EXIT_ON_ITER];
	spartrace.boundtime2 = new double[EXIT_ON_ITER];
	spartrace.nzmean1 = new double[EXIT_ON_ITER];
	spartrace.nzmean2 = new double[EXIT_ON_ITER];
	spartrace.nzmin1 = new double[EXIT_ON_ITER];
	spartrace.nzmin2 = new double[EXIT_ON_ITER];
	spartrace.nzmax1 = new double[EXIT_ON_ITER];
	spartrace.nzmax2 = new double[EXIT_ON_ITER];
	spartrace.decomp1 = new int[EXIT_ON_ITER];
	spartrace.decomp2 = new int[EXIT_ON_ITER];
	spartrace.single_column_sparsity_mean1 = new double[EXIT_ON_ITER];
	spartrace.single_column_sparsity_mean2 = new double[EXIT_ON_ITER];
	spartrace.single_column_sparsity_max1 = new int[EXIT_ON_ITER];
 	spartrace.single_column_sparsity_max2 = new int[EXIT_ON_ITER];
	spartrace.single_column_sparsity_min1 = new int[EXIT_ON_ITER];
	spartrace.single_column_sparsity_min2 = new int[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_mean1 = new double[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_mean2 = new double[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_max1 = new int[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_max2 = new int[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_min1 = new int[EXIT_ON_ITER];
	spartrace.column_pair_sparsity_min2 = new int[EXIT_ON_ITER];
	spartrace.top_cuts_mean_violation1 = new double[EXIT_ON_ITER];
	spartrace.top_cuts_mean_violation2 = new double[EXIT_ON_ITER];
	spartrace.bounds1 = new double[EXIT_ON_ITER];
	spartrace.times1 = new double[EXIT_ON_ITER];
	spartrace.bounds2 = new double[EXIT_ON_ITER];
	spartrace.times2 = new double[EXIT_ON_ITER];
	spartrace.iterations = new int[1];
#endif

}
/************************************************************************/
// destructor
CutGen::~CutGen () {
	delete seed_;

	delete heuristics_;

	delete [] b_;
	delete [] c_;
	for(int i=0;i<n_;i++)
		delete [] Q_[i];
	delete [] Q_;
	for (int i=0;i<cons_;i++)
		delete [] origMat_[i];
	delete [] origMat_;
	delete [] origRhs_;
	delete [] origSense_;
	delete [] xlb_;
	delete [] xub_;
	delete [] ylb_;
	delete [] yub_;


	delete [] spartrace.generated_cuts1;
	delete [] spartrace.generated_cuts2;
	delete [] spartrace.duplicate1;
	delete [] spartrace.duplicate2;
	delete [] spartrace.sparsifytime1;
	delete [] spartrace.sparsifytime2;
	delete [] spartrace.boundtime1;
	delete [] spartrace.boundtime2;
	delete [] spartrace.nzmean1;
	delete [] spartrace.nzmean2;
	delete [] spartrace.nzmin1;
	delete [] spartrace.nzmin2;
	delete [] spartrace.nzmax1;
	delete [] spartrace.nzmax2;
	delete [] spartrace.decomp1;
	delete [] spartrace.decomp2;
	delete [] spartrace.single_column_sparsity_mean1;
	delete [] spartrace.single_column_sparsity_mean2;
	delete [] spartrace.single_column_sparsity_max1;
 	delete [] spartrace.single_column_sparsity_max2;
	delete [] spartrace.single_column_sparsity_min1;
	delete [] spartrace.single_column_sparsity_min2;
	delete [] spartrace.column_pair_sparsity_mean1;
	delete [] spartrace.column_pair_sparsity_mean2;
	delete [] spartrace.column_pair_sparsity_max1;
	delete [] spartrace.column_pair_sparsity_max2;
	delete [] spartrace.column_pair_sparsity_min1;
	delete [] spartrace.column_pair_sparsity_min2;
	delete [] spartrace.top_cuts_mean_violation1;
	delete [] spartrace.top_cuts_mean_violation2;
	delete [] spartrace.bounds1;
	delete [] spartrace.times1;
	delete [] spartrace.bounds2;
	delete [] spartrace.times2;
	delete [] spartrace.iterations;
}

/***********************************************************************/
void CutGen::myremoveBestOneRowCol(double *matrix, int n, int running_n, int min_nz,bool *del_idx, double **sparse_v_mat, int *card_v_mat, int *evdec_num) const {
	
	double best_val=1;
	int best_idx=-1;

	if(running_n==1) 
		return;

	double *matrixCopy = new double[(running_n)*(running_n)];
	for (int i=0;i<(running_n)*(running_n);i++)
		matrixCopy[i] = matrix[i];

	double *T = new double[(running_n-1)*(running_n-1)];
	double *Tcopy = new double[(running_n - 1)*(running_n - 1)];
	double *Tbest = new double[(running_n - 1)*(running_n - 1)];
	double *zbest = new double[(running_n - 1)*(running_n - 1)];
	double *wbest = new double[running_n - 1];
	double *w,*z;
	int card_ev_best;

	for(int k=0;k<running_n;k++) {
		int ii,jj;
		ii=0;
		for(int i=0;i<running_n;i++) {
			if(i==k) continue;
			jj=0;
			for(int j=0;j<running_n;j++) {
				if(j==k) continue;
				T[(running_n-1)*ii+jj]=matrixCopy[running_n*i+j];
				Tcopy[(running_n-1)*ii+jj]=matrixCopy[running_n*i+j];
				jj++;
			}
			ii++;
		}
		int card_ev;
		w=NULL;
		z=NULL;

		(*evdec_num)++;
		if (running_n-1 == min_nz)
			dsyevx_wrapper_only_negative(running_n - 1,T,card_ev,w,z,tracer_);
		else
			dsyevx_wrapper_only_most_neg(running_n - 1,T,card_ev,w,z,tracer_);

		double val=w[0];

		if(val<0 && val<best_val) {
			best_val=val;
			best_idx=k;
			for(int i=0;i<(running_n-1)*(running_n-1);i++) {
				Tbest[i] = Tcopy[i];
				zbest[i] = z[i];
			}
			for(int i=0;i<(running_n-1);i++) {
				wbest[i] = w[i];
			}
			card_ev_best = card_ev;
		}
		delete [] z;
		delete [] w;
	}
	delete [] T;
	delete [] Tcopy;
	delete [] matrixCopy;


	if(best_idx>=0) {
		
		if (del_idx == NULL) {
			del_idx = new bool[n];
			for (int i=0;i<n;i++)
				del_idx[i]=false;
		}
		int cnt_idx_orig = 0;
		int cnt_idx_minor = 0;
		while (cnt_idx_minor < running_n) {
			if (del_idx[cnt_idx_orig] == false) {
				if (cnt_idx_minor == best_idx) {
					del_idx[cnt_idx_orig] = true;
					break;
				}
				else 
					cnt_idx_minor++;
			}
			cnt_idx_orig++;
		}

		if (running_n-1 == min_nz) {
			for(int i=0;i<card_ev_best;i++) {
				if (wbest[i] < 0) {
					double *curr_ev = zbest + (i*(running_n-1));
					
					for(int j=0;j<n;j++)
						sparse_v_mat[i][j]=0.0;
					
					int idx_orig = 0;
					int idx_minor = 0;

					while (idx_orig < n) {
						if (!(del_idx[idx_orig])) {
							sparse_v_mat[i][idx_orig] = curr_ev[idx_minor];
							idx_minor++;
						}
						idx_orig++;
					}

					(*card_v_mat)++;
				}
				else //no more negative eigenvalues
					break;
			}
			delete [] del_idx;
		}
		else {
			myremoveBestOneRowCol(Tbest, n, running_n-1,min_nz,del_idx,sparse_v_mat,card_v_mat,evdec_num);
		}
	}
	
	delete [] Tbest;
	delete [] zbest;
	delete [] wbest;

}// myremoveBestOneRowCol()
/************************************************************************/
void CutGen::sparsify2(const int n,
			 const double *sol, double **sparse_v_mat,
			 int *card_v_mat, int min_nz, int *evdec_num) const {
	
	int np = n+1;
	double *matrix = new double[np*np];

	matrix[0] = 1;
	for (int i=0;i<n;i++)
		matrix[np*(i+1)] = sol[i];
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++)
			matrix[(j+1)*np+(i+1)] = sol[indexQ (i, j, n)];
	}

	myremoveBestOneRowCol(matrix, np, np,min_nz,NULL,sparse_v_mat,card_v_mat,evdec_num);

	delete [] matrix;
}// sparsify2()
/************************************************************************/
void CutGen::additionalSDPcuts(const OsiSolverInterface &si,OsiCuts &cs, int np, const double *A, const double *vector, int *duplicate_cuts) const{

	int *indices;
	indices = new int[np];
	int cnt = 0;
	for(int i=0;i<np;i++) {
		if (vector[i] != 0.0)
			indices[i] = cnt++;
		else
			indices[i] = -1;
	}
	
	double *subA = new double[cnt*cnt];

	for (register int i=0; i<np; i++) {
		if (indices[i] >= 0) {
			for (register int j=0; j<np; j++) {
				if (indices[j] >= 0)
					subA[cnt*indices[j]+indices[i]] = A[np*j+i];
			}
		}
	}

	double *w = NULL, *z = NULL;
	int m;
	dsyevx_wrapper_only_negative (cnt, subA, m, w, z,tracer_);

	double *v = new double[np];
	double *newv = new double[np];


	for (int k=0; k<m; k++) {
	
#ifdef ONLY_NEG_EIGENV
		if(w [k] > 0) {
			break;
		}
#endif

		double *zbase = z + k * cnt;
		for (int j=0; j<cnt; j++) {
			v [j] = *zbase++;
		}

		for(int j=0;j<np;j++) {
			if (indices[j] >= 0)
				newv[j] = v[indices[j]];
			else
				newv[j] = 0;
		}

		genSDPcut (si, cs, newv, newv,removeduplicates_,duplicate_cuts);
	}

	delete [] v;
	delete [] newv;

	delete [] w;
	delete [] z;

	delete [] subA;
	delete [] indices;
} // additionalSDPcuts







/************************************************************************/
void CutGen::update_sparsify_structures(const int np, const double *sol, double *v,double* margin, double** mat, double *lhs, const int *zeroed, int evidx, bool decompose, int *evdec_num) const {

	// copy sol[] in mat[][]
	mat[0][0] = 1;
	for(int i=1; i<np; i++) {
		mat[0][i] = sol[i-1];
		mat[i][0] = sol[i-1];
	}
	for(int i=1; i<np; i++) {
		for(int j=i; j<np; j++) {
			int ind = indexQ(i-1, j-1, np-1);
			mat[i][j] = sol[ind];
			mat[j][i] = sol[ind];
		}
	}

	int minor_n = np;
	if (zeroed != NULL) {
		for(int i=0;i<np;i++)
			if (zeroed[i] == 0)
				minor_n--;
	}

	if ((decompose)  && (minor_n > 2)) {

/*
	if (minor_n < np) {
			add_v_cut(np, loc_selected, loc_card_selected, locv, 
				init_card_selected, &has_init_vect,
				selected, &card_selected, &new_selected, 
				trace_bin, trace_bin_size,
				sparse_v_mat, card_v_mat);
	}
*/
		decomposition_counter++;
		(*evdec_num)++;
		double *minor_A = new double[np*np];
		double *minor_w = new double[np];
		double *minor_z = new double[np*np];

		//prepare active submatrix (induced by zeroed vector)

		int ii = 0;
		int jj = 0;
		for (int i=0;i<np;i++) {
			if (zeroed[i] == 0)
				continue;
			jj = 0;
			for (int j=0;j<np;j++) {
				if (zeroed[j] == 0)
					continue;
				minor_A[(minor_n*ii) + jj] = mat[i][j];
				jj++;
			}
			ii++;
		}
		//eigendecomposition
		int m;
//		dsyevx_wrapper_first_p (minor_n, minor_A, m, minor_w, minor_z,evidx+1,tracer_);
		dsyevx_wrapper_only_most_neg (minor_n, minor_A, m, minor_w, minor_z,tracer_);

		//update v (reindex the evidx-th eigenvector entries)
		ii = 0;
		for (int i=0;i<np;i++) {
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

	for(int i=0; i<np; i++) {
		for(int j=0; j<np; j++) {
			mat[i][j] *= v[i] * v[j];
			if ((zeroed != NULL) && (zeroed[j] == 0)) {
				mat[i][j] = 0;
				mat[j][i] = 0;
			}
		}
	}

	(*lhs) = 0;
	for(int i=0; i<np; i++) {
		margin[i] = 0;
		for(int j=0; j<np; j++) {
			margin[i] += mat[i][j];	
		}
		(*lhs) += margin[i];
	}
}
/************************************************************************/
void CutGen::zero_comp(const int ind_i, const double delta,
			  const int np, const int *selected,
			  int *loc_selected, 
			  int *ploc_card_selected, int *ploc_card_new_selected, 
			  double *ploc_lhs, 
			  double *locmargin, double **locmat, 
			  const double *sol, double *locv, 
			  const int evidx, bool wise, int *evdec_num, double *recomp_gap, double *threshold) const {

double  curr_lhs = (*ploc_lhs);
static int zerocount;
bool local_wise = false;
if ((wise) && ((*ploc_lhs)-delta > (*threshold))) {
	(*threshold) = (*ploc_lhs)-delta + (*recomp_gap);
	local_wise = true;
}


zerocount++;


  loc_selected[ind_i] = 0;
  (*ploc_card_selected)--;
  
  if(selected[ind_i] != 1) {
    (*ploc_card_new_selected)--;
  }
  (*ploc_lhs) -= delta;

  update_sparsify_structures(np,sol,locv,locmargin,locmat,ploc_lhs, loc_selected, evidx, local_wise, evdec_num);

} /* zero_comp */

/************************************************************************/
void CutGen::zero_valid_delta(const int np, const int *order,
				 const int * selected,
				 const int min_card_new_selected,
				 const double min_delta, const int start_point, 
				 const int curr_i, 
				 int *loc_selected, 
				 int *ploc_card_selected, 
				 int *ploc_card_new_selected, 
				 double *ploc_lhs, 
				 double *locmargin, double **locmat, 
				 int *pnchanged, 
				 const double *sol, double *locv, 
				 const int evidx, bool wise,double *recomp_gap, double *threshold,
			         int *pcard_selected,
			         int *pnew_selected,
			         int *trace_bin, const int trace_bin_size,
			         double **sparse_v_mat,
			         int *pcard_v_mat,
			         const int init_card_selected, int *has_init_vect,
			         int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged = 0);
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];
    int skip = 0;

    if((selected[ind_i] == 0) && 
       (min_card_new_selected >= *ploc_card_new_selected)) {
      skip = 1;
    }

    if((skip) || (curr_ind == start_point) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(*ploc_lhs - delta < min_delta) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx, wise, evdec_num , recomp_gap,threshold);
      (*pnchanged)++;

    }
  }
} /* zero_valid_delta */

/************************************************************************/
void CutGen::zero_selected(const int np, const int *order,
			      const int *selected,
			      const int min_card_new_selected,
			      const double min_delta, const int start_point,
			      const int curr_i, 
			      int *loc_selected, int *ploc_card_selected, 
			      int *ploc_card_new_selected, 
			      double *ploc_lhs, 
			      double *locmargin, double **locmat, 
			      int *pnchanged, 
			      const double *sol, double *locv, 
			      const int evidx, bool wise,double *recomp_gap, double *threshold,
			      int *pcard_selected,
			      int *pnew_selected,
			      int *trace_bin, const int trace_bin_size,
			      double **sparse_v_mat,
			      int *pcard_v_mat,
			      const int init_card_selected, int *has_init_vect,
			      int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged = 0);
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];

    if((selected[ind_i] == 0) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(*ploc_lhs - delta < min_delta) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx,wise,evdec_num,recomp_gap,threshold);
      (*pnchanged)++;
    } 
  }
} /* zero_selected */

/************************************************************************/
void CutGen::zero_pos_delta(const int np, const int *order,
			       const int *selected,
			       const int min_card_new_selected,
			       const int start_point, const int curr_i, 
			       int *loc_selected, int *ploc_card_selected, 
			       int *ploc_card_new_selected, 
			       double *ploc_lhs, 
			       double *locmargin, double **locmat, 
			       int *pnchanged, 
			       const double *sol, double *locv, 
			       const int evidx, bool wise, double *recomp_gap, double *threshold,
			       int *pcard_selected,
			       int *pnew_selected,
			       int *trace_bin, const int trace_bin_size,
			       double **sparse_v_mat,
			       int *pcard_v_mat,
			       const int init_card_selected, int *has_init_vect,
			       int *evdec_num) const {

  int curr_ind = curr_i;

  (*pnchanged) = 0;
  for(int i=0; i<np; i++) {
    
    curr_ind++;
    if(curr_ind == np) {
      curr_ind = 0;
    }
    
    int ind_i = order[curr_ind];
    
    int skip = 0;

    if((selected[ind_i] == 0) && 
       (min_card_new_selected >= *ploc_card_new_selected)) {
      skip = 1;
    }

    if((skip) || (curr_ind == start_point) || (loc_selected[ind_i] == 0)) {
      continue;
    }
    
    double delta = 2 * locmargin[ind_i] - locmat[ind_i][ind_i];
    if(delta > 0) {

      zero_comp(ind_i, delta, np, selected, loc_selected, 
		ploc_card_selected, ploc_card_new_selected, 
		ploc_lhs, locmargin, locmat, sol, locv, evidx,wise,evdec_num,recomp_gap,threshold);
      (*pnchanged)++;

    } 
  }
} /* zero_pos_delta */

/************************************************************************/
void CutGen::add_v_cut(const int np,
			  const int *loc_selected, 
			  const int loc_card_selected,
			  const double *locv,
			  const int init_card_selected, int *has_init_vect,
			  int *selected, int *pcard_selected,
			  int *pnew_selected,
			  int *trace_bin, const int trace_bin_size,
			  double **sparse_v_mat,
			  int *pcard_v_mat) const {

  (*pnew_selected) = 0;

  for(int i=0; i<np; i++) {
    if(loc_selected[i]) {
      sparse_v_mat[*pcard_v_mat][i] = locv[i];
      if(selected[i] == 0) {
	selected[i] = 1;
	(*pcard_selected)++;
	(*pnew_selected)++;
      }
    }
    else {
      sparse_v_mat[*pcard_v_mat][i] = 0;
    }
  }


#ifdef NORMALIZE_SPARSE_CUTS
//normalization (setting vector norm to 1)
double curr_norm = 0.0;
for (int i=0;i<np;i++) {
	curr_norm += fabs(sparse_v_mat[*pcard_v_mat][i]);
}
for (int i=0;i<np;i++) {
	if (sparse_v_mat[*pcard_v_mat][i] != 0.0)
		sparse_v_mat[*pcard_v_mat][i] = sparse_v_mat[*pcard_v_mat][i]/curr_norm;
}
#endif

/*
  if(loc_card_selected + init_card_selected == np) {
    if(*has_init_vect == 1) {

#ifdef TRACE_ALL
      printf("SdpCutGen::add_v_cut(): repeat of original cut skipped\n");
#endif

      return;
    }
    else {
      (*has_init_vect) = 1;
    }
  }
*/
	
#ifdef TRACE_ALL
  printf("SdpCutGen::add_v_cut(): loc_card_selected: %d  new_selected: %d\n", 
	 loc_card_selected, *pnew_selected);
#endif
	
  (*pcard_v_mat)++;
  
#ifdef TRACE_ALL
  trace_bin[loc_card_selected / trace_bin_size] += 1;
#endif
  
#ifdef TRACE_ALL
  cpp_printvecDBL("SdpCutGen::add_v_cut(): sparse vector", 
		  sparse_v_mat[(*pcard_v_mat) - 1], np);
#endif
} /* add_v_cut */

/************************************************************************/
void CutGen::sparsify(const int evidx, const double eigen_val, 
			 const double *v, const int n,
			 const double *sol, double **sparse_v_mat,
			 int *card_v_mat, double *work_ev,bool wise,int *evdec_num) const {

	int i, j, np = n+1, nchanged = 0;
	double sq_np = sqrt((double)np);

	double min_delta;
	double is_zero = 1/(10 * sq_np);
	int min_number_new_per_cut = 1;
	
	int *selected = new int[np], card_selected = 0;
	int *loc_selected = new int[np], loc_card_selected = 0;
	int loc_card_new_selected = 0;
	
	double lhs = 0, loc_lhs = 0;
	double *margin = new double[np];
	double *locv = new double[np];
	double *locv_orig = new double[np];
	double *locmargin = new double[np];
	double **mat = new double*[np];
	double **locmat = new double*[np];
	
	int seed = 225535;
	int *order = new int[np];
	double *rand_val = new double[np];

	*card_v_mat = 0;
	
	for (i=0; i<np; i++) {
		selected[i] = 0;
		mat[i] = new double[np];
		locmat[i] = new double[np];
		order[i] = i;
		rand_val[i] = cpp_genalea(&seed);
		
		
		// zero small components in v
		if(fabs(v[i]) < is_zero) {

#ifdef TRACE_ALL
			printf("zero: ind: %d  value: %8.6f\n", i, v[i]);
#endif

			locv_orig[i] = 0;
			selected[i] = -1; // -1: ind will be set to 0 in loc_selected
			card_selected++;
		} else {
			locv_orig[i] = v[i];
		}
	}


	for (int i=0;i<np;i++)
		work_ev[i] = locv_orig[i];

	// get random ordering
	cpp_quicksort_dec(0, np, order, rand_val);
	

	update_sparsify_structures(np,sol,locv_orig,margin,mat,&lhs, NULL, evidx, false,evdec_num);

	int init_card_selected = card_selected; // to recognize if cut from original
						// vector is generated
	int has_init_vect = 0;
	
	min_delta = lhs * SPARSIFY_OLD_DELTA; // do not weaken the cut too much
	int start_point = -1; // order[start_point]: index that should not be removed
	
	int trace_bin_size = 0;
	int *trace_bin = NULL;

#ifdef TRACE_ALL
	trace_bin_size = 5;
	int card_trace_bin = np / trace_bin_size + 1;
	trace_bin = new int[card_trace_bin];
	for(i=0; i<card_trace_bin; i++) {
		trace_bin[i] = 0;
	}
#endif


	while(card_selected < np) {
		for(i=0; i<np; i++) {
			if(selected[order[i]] == 0) {
				start_point = i;
				break;
			}
		}
    
		loc_card_selected = np;
		loc_card_new_selected = np;
		loc_lhs = lhs;
		double recomp_gap = fabs(lhs*WISE_SPARSIFY_GAP);
		double threshold = lhs + recomp_gap;

		// restore locv (might have been changed by WISE_SPARSIFY during the generation of the last sparse cut)
		for(i=0;i<np;i++)
			locv[i] = locv_orig[i];

		for(i=0; i<np; i++) {
			if(selected[i] == -1) {
				loc_selected[i] = 0;
				loc_card_selected--;
				loc_card_new_selected--;
			} else {
				loc_selected[i] = 1;
				
				if(selected[i] == 1) {
					loc_card_new_selected--;
				}
			}
			locmargin[i] = margin[i];
			for(j=0; j<np; j++) {
				locmat[i][j] = mat[i][j];
			}
		}

		if(loc_lhs < min_delta) {

			int changed = 1;
			
			while(changed) {

#ifdef TRACE_ALL
				printf("SdpCutGen::sparsify(): loc_lhs: %8.6f\n", loc_lhs);
				cpp_printmatDBL("locmat", locmat, np, np);
				cpp_printvecDBL("locmargin", locmargin, np);
				cpp_printvecINT("loc_selected", loc_selected, np);
#endif


				int curr_i = start_point;
				
				changed = 0;
				
				int sel_nchanged = -1; 

				while(sel_nchanged != 0) {
					int new_selected = 0;
					zero_selected(np, order, selected, min_number_new_per_cut,
						min_delta, start_point,
						curr_i, loc_selected, 
						&loc_card_selected, &loc_card_new_selected, 
						&loc_lhs, locmargin, locmat, 
						&sel_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
						&card_selected, &new_selected, 
						trace_bin, trace_bin_size,
						sparse_v_mat, card_v_mat,
						init_card_selected, &has_init_vect,evdec_num);
					
					if(sel_nchanged) {
						nchanged += sel_nchanged;
					//	    changed = 1;
					}
				} // while(sel_nchanged != 0) 

				int pos_nchanged = -1;
	
				while(pos_nchanged != 0) {
					int new_selected = 0;
					zero_pos_delta(np, order, selected, min_number_new_per_cut,
							start_point, start_point, loc_selected, 
							&loc_card_selected, &loc_card_new_selected, 
							&loc_lhs, locmargin, locmat, 
							&pos_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
							&card_selected, &new_selected, 
							trace_bin, trace_bin_size,
							sparse_v_mat, card_v_mat,
							init_card_selected, &has_init_vect,evdec_num);
				
					if(pos_nchanged) {
						nchanged += pos_nchanged;
						changed = 1;
					}
				} /* while(pos_nchanged != 0) */

				if(changed) {
					continue;
				}
				

				curr_i = start_point;
				
				int val_nchanged = -1;

				if(val_nchanged) {
					int new_selected = 0;
					zero_valid_delta(np, order, selected, min_number_new_per_cut,
							min_delta, start_point,
							curr_i, loc_selected, 
							&loc_card_selected, &loc_card_new_selected, 
							&loc_lhs, locmargin, locmat, 
							&val_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
							&card_selected, &new_selected, 
							trace_bin, trace_bin_size,
							sparse_v_mat, card_v_mat,
							init_card_selected, &has_init_vect,evdec_num);
					
					if(val_nchanged) {
					nchanged += val_nchanged;
					changed = 1;
					}
				}


			} /* while(changed) */

			if((loc_card_selected < np * SPARSIFY_OLD_NZ_THRESHOLD) || (*card_v_mat == 0)) {
				
				int new_selected = 0;
				
				add_v_cut(np, loc_selected, loc_card_selected, locv, 
					init_card_selected, &has_init_vect,
					selected, &card_selected, &new_selected, 
					trace_bin, trace_bin_size,
					sparse_v_mat, card_v_mat);
			} else {
				selected[order[start_point]] = 1;
				card_selected++;
			}
		} else {
		// loc_lhs >= min_delta  use vector as is

			card_selected = np;

#ifdef TRACE_ALL
			printf("SdpCutGen::sparsify(): lhs: %8.6f  too large. No sparsification\n", lhs);
#endif

#ifndef ONLY_NEG_EIGENV
			int new_selected = 0;
			
			add_v_cut(np, loc_selected, loc_card_selected, locv, 
					init_card_selected, &has_init_vect,
					selected, &card_selected, &new_selected, 
					trace_bin, trace_bin_size,
					sparse_v_mat, card_v_mat);
#endif

		}
	} /* while(card_selected < np) */


#ifdef TRACE_ALL
	printf("SdpCutGen::sparsify(): bin size: %d\n", trace_bin_size);
	cpp_printvecINT("trace_bin", trace_bin, card_trace_bin);
	delete[] trace_bin;
#endif

	delete[] order;
	delete[] rand_val;
	
	for (i=0; i<np; i++) {
		delete [] mat[i];
		delete [] locmat[i];
	}
	delete [] mat;
	delete [] locmat;

	delete[] locv;
	delete[] locv_orig;
	delete[] margin;
	delete[] locmargin;
	
	delete[] selected;
	delete[] loc_selected;
} // sparsify
/************************************************************************/
void CutGen::sparsify_new(const int evidx, const double eigen_val, 
			 const double *v, const int n,
			 const double *sol, double **sparse_v_mat,
			 int *card_v_mat, double *work_ev, bool wise, int *evdec_num) const {

	int i, j, np = n+1, nchanged = 0;
	double sq_np = sqrt((double)np);

	double min_delta;
	double is_zero = 1/(10 * sq_np);
	int min_number_new_per_cut = 1;
	
	int *selected = new int[np], card_selected = 0;
	int *loc_selected = new int[np], loc_card_selected = 0;
	int loc_card_new_selected = 0;
	
	double lhs = 0, loc_lhs = 0;
	double *margin = new double[np];
	double *locv = new double[np];
	double *locv_orig = new double[np];
	double *locmargin = new double[np];
	double **mat = new double*[np];
	double **locmat = new double*[np];
	
	int seed = 225535;
	int *order = new int[np];
	double *rand_val = new double[np];

	*card_v_mat = 0;
	
	for (i=0; i<np; i++) {
		selected[i] = 0;
		mat[i] = new double[np];
		locmat[i] = new double[np];
		order[i] = i;
		rand_val[i] = cpp_genalea(&seed);
		
		
		// zero small components in v
		if(fabs(v[i]) < is_zero) {

#ifdef TRACE_ALL
			printf("zero: ind: %d  value: %8.6f\n", i, v[i]);
#endif

			locv_orig[i] = 0;
			selected[i] = -1; // -1: ind will be set to 0 in loc_selected
			card_selected++;
		} else {
			locv_orig[i] = v[i];
		}
	}


	for (int i=0;i<np;i++)
		work_ev[i] = locv_orig[i];

	// get random ordering
	cpp_quicksort_dec(0, np, order, rand_val);
	
	update_sparsify_structures(np,sol,locv_orig,margin,mat,&lhs, NULL, evidx, false, evdec_num);

	int init_card_selected = card_selected; // to recognize if cut from original
						// vector is generated
	int has_init_vect = 0;
	
min_delta = lhs * SPARSIFY_NEW_DELTA; // do not weaken the cut too much
	int start_point = -1; // order[start_point]: index that should not be removed
	
	int trace_bin_size = 0;
	int *trace_bin = NULL;

#ifdef TRACE_ALL
	trace_bin_size = 5;
	int card_trace_bin = np / trace_bin_size + 1;
	trace_bin = new int[card_trace_bin];
	for(i=0; i<card_trace_bin; i++) {
		trace_bin[i] = 0;
	}
#endif


	while(card_selected < np) {
		for(i=0; i<np; i++) {
			if(selected[order[i]] == 0) {
				start_point = i;
				break;
			}
		}
    
		loc_card_selected = np;
		loc_card_new_selected = np;
		loc_lhs = lhs;
		double recomp_gap = fabs(lhs*WISE_SPARSIFY_GAP);
		double threshold = lhs + recomp_gap;

		// restore locv (might have been changed by WISE_SPARSIFY during the generation of the last sparse cut)
		for(i=0;i<np;i++)
			locv[i] = locv_orig[i];

		for(i=0; i<np; i++) {
			if(selected[i] == -1) {
				loc_selected[i] = 0;
				loc_card_selected--;
				loc_card_new_selected--;
			} else {
				loc_selected[i] = 1;
				
				if(selected[i] == 1) {
					loc_card_new_selected--;
				}
			}
			locmargin[i] = margin[i];
			for(j=0; j<np; j++) {
				locmat[i][j] = mat[i][j];
			}
		}

		if(loc_lhs < min_delta) {

			int changed = 1;
			while(changed) {

#ifdef TRACE_ALL
				printf("SdpCutGen::sparsify(): loc_lhs: %8.6f\n", loc_lhs);
				cpp_printmatDBL("locmat", locmat, np, np);
				cpp_printvecDBL("locmargin", locmargin, np);
				cpp_printvecINT("loc_selected", loc_selected, np);
#endif
	
				int curr_i = start_point;
				
				changed = 0;
				
				int sel_nchanged = -1; 

				while(sel_nchanged != 0) {
					int new_selected = 0;
					zero_selected(np, order, selected, min_number_new_per_cut,
						min_delta, start_point,
						curr_i, loc_selected, 
						&loc_card_selected, &loc_card_new_selected, 
						&loc_lhs, locmargin, locmat, 
						&sel_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
						&card_selected, &new_selected, 
						trace_bin, trace_bin_size,
						sparse_v_mat, card_v_mat,
						init_card_selected, &has_init_vect,evdec_num);
					
					if(sel_nchanged) {
						nchanged += sel_nchanged;
					//	    changed = 1;
					}
				} // while(sel_nchanged != 0) 

				int pos_nchanged = -1;
	
				while(pos_nchanged != 0) {
					int new_selected = 0;
					zero_pos_delta(np, order, selected, min_number_new_per_cut,
							start_point, start_point, loc_selected, 
							&loc_card_selected, &loc_card_new_selected, 
							&loc_lhs, locmargin, locmat, 
							&pos_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
							&card_selected, &new_selected, 
							trace_bin, trace_bin_size,
							sparse_v_mat, card_v_mat,
							init_card_selected, &has_init_vect,evdec_num);
				
					if(pos_nchanged) {
						nchanged += pos_nchanged;
						changed = 1;
					}
				} /* while(pos_nchanged != 0) */
	
				if(changed) {
					continue;
				}

				curr_i = start_point;
				
				int val_nchanged = -1;

				if(val_nchanged) {
					int new_selected = 0;
					zero_valid_delta(np, order, selected, min_number_new_per_cut,
							min_delta, start_point,
							curr_i, loc_selected, 
							&loc_card_selected, &loc_card_new_selected, 
							&loc_lhs, locmargin, locmat, 
							&val_nchanged,sol,locv,evidx,wise,&recomp_gap,&threshold,
							&card_selected, &new_selected, 
							trace_bin, trace_bin_size,
							sparse_v_mat, card_v_mat,
							init_card_selected, &has_init_vect,evdec_num);
					
					if(val_nchanged) {
					nchanged += val_nchanged;
					changed = 1;
					}
				}
			} /* while(changed) */

			if((loc_card_selected < np * SPARSIFY_NEW_NZ_THRESHOLD) || (*card_v_mat == 0)) {
				
				int new_selected = 0;
				
				add_v_cut(np, loc_selected, loc_card_selected, locv, 
					init_card_selected, &has_init_vect,
					selected, &card_selected, &new_selected, 
					trace_bin, trace_bin_size,
					sparse_v_mat, card_v_mat);
			} else {
				selected[order[start_point]] = 1;
				card_selected++;
			}
		} else {
		// loc_lhs >= min_delta  use vector as is
			card_selected = np;

#ifdef TRACE_ALL
			printf("SdpCutGen::sparsify(): lhs: %8.6f  too large. No sparsification\n", lhs);
#endif

#ifndef ONLY_NEG_EIGENV
			int new_selected = 0;
			
			add_v_cut(np, loc_selected, loc_card_selected, locv, 
					init_card_selected, &has_init_vect,
					selected, &card_selected, &new_selected, 
					trace_bin, trace_bin_size,
					sparse_v_mat, card_v_mat);
#endif

		}
	} /* while(card_selected < np) */


#ifdef TRACE_ALL
	printf("SdpCutGen::sparsify(): bin size: %d\n", trace_bin_size);
	cpp_printvecINT("trace_bin", trace_bin, card_trace_bin);
	delete[] trace_bin;
#endif

	delete[] order;
	delete[] rand_val;
	
	for (i=0; i<np; i++) {
		delete [] mat[i];
		delete [] locmat[i];
	}
	delete [] mat;
	delete [] locmat;

	delete[] locv;
	delete[] locv_orig;
	delete[] margin;
	delete[] locmargin;
	
	delete[] selected;
	delete[] loc_selected;
} // sparsify_new


