/* $Id$
 *
 * Name:    disjunctive_cuts.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <disjunctive_cuts.hpp>

#include <CoinPackedVector.hpp>
#include <CglCutGenerator.hpp>
#include <OsiSolverInterface.hpp>
#include <AsxBmqsTest.hpp>
#include <AsxBmqsUtil.hpp>
#include <AsxUtil.hpp>
#include <AsxBpTest.hpp>
#include <AsxBpFracData.hpp>
#include <AsxUtil.hpp>
#include <OsiClpSolverInterface.hpp>
#include <AsxBpDisjunction.hpp>
#include <AsxBpCGLP.hpp>
#include <cplex.h>
#include <dsyevx_wrapper.hpp>
#include <OsiXxxSolverInterface.hpp>
#include <sdpcuts.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>


void disjunctiveCutGen(const OsiSolverInterface &si, OsiCuts &cs, const double *sol, int n, Tracer *tracer) {
	Timer disjcuts_timer;
	disjcuts_timer.start();
	int origcuts = cs.sizeCuts();

	//prepare matrix X-xxT for eigendecomposition
	double *A = new double[n*n];
	double *z = new double[n*n];
	double *w = new double[n];
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++) {
			A[(i*n) + j] = sol[indexQ(i,j,n)] - (sol[i] * sol[j]);
			A[(j*n) + i] = sol[indexQ(i,j,n)] - (sol[i] * sol[j]);
		}
	}

	double *coeff = new double[n];
	int *ind = new int[n];
	// create indices Xind
	int **Xind = new int*[n];
	for (int i=0;i<n;i++) {
		Xind[i] = new int[n];
	}
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++) {
			Xind[i][j] = indexQ(i,j,n);
			Xind[j][i] = indexQ(i,j,n);
		}
	}

	int m;
	dsyevx_wrapper_only_positive(n,A,m,w,z,tracer);

#ifdef TRACE_DISJUNCTIVE_CUTS
	si.writeMps("outerappr.mps");
	
	printf("X-xxT positive ev: ");
	for (int i=0;i<m;i++) {
	printf("%.5f ",w[i]);
	}
	printf("\n");
#endif

	// Cplex environment for AsxBp
	CPXENVptr asxenv;
	util_open_cpxenv(asxenv);

	// TODO:probably is more efficient if we copy the current solution and the current objective function and we use the same solver iterface to generate min and max, and then we restore current solution and obj functions for warm restart

	// clone the current outer-approximation
	OsiSolverInterface *sifull = si.clone(true);
	//reset current objective function
	for (int i=0;i<sifull->getNumCols();i++) {
		sifull->setObjCoeff(i,0.0);
	}


	// for each positive eigenvector we try to generate a disjunctive cut
	for (int k=0;k<m;k++) {
if (w[k] < 10e-8)
continue;
		double *curr_ev = z + k *n;
		for (int i=0;i<n;i++) {
			sifull->setObjCoeff(i,curr_ev[i]);
		}
		// compute lb
		sifull->setObjSense(1.0); // 1.0 sets min, -1.0 sets max
		sifull->resolve();
		solver_status(sifull);
		double lb = sifull->getObjValue();
		//compute ub
		sifull->setObjSense(-1.0); // 1.0 sets min, -1.0 sets max
		sifull->resolve();
		solver_status(sifull);
		double ub = sifull->getObjValue();
		
		
		// create CoinPackedVector
		int cnt = 0;
		for(int i=0;i<n;i++) {
			if(curr_ev[i] != 0.0) {
				coeff[cnt] = curr_ev[i];
				ind[cnt++] = i;
			}
		}
		CoinPackedVector *curr_ev_cpv = new CoinPackedVector(cnt, ind,coeff, false);

#ifdef TRACE_DISJUNCTIVE_CUTS
		printf("ev[%d]=%.5f   lb=%.5f   ub=%.5f\n",k,w[k],lb,ub);
#endif

		//call to the disjunctive cut generator function
		AsxBpDisjunction* dsj=
			AsxBmqs_compute_disjunction(curr_ev_cpv,lb,ub,Xind,sol,si.getNumCols());

#ifdef TRACE_DISJUNCTIVE_CUTS
		dsj->fprint("dsj.txt");
#endif

		// generate disjunctive cut
		AsxBpFracData* fdata=AsxBpFracData::instance(&si,sol);
		AsxBpCGLP* cglp=AsxBpCGLP::instance(fdata, fdata, dsj, asxenv);

		cglp->generate_cuts(&cs, false);
		
		delete curr_ev_cpv;
		delete cglp;
		delete fdata;
		delete dsj;
	}

	delete sifull;

	delete [] coeff;
	delete [] ind;  

	for (int i=0;i<n;i++)
		delete [] Xind[i];
	delete [] Xind;

	delete [] A;
	delete [] z;
	delete [] w;

	tracer->setDisjunctiveCutsTime(disjcuts_timer.time());
	tracer->setDisjunctiveCutsTotalCuts(cs.sizeCuts() - origcuts);
}



