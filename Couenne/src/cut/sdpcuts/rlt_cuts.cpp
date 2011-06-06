/* $Id$
 *
 * Name:    rlt_cuts.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <string.h>
#include <OsiCuts.hpp>
#include <OsiSolverInterface.hpp>
#include <populate.hpp>
#include <rlt_cuts.hpp>
#include <tracer.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

#define RLT_CUTS_TOL 0.00000001

/// TODO: missing calls to Tracer

void rltCutsGen(const double *sol, int n, OsiCuts &cs, double *lb, double *ub, int m, Tracer *tracer) {
	int idx;
	double lhs = 0.0;
	//McKormick cuts 
	//for square terms
	for(int i=0;i<n;i++) {
		idx = indexQ(i,i,n);

		if (lb[i] != 0.0){
			lhs = -1.0*sol[idx] + 2.0 * lb[i] * sol[i];
			if (lhs - RLT_CUTS_TOL > lb[i]*lb[i])
				createCut(cs,lb[i]*lb[i], -1.0,idx,-1.0,i,2.0 * lb[i]); 
		}
			// -Xii <= - 2 li xi + li^2
		// else bound -Xii <= 0 already present
		if (ub[i] != 0.0) {
			lhs = -1.0*sol[idx] + 2.0 * ub[i] * sol[i];
			if (lhs - RLT_CUTS_TOL > ub[i]*ub[i])	
				createCut(cs,ub[i]*ub[i], -1.0,idx,-1.0,i,2.0 * ub[i]); 
		}
			// -Xii <= - 2 ui xi + ui^2
		// else bound -Xii <= 0 already present
		
		lhs = 1.0*sol[idx] + (- lb[i]-ub[i]) * sol[i];
		if (lhs - RLT_CUTS_TOL > -lb[i]*ub[i])	
			createCut(cs,-lb[i]*ub[i],-1.0,idx, 1.0,i, - lb[i]-ub[i] );
			//  Xii <= (li+ui) xi - li ui
	}
	//for i!=j
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++) {
			idx = indexQ(i,j,n);
			lhs = -1.0*sol[idx] + lb[j] * sol[i] + lb[i] * sol[j];
			if (lhs - RLT_CUTS_TOL > lb[i]*lb[j])
				createCut(cs,lb[i]*lb[j],-1,idx,-1.0,i,lb[j],j,lb[i]); 
				// -Xij <= - lj xi - li xj + li lj

			lhs = -1.0*sol[idx] + ub[j] * sol[i] + ub[i] * sol[j];
			if (lhs - RLT_CUTS_TOL > ub[i]*ub[j])
				createCut(cs,ub[i]*ub[j],-1,idx,-1.0,i,ub[j],j,ub[i]);
				// -Xij <= - uj xi - ui xj + ui uj

			lhs = 1.0*sol[idx] - ub[j] * sol[i] - lb[i] * sol[j];
			if (lhs - RLT_CUTS_TOL > -lb[i]*ub[j])
				createCut(cs,-lb[i]*ub[j],-1,idx,1.0,i,-ub[j],j,-lb[i]);
				//  Xij <=   uj xi + li xj - li uj

			lhs = 1.0*sol[idx] - lb[j] * sol[i] - ub[i] * sol[j];
			if (lhs - RLT_CUTS_TOL > -ub[i]*lb[j])
				createCut(cs,-ub[i]*lb[j],-1,idx,1.0,i,-lb[j],j,-ub[i]);
				//  Xij <=   lj xi + ui xj - ui lj
		}
}



