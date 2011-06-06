/* $Id$
 *
 * Name:    OsiXxxSolverInterface.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef OSIXXXSOLVERINTERFACE_HPP
#define OSIXXXSOLVERINTERFACE_HPP

#ifdef COIN_HAS_CPX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"

class OsiXxxSolverInterface: public OsiCpxSolverInterface {
public:
	
	void myResolve() {
		int status;

		CPXsetintparam( getEnvironmentPtr(), CPX_PARAM_SIMDISPLAY, 0 );
		if ( (status = CPXdualopt( getEnvironmentPtr(), getLpPtr (FREECACHED_RESULTS))))
			printf("ERROR WITH myResolve returncode=%d\n",status);
	}

	double myGetObjVal() {
		int status;
		double objval=0.0;
		if ((status = CPXgetobjval( getEnvironmentPtr(), getLpPtr (FREECACHED_RESULTS), &objval )))
			printf("ERROR WITH myGetObjVal returncode=%d\n",status);
		return objval;
	}

	void myWriteMps(const char *name) {
	 	int status;
		if ((status = CPXwriteprob (getEnvironmentPtr(), getLpPtr (FREECACHED_RESULTS), name, NULL)))
			printf("ERROR WITH myWriteMps returncode=%d\n",status);
		
	}
	
	// NOTE: only modifies cplex problem data NOT OsiSolverInterface cache !
	void XxxModifyCoefficient(int row,int col,double value) {
		if (CPXchgcoef( getEnvironmentPtr (), getLpPtr (KEEPCACHED_RESULTS), row, col, value)) 
			printf("ERROR WITH XxxModifyCoefficient\n");
	}

	void XxxInitialSolveBaropt() {
		CPXCENVptr env = getEnvironmentPtr();
		CPXLPptr lp = getLpPtr(KEEPCACHED_RESULTS);	
		CPXhybbaropt(env,lp,CPX_ALG_NONE);
	}
	
	void XxxResolveBaropt() {
		CPXCENVptr env = getEnvironmentPtr();
		CPXLPptr lp = getLpPtr(KEEPCACHED_RESULTS);
		CPXhybbaropt(env,lp,CPX_ALG_NONE);
	}
};


#else

#include "OsiClpSolverInterface.hpp"
class OsiXxxSolverInterface: public OsiClpSolverInterface {
public:
	void XxxModifyCoefficient(int row,int col,double value) {
		modifyCoefficient (row,col,value,false);
	}
};

#endif

#endif

