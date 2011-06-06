/* $Id$
 *
 * Name:    populate.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <string.h>
#include <OsiCuts.hpp>
#include <OsiSolverInterface.hpp>
#include <CoinMpsIO.hpp>
#include <populate.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))


int getFirstIdx(const char *colName);
int getSecondIdx(const char *colName);
double _mul(double a, double b, double infinity);


int populateProblem (
	const char *filename,
	int *nptr,
	int *tptr,
	int *consptr,
	double **bptr,
	double **cptr,
	double ***Qptr,
	double *constantptr,
	double ***origmatptr,
	double **origrhsptr,
	char **origsenseptr,
	double **xlbptr,
	double **xubptr,
	double **ylbptr,
	double **yubptr,
	OsiSolverInterface *si ) {
	//this function assumes that the MPS file uses the following
	//convention for the variables' names:
	//  x_i represents an original variable of the problem.
	//  x_i_j represents a bilinear term of the form x_i * x_j in the 
	//        original formulation

	double _constant = 0.0;
	char tempstring[1024];
	FILE *mpsfile = NULL;
	mpsfile = fopen(filename,"r");
 	if (mpsfile != NULL)  {
		fscanf(mpsfile,"%s",tempstring); //NAME
		fscanf(mpsfile,"%s",tempstring); //@objconst=<value>
		if (char *tempstring2 = strstr(tempstring,"@objconst=")) {
			_constant= - atof(tempstring2+10);
		}
		fclose(mpsfile);
	} else {
		printf("File %s does not exist !\n",filename);
		return 1;
	}
	

	*constantptr = _constant;

	CoinMpsIO cmpsio;
	cmpsio.messageHandler()->setLogLevel(0);
	cmpsio.readMps(filename,"mps");

	const CoinPackedMatrix *problemMatrixByCol = cmpsio.getMatrixByCol();

	const double *rowlb,*rowub,*objcoeff,*collb,*colub;

	rowlb = cmpsio.getRowLower();
	rowub = cmpsio.getRowUpper();
	objcoeff = cmpsio.getObjCoefficients();
	collb = cmpsio.getColLower();
	colub = cmpsio.getColUpper();

	
	//find out what's the size of the original problem
	int origVars=0;		
	for (int i=0;i<cmpsio.getNumCols();i++) {
		if (getSecondIdx(cmpsio.columnName(i)) == -1)
			origVars++; //if variable's name is of the form x_i (original variable)
	}

	//determine which original variables appear in quadratic terms
	bool *inBilinear;
	inBilinear = new bool[origVars];
	for (int i=0;i<origVars;i++)
		inBilinear[i] = false;

	int origBilinearVars = 0;
	int elimVars = 0;
	for(int i=0;i<cmpsio.getNumCols();i++) {
		if (getSecondIdx(cmpsio.columnName(i)) >= 0) {
			int firstIdx = getFirstIdx(cmpsio.columnName(i));
			int secondIdx = getSecondIdx(cmpsio.columnName(i));
			if ((problemMatrixByCol->getVector(i).getNumElements() > 0) || 
			    (objcoeff[i] != 0.0)) {
				//column is non empty OR objective function for x_i_j != 0
				inBilinear[firstIdx] = true;
				inBilinear[secondIdx] = true;
				origBilinearVars++;
			}
			else
				elimVars++;
		}
	}

	int n = 0;
	int t = 0; //are those variables that do not appear in any bilinear term
	for(int i=0;i<origVars;i++){
		if (inBilinear[i])
			n++;
		else
			t++;
	}

	int N = n*(n+3)/2;
	int totVars = N + t;
	assert(origVars = (elimVars + n + t +origBilinearVars));

	
	#ifdef TRACE
	printf("populate:: MPS file : vars=%d cons=%d\n",cmpsio.getNumCols(),cmpsio.getNumRows());
	if (elimVars)
	printf("populate:: eliminated %d bilinear vars (0.0 obj coefficient AND empty column)\n",elimVars);
	printf("populate:: original problem : n=%d  t=%d  bilinearVars=%d  total=%d\n",
	n,t,origBilinearVars, n+t+origBilinearVars);
	printf("populate:: sdp formulation : n=%d  N=n*(n+3)/2=%d  t=%d  totVars=(N+t)=%d\n",n,N,t,totVars);
	#endif

	//first we do a remapping of the indices (lookupIdx) of the new x_i and y_i
	int *lookupIdx = new int[n+t];
	for(int i=0;i<n+t;i++)
		lookupIdx[i] = -1;
	int xcnt = 0;
	int ycnt = 0;
	for(int i=0;i<cmpsio.getNumCols();i++) {
		int firstIdx = getFirstIdx(cmpsio.columnName(i));
		int secondIdx = getSecondIdx(cmpsio.columnName(i));
		if (secondIdx == -1){
			if (inBilinear[firstIdx]) {
				lookupIdx[firstIdx] = xcnt;
				xcnt++;
			}
			else { //not in any bilinear term, it will be a y_i variable
				lookupIdx[firstIdx] = ycnt;
				ycnt++;
			}
		}
	}
	#ifdef CHECK
	for(int i=0;i<n+t;i++)
		if(lookupIdx[i] == -1)
			printf("populate::WARNING: variable x_%d not found in the original problem\n",i);
	#endif


	double *lb,*ub;
	lb = new double[totVars];
	ub = new double[totVars];
	for(int i=0;i<totVars;i++) {
		lb[i] = - si->getInfinity();
		ub[i] = si->getInfinity();
	}

	// add all n + n(n+1)/2 + t variables
	for (int i=totVars; i--;)
		si -> addCol (0, NULL, NULL,  - si->getInfinity(),  si->getInfinity(), 0);


	int *indices;
	indices = new int[N+t];
	//set OBJ_FUNCTION_MULTIPLIER to 1.0 to keep a minimization pb, -1.0 for maximization
	#define OBJ_FUNCTION_MULTIPLIER -1.0
	for(int i=0;i<cmpsio.getNumCols();i++) {
		int firstIdx = getFirstIdx(cmpsio.columnName(i));
		int secondIdx = getSecondIdx(cmpsio.columnName(i));
		if (secondIdx >= 0){
			if ((inBilinear[firstIdx]) && (inBilinear[secondIdx])) {
				indices[i] = indexQ(lookupIdx[firstIdx],lookupIdx[secondIdx],n);
			}
			else
				indices[i] = -1;
		}
		else {
			if (inBilinear[firstIdx]) {
				indices[i] = lookupIdx[i];
			}
			else {
				indices[i] = N + lookupIdx[i];
			}
		}
		if (indices[i] >= 0) {
			lb[indices[i]] = collb[i];
			ub[indices[i]] = colub[i];
			si->setColLower(indices[i],collb[i]);
			si->setColUpper(indices[i],colub[i]);
			si->setObjCoeff(indices[i], OBJ_FUNCTION_MULTIPLIER * objcoeff[i]);
		}
	}



	//setting lb and ub for extended variables x_i_j
	double min;
	double max;
	for(int i=0;i<n;i++)
		for(int j=i;j<n;j++) {
			min = + si->getInfinity();
			max = - si->getInfinity();
			double tempval;
			tempval = _mul(lb[i],lb[j],si->getInfinity());
			if (tempval < min)
				min = tempval;
			tempval = _mul(lb[i],ub[j],si->getInfinity());
			if (tempval < min)
				min = tempval;
			tempval = _mul(ub[i],lb[j],si->getInfinity());
			if (tempval < min)
				min = tempval;
			tempval = _mul(ub[i],ub[j],si->getInfinity());
			if (tempval < min)
				min = tempval;
			
			tempval = _mul(lb[i],lb[j],si->getInfinity());
			if (tempval > max)
				max = tempval;
			tempval = _mul(lb[i],ub[j],si->getInfinity());
			if (tempval > max)
				max = tempval;
			tempval = _mul(ub[i],lb[j],si->getInfinity());
			if (tempval > max)
				max = tempval;
			tempval = _mul(ub[i],ub[j],si->getInfinity());
			if (tempval > max)
				max = tempval;

			//sets bounds for X_i_j :
			// [ min {l_i*l_j, l_i*u_j, u_i*l_j, u_i*u_j},
			//   max {l_i*l_j, l_i*u_j, u_i*l_j, u_i*u_j} ]
			//for X_i_i we will have the bounds
			// [ max {0 , min {l_i*l_i, l_i*u_i, u_i*l_i, u_i*u_i} },
			//   max {l_i*l_i, l_i*u_i, u_i*l_i, u_i*u_i} ]
#ifdef ZERO_LB_DIAGONAL
			if (i==j) {
				if (min < 0)
					min = 0;
			}

#endif
			si->setColLower(indexQ(i,j,n),min);
			si->setColUpper(indexQ(i,j,n),max);
		}

	OsiCuts cs;


	//original linear cuts
	const CoinPackedMatrix *problemMatrixByRow = cmpsio.getMatrixByRow();
	for(int i=0;i<cmpsio.getNumRows();i++) {
		CoinPackedVector currRow = (CoinPackedVector) problemMatrixByRow->getVector(i);
		int *currRowIdx = (int*) currRow.getIndices();
		double *currRowElem = (double*) currRow.getElements();
		OsiRowCut *cut   = new OsiRowCut;
		
		double *newRowElem = new double[currRow.getNumElements()];
		int *newRowIdx = new int[currRow.getNumElements()];

		for (int j=0;j<currRow.getNumElements();j++) {
			assert(indices[currRowIdx[j]] != -1);
			newRowIdx[j]=indices[currRowIdx[j]];
			newRowElem[j]=currRowElem[j];
		}
		cut->setRow (currRow.getNumElements(), newRowIdx, newRowElem);
		cut->setLb(rowlb[i]);
		cut->setUb(rowub[i]);
		delete [] newRowElem;
		delete [] newRowIdx;

		cs.insert(cut);
	}



	int idx;
#ifdef RLT_CUTS
#ifndef RLT_SEPARATION
	//McKormick cuts 
	//for square terms
	for(int i=0;i<n;i++) {
		idx = indexQ(i,i,n);
//			if (origterm[ind]==true) { //create cuts only for terms x_ij in the original formulation
			if (lb[i] != 0.0)
				createCut(cs,lb[i]*lb[i], -1,idx,-1.0,i,2 * lb[i]); 
				// -Xii <= - 2 li xi + li^2
			// else bound -Xii <= 0 already present
			
			if (ub[i] != 0.0)
				createCut(cs,ub[i]*ub[i], -1,idx,-1.0,i,2 * ub[i]); 
				// -Xii <= - 2 ui xi + ui^2
			// else bound -Xii <= 0 already present
			
			createCut(cs,-lb[i]*ub[i],-1,idx, 1.0,i, - lb[i]-ub[i] );
			//  Xii <= (li+ui) xi - li ui
//			}
	}
	//for i!=j
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++) {
			idx = indexQ(i,j,n);
			//creates cuts only for terms x_ij in the original formulation
				createCut(cs,lb[i]*lb[j],-1,idx,-1.0,i,lb[j],j,lb[i]); 
				// -Xij <= - lj xi - li xj + li lj
				createCut(cs,ub[i]*ub[j],-1,idx,-1.0,i,ub[j],j,ub[i]);
				// -Xij <= - uj xi - ui xj + ui uj
				createCut(cs,-lb[i]*ub[j],-1,idx,1.0,i,-ub[j],j,-lb[i]);
				//  Xij <=   uj xi + li xj - li uj
				createCut(cs,-ub[i]*lb[j],-1,idx,1.0,i,-lb[j],j,-ub[i]);
				//  Xij <=   lj xi + ui xj - ui lj
		}
#endif // #RLT_SEPARATION
#endif //#ifdef RLT_CUTS


	// add constraints
	si -> applyCuts (cs);

	si -> setObjSense(-1); //maximization (the mps are minimization problems, but before we took -1.0 * obj)

	#ifdef TRACE_RLT_MPS
	//write mps formulation

	char **varNames = new char*[N+t];
	for(int j=0;j<N+t;j++)
		varNames[j] = new char[20];

	for (int i=0; i<n; i++) {
		sprintf(varNames[i],"x_%d",i);
	}
	for (int i=0; i<n; i++) {
		sprintf(varNames[indexQ (i,i,n)],"x_%d_%d",i,i);
		for (int j=i+1; j<n; j++) {
			sprintf(varNames[indexQ (i,j,n)],"x_%d_%d",i,j);
		}
	}
	for(int i=0; i<t; i++)
		sprintf(varNames[N+i],"y_%d",i);

	char **conNames = new char*[si->getNumRows()];
	for (int j=0;j<si->getNumRows();j++)
		conNames[j] = new char[20];

	for(int j=0;j<cmpsio.getNumRows();j++)
		sprintf(conNames[j],"orig_%d",j);
	for(int j=cmpsio.getNumRows();j<si->getNumRows();j++)
		sprintf(conNames[j],"rlt_%d",j-cmpsio.getNumRows());
		
	si->writeMpsNative("rlt.mps",
		const_cast <const char **> (conNames), 
		const_cast <const char **> (varNames),
		0,1,1.0);

	for(int j=0;j<N+t;j++)
		delete [] varNames[j];
	delete [] varNames;
	for(int j=0;j<si->getNumRows();j++)
		delete [] conNames[j];
	delete [] conNames;
	#endif

	delete [] lookupIdx;
	delete [] inBilinear;

	delete [] lb;
	delete [] ub;

	delete [] indices;



	// save problem data into Q , origmat, b, c, xlb, xub, ylb, yub
	//objective function: max b x + Q X + c y
	//set nptr, tptr, Qptr, bptr,cptr, origmat,xlb...
	double *b;
	double *c;
	double **Q;
	double **origmat;
	double *origrhs;
	char *origsense;
	double *xlb;
	double *xub;
	double *ylb;
	double *yub;
	int cons;
	cons = cmpsio.getNumRows();
	*consptr = cons;

	*nptr = n;
	*tptr = t;

	const double *tmp_objcoeff = si->getObjCoefficients();
	const double *tmp_lb = si->getColLower();
	const double *tmp_ub = si->getColUpper();
	const double *tmp_rhs = si->getRightHandSide();
	const char   *tmp_sense = si->getRowSense();
	
	xlb = new double[n+t];
	*xlbptr = xlb;
	xub = new double[n+t];
	*xubptr = xub;
	ylb = new double[n+t];
	*ylbptr = ylb;
	yub = new double[n+t];
	*yubptr = yub;

	b = new double[n];
	*bptr = b;
	for (int i=0; i<n; i++) {
		b[i]=tmp_objcoeff[i];
		xlb[i] = tmp_lb[i];
		xub[i] = tmp_ub[i];
	}

	Q = new double*[n];
	for (int i=0; i<n; i++)
		Q[i] = new double[n];
	*Qptr = Q;
	for (int i=0; i<n; i++)
		for (int j=i; j<n; j++) {
			Q[i][j] = tmp_objcoeff[indexQ(i,j,n)] / 2;
			Q[j][i] = tmp_objcoeff[indexQ(i,j,n)] / 2;
			//should we multiply by two the entry Q[i][i] ?
			//check heuristics... we always add Q[i][j] and Q[j][i]
			//if (i==j)
			//	Q[i][i] *= 2.0;
		}


	c = new double[t];
 	*cptr = c;
	for (int i=0;i<t;i++) {
		c[i]=tmp_objcoeff[i+N];
		ylb[i] = tmp_lb[i+N];
		yub[i] = tmp_ub[i+N];
	}


	origmat = new double*[cmpsio.getNumRows()];
	for (int i=0;i<cmpsio.getNumRows();i++)
		origmat[i] = new double[N + t];
	for (int i=0;i<cmpsio.getNumRows();i++)
		for(int j=0;j<N+t;j++)
			origmat[i][j] = 0.0;
	*origmatptr = origmat;

	origrhs = new double[cmpsio.getNumRows()];
	origsense = new char[cmpsio.getNumRows()];
	for (int i=0;i<cmpsio.getNumRows();i++) {
		origrhs[i]= tmp_rhs[i];
		origsense[i] = tmp_sense[i];
	}
	*origrhsptr = origrhs;
	*origsenseptr = origsense;

	const CoinPackedMatrix* tmp_matbyrow = si->getMatrixByRow ();
	for(int i=0;i<cons;i++) {
		CoinPackedVector tmp_currRow = (CoinPackedVector) tmp_matbyrow->getVector(i);
		int *tmpRowIdx = (int*) tmp_currRow.getIndices();
		double *tmpRowElem = (double*) tmp_currRow.getElements();

		for(int j=0;j<tmp_currRow.getNumElements();j++) {
			origmat[i][tmpRowIdx[j]] = tmpRowElem[j];
		}
	}

	return 0;
}


/***********************************************************************/
// creates a cut. return 1 if cut was inserted into cs, 0 otherwise
int createCut (OsiCuts &cs,           // set of cuts
	       double rhs,            // rhs 
	       int sign,              // -1: $\le$, +1: $\ge$, 0: =
	                              // indices, coefficients 
	                              // (index == -1 means don't care)
	       int i1, double c1,     // first 
	       int i2, double c2,     // second
	       int i3, double c3,     // third
	       int i4, double c4,     // forth
	       bool is_global) {      // is this cut global (true) or local (false)?

	int nterms = 0;

	if ((i1 >= 0) && (c1 != 0.0)) nterms++; else c1 =0;
	if ((i2 >= 0) && (c2 != 0.0)) nterms++; else c2 =0;
	if ((i3 >= 0) && (c3 != 0.0)) nterms++; else c3 =0;
	if ((i4 >= 0) && (c4 != 0.0)) nterms++; else c4 =0;
	
	if (!nterms)
		return 0; // nonsense cut
	int tmpcnt = 0;
	
	double    *coeff = new double [nterms]; 
	int       *index = new int    [nterms];
	OsiRowCut *cut   = new OsiRowCut;
	
	if ((i1 >= 0) && (c1 != 0.0))
	{
		coeff[tmpcnt] = c1;
		index[tmpcnt] = i1;
		tmpcnt++;
	}
	if ((i2 >= 0) && (c2 != 0.0))
	{
		coeff[tmpcnt] = c2;
		index[tmpcnt] = i2;
		tmpcnt++;
	}
	if ((i3 >= 0) && (c3 != 0.0))
	{
		coeff[tmpcnt] = c3;
		index[tmpcnt] = i3;
		tmpcnt++;
	}
	if ((i4 >= 0) && (c4 != 0.0))
	{
		coeff[tmpcnt] = c4;
		index[tmpcnt] = i4;
		tmpcnt++;
	}

	if (sign <= 0) cut -> setUb (rhs);
	if (sign >= 0) cut -> setLb (rhs);
	
	cut -> setRow (nterms, index, coeff);
	
	delete [] coeff;
	delete [] index;
	
	cut -> setGloballyValid (is_global); // global?

//TODO: check the precision here
	cs.insertIfNotDuplicate(*cut,CoinAbsFltEq(1.0e-12));
	//cs.insert (cut); //old

	delete cut;
	
	return 1;
}



/***********************************************************************/
int getFirstIdx(const char *colName) {
// if colName = x_i or x_i_j returns i
	char *temp = new char[strlen(colName)+1];
	char *idx1,*idx2;
	int idx1int;
	strcpy(temp,colName);
	idx1 = temp;
	idx1 += 2; // gets rid of 'x_' at the beginning of the name
	idx2 = strchr(idx1,'_');
	if (idx2 != NULL) //variable of the form x_i_j
		*idx2 = '\0';
	idx1int = atoi(idx1);
	delete [] temp;
	return idx1int;
}


/***********************************************************************/
int getSecondIdx(const char *colName) {
// if colName = x_i returns -1
// if colNamr = x_i_j returns j
	char *temp = new char[strlen(colName)+1];
	char *idx1,*idx2;
	int idx2int =-1;
	strcpy(temp,colName);
	idx1 = temp;
	idx1 += 2; // gets rid of 'x_' at the beginning of the name
	idx2 = strchr(idx1,'_');
	if (idx2 != NULL){ //variable of the form x_i_j
		*idx2 = '\0';
		idx2++;
		idx2int = atoi(idx2);
	}
	else
		idx2int = -1;
	delete [] temp;
	return idx2int;
}

/***********************************************************************/
double _mul(double a, double b, double infinity) {
	if ((a==0.0) || (b==0.0))
		return 0.0;

	if ((a >= infinity) && (b > 0))
		return infinity;
	if ((a > 0) && (b >= infinity))
		return infinity;
	if ((a >= infinity) && (b < 0))
		return - infinity;
	if ((a < 0) && (b >= infinity))
		return - infinity;

	if ((a <= -infinity) && (b > 0))
		return - infinity;
	if ((a > 0) && (b <= -infinity))
		return - infinity;
	if ((a <= -infinity) && (b < 0))
		return  infinity;
	if ((a < 0) && (b <= -infinity))
		return - infinity;

	return a*b;
}

