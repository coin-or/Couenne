/* */
// (C) Copyright International Business Machines Corporation (IBM) 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pietro Belotti, Carnegie Mellon University
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006


#include "BonCouenneInterface.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

/** Default constructor. */
CouenneInterface::CouenneInterface():
  AmplInterface(),
  have_nlp_solution_ (false)
{}

/** Copy constructor. */
CouenneInterface::CouenneInterface(const CouenneInterface &other):
  AmplInterface(other),
  have_nlp_solution_ (false)
{}

/** virtual copy constructor. */
CouenneInterface * CouenneInterface::clone(bool CopyData){
  return new CouenneInterface(*this);
}

/** Destructor. */
CouenneInterface::~CouenneInterface(){
}

#ifdef COUENNEINTERFACE_FROM_ASL
void
CouenneInterface::readAmplNlFile(char **& argv, Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                                 Ipopt::SmartPtr<Ipopt::OptionsList> options,
                                 Ipopt::SmartPtr<Ipopt::Journalist> journalist){
  //  if (!IsValid (app_))
  //createApplication (roptions, options, journalist, "couenne.");
  AmplInterface::readAmplNlFile(argv, roptions, options, journalist);
}
#endif

/** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   *
   * Solve the continuous relaxation and takes first-order
   * outer-approximation constraints at the optimum.  Then put
   * everything in an OsiSolverInterface.
   *
   * The OsiSolverInterface si is empty and has to be populated with
   * the initial linear relaxation.
   */

void
CouenneInterface::extractLinearRelaxation
(OsiSolverInterface &si, CouenneCutGenerator & couenneCg, bool getObj, bool solveNlp) {

  {
    int nlpLogLevel;
    options () -> GetIntegerValue ("nlp_log_level", nlpLogLevel, "couenne.");
    messageHandler () -> setLogLevel (nlpLogLevel);
  }

  CouenneProblem *p = couenneCg.Problem ();
  bool is_feasible = true;

  if (solveNlp) {

    int nvars = p -> nVars();

    if (p -> doFBBT ()) {

      // include the rhs of auxiliary-based constraints into the FBBT
      // (should be useful with Vielma's problems, for example)

      for (int i=0; i < p -> nCons (); i++) {

	// for each constraint
	CouenneConstraint *con = p -> Con (i);

	// (which has an aux as its body)
	int index = con -> Body () -> Index ();

	if ((index >= 0) && (con -> Body () -> Type () == AUX)) {

	  // if there exists violation, add constraint
	  CouNumber
	    l = con -> Lb () -> Value (),	
	    u = con -> Ub () -> Value ();

	  // tighten bounds in Couenne's problem representation
	  p -> Lb (index) = CoinMax (l, p -> Lb (index));
	  p -> Ub (index) = CoinMin (u, p -> Ub (index));
	}
      }

      t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

      for (int i=0; i<nvars; i++) {
	chg_bds [i].setLower(t_chg_bounds::CHANGED);
	chg_bds [i].setUpper(t_chg_bounds::CHANGED);
      }

      if (!(p -> boundTightening (chg_bds, CglTreeInfo (), NULL))) {
	is_feasible = false;
	//*messageHandler() << "Couenne: Warning, tightened NLP is infeasible" << CoinMessageEol;
      }

      delete [] chg_bds;

      const double
	*nlb = getColLower (),
	*nub = getColUpper ();

      for (int i=0; i < p -> nOrigVars () - p -> nDefVars (); i++)
	if (p -> Var (i) -> Multiplicity () > 0) {
	  /*printf ("---- %4d [%g,%g] [%g,%g]\n", i,
		  nlb [i], nub [i],
		  p -> Lb (i), p -> Ub (i));*/

	  double
	    lower = nlb [i],
	    upper = nub [i];

	  if (lower       > upper)       CoinSwap (lower,       upper);
	  if (p -> Lb (i) > p -> Ub (i)) CoinSwap (p -> Lb (i), p -> Ub (i));

	  if (lower < p -> Lb (i) - COUENNE_EPS) setColLower (i, CoinMin (nub[i], p -> Lb (i)));
	  if (upper > p -> Ub (i) + COUENNE_EPS) setColUpper (i, CoinMax (nlb[i], p -> Ub (i)));

	} else {
	  // if not enabled, fix them in the NLP solver
	  setColLower (i, -COIN_DBL_MAX);
	  setColUpper (i,  COIN_DBL_MAX);
	}
    } // ends FBBT part

    if (is_feasible) {
      try {
	options () -> SetNumericValue ("max_cpu_time", CoinMax (0.1, (p -> getMaxCpuTime () - CoinCpuTime ()) / 2));

	initialSolve ();

	if (isDualObjectiveLimitReached() &&
	    (getNumIntegers () == 0))
	  *messageHandler () << "Couenne: Warning, NLP is unbounded" << CoinMessageEol;
      }
      catch (Bonmin::TNLPSolver::UnsolvedError *E) {
	// wrong, if NLP has problems this is not necessarily true...
	//is_feasible = false;
      }
    }

    if (!is_feasible) {
      OsiAuxInfo * auxInfo = si.getAuxiliaryInfo ();
      Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (auxInfo);

      if (babInfo)
	babInfo -> setInfeasibleNode ();
    }

    if (is_feasible && isProvenOptimal ()) {

      CouNumber obj             = getObjValue    ();
      const CouNumber *solution = getColSolution ();

      if (getNumIntegers () > 0) {

	int
	  norig = p -> nOrigVars () - p -> nDefVars (),
	  nvars = p -> nVars ();

	bool fractional = false;

	// if problem is integer, check if any integral variable is
	// fractional. If so, round them and re-optimize

	for (int i=0; i<norig; i++)
	  if ((p -> Var (i) -> Multiplicity () > 0) &&
	      p  -> Var (i) -> isDefinedInteger () &&
	      (!::isInteger (solution [i]))) {
	    fractional = true;
	    break;
	  }

	if (fractional) { // try again if solution found by Ipopt is fractional

	  double
	    *lbSave = new double [norig],
	    *ubSave = new double [norig],

	    *lbCur  = new double [nvars],
	    *ubCur  = new double [nvars],

	    *Y      = new double [nvars];

	  CoinCopyN (getColLower (), norig, lbSave);
	  CoinCopyN (getColUpper (), norig, ubSave);

	  CoinFillN (Y,     nvars, 0.);
	  CoinFillN (lbCur, nvars, -COUENNE_INFINITY);
	  CoinFillN (ubCur, nvars,  COUENNE_INFINITY);

	  CoinCopyN (getColLower (), norig, lbCur);
	  CoinCopyN (getColUpper (), norig, ubCur);

	  if (p -> getIntegerCandidate (solution, Y, lbCur, ubCur) >= 0) {

	    for (int i = getNumCols (); i--;) {

	      if (lbCur [i] > ubCur [i])
		CoinSwap (lbCur [i], ubCur [i]);

	      if      (Y [i] < lbCur [i]) Y [i] = lbCur [i];
	      else if (Y [i] > ubCur [i]) Y [i] = ubCur [i];
	    }

	    for (int i=0; i<norig; i++)
	      if ((p -> Var (i) -> Multiplicity () > 0) &&
		  p  -> Var (i) -> isDefinedInteger ()) {

		setColLower (i, lbCur [i]);
		setColUpper (i, ubCur [i]);
	      }

	    setColSolution (Y); // use initial solution given

	    try {
	      options () -> SetNumericValue ("max_cpu_time", CoinMax (0.1, p -> getMaxCpuTime () - CoinCpuTime ()));

	      resolve (); // solve with integer variables fixed

	      if (isDualObjectiveLimitReached() &&
		  (getNumIntegers () == 0))
		*messageHandler () << "Couenne: Warning, NLP is is unbounded" << CoinMessageEol;
	    }
	    catch (Bonmin::TNLPSolver::UnsolvedError *E) {
	    }

	    //resolve ();

	    obj      = getObjValue ();
	    solution = getColSolution ();

	    // restore previous bounds on integer variables
	    for (int i=0; i<norig; i++)
	      if ((p -> Var (i) -> Multiplicity () > 0) &&
		  p  -> Var (i) -> isDefinedInteger ()) {

		if (lbSave [i] > ubSave [i])
		  CoinSwap (lbSave [i], ubSave [i]);

		setColLower (i, lbSave [i]);
		setColUpper (i, ubSave [i]);
	      }
	  }

	  delete [] Y;
	  delete [] lbSave;
	  delete [] ubSave;
	  delete [] lbCur;
	  delete [] ubCur;
	}
      }

      // re-check optimality in case resolve () was called
      if (isProvenOptimal () &&
	  //	  (obj < p -> getCutOff ()) && // check #1 (before re-computing) -- BUG. What if real object is actually better?

#ifdef FM_CHECKNLP2
	  (p->checkNLP2(solution, 0, false, true, false, p->getFeasTol())) &&
	  (p->getRecordBestSol()->getModSolVal() < p->getCutOff())
#else
	  p -> checkNLP (solution, obj, true) && // true for recomputing obj
	  (obj < p -> getCutOff ())
#endif
	  ) {           // check #2 (real object might be different)

	// tell caller there is an initial solution to be fed to the initHeuristic
	have_nlp_solution_ = true;

	// set cutoff to take advantage of bound tightening

#ifdef FM_CHECKNLP2
	obj = p->getRecordBestSol()->getModSolVal();
#endif

	p -> setCutOff (obj, solution);

	OsiAuxInfo * auxInfo = si.getAuxiliaryInfo ();
	Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (auxInfo);

	if (babInfo) {

#ifdef FM_CHECKNLP2
	  babInfo -> setNlpSolution (p->getRecordBestSol()->modSol,
				     getNumCols(), obj);
#else
	  babInfo -> setNlpSolution (solution, getNumCols (), obj);
#endif
	  babInfo -> setHasNlpSolution (true);
	}

#ifdef FM_TRACE_OPTSOL
#ifdef FM_CHECKNLP2
	p->getRecordBestSol()->update();
#else
	p->getRecordBestSol()->update(solution, getNumCols(),
				      obj, p->getFeasTol());
#endif
#endif

      }
    } else {


    }
  }

  if (!is_feasible) // nothing else to do, problem infeasible
    return;

  int
    numcols     = p -> nOrigVars (), // # original               variables
    numcolsconv = p -> nVars     (); // # original + # auxiliary variables

  const double
    *lb = p -> Lb (), //getColLower (),
    *ub = p -> Ub (); //getColUpper ();

   // add original and auxiliary variables to the new problem
   for (int i=0; i<numcols; i++)
     if (p -> Var (i) -> Multiplicity () > 0) si.addCol (0, NULL,NULL, lb [i],       ub [i],      0);
     else                                     si.addCol (0, NULL,NULL, -COIN_DBL_MAX,COIN_DBL_MAX,0);
   for (int i=numcols; i<numcolsconv; i++)    si.addCol (0, NULL,NULL, -COIN_DBL_MAX,COIN_DBL_MAX,0);

   // get initial relaxation
   OsiCuts cs;
   couenneCg.generateCuts (si, cs);

   // store all (original + auxiliary) bounds in the relaxation
   CouNumber * colLower = new CouNumber [numcolsconv];
   CouNumber * colUpper = new CouNumber [numcolsconv];

   CouNumber *ll = p -> Lb ();
   CouNumber *uu = p -> Ub ();

   // overwrite original bounds, could have improved within generateCuts
   for (int i = numcolsconv; i--;)
     if (p -> Var (i) -> Multiplicity () > 0) {
       colLower [i] = (ll [i] > - COUENNE_INFINITY) ? ll [i] : -COIN_DBL_MAX;
       colUpper [i] = (uu [i] <   COUENNE_INFINITY) ? uu [i] :  COIN_DBL_MAX;
     } else {
       colLower [i] = -COIN_DBL_MAX;
       colUpper [i] =  COIN_DBL_MAX;
     }

   int numrowsconv = cs.sizeRowCuts ();

   // create matrix and other stuff
   CoinBigIndex * start = new CoinBigIndex [numrowsconv + 1];

   int    * length   = new int    [numrowsconv];
   double * rowLower = new double [numrowsconv];
   double * rowUpper = new double [numrowsconv];

   start[0] = 0;
   int nnz = 0;
   /* fill the four arrays. */
   for(int i = 0 ; i < numrowsconv ; i++)
   {
     OsiRowCut * cut = cs.rowCutPtr (i);

     const CoinPackedVector &v = cut->row();
     start[i+1] = start[i] + v.getNumElements();
     nnz += v.getNumElements();
     length[i] = v.getNumElements();

     rowLower[i] = cut->lb();
     rowUpper[i] = cut->ub();
   }

   assert (nnz == start [numrowsconv]);
   /* Now fill the elements arrays. */
   int * ind = new int[start[numrowsconv]];
   double * elem = new double[start[numrowsconv]];
   for(int i = 0 ; i < numrowsconv ; i++) {

     OsiRowCut * cut = cs.rowCutPtr (i);

     const CoinPackedVector &v = cut->row();

     if(v.getNumElements() != length[i])
       std::cout<<"Empty row"<<std::endl;
     //     cut->print();
     CoinCopyN (v.getIndices(),  length[i], ind  + start[i]);
     CoinCopyN (v.getElements(), length[i], elem + start[i]);
   }

   // Ok everything done now create interface
   CoinPackedMatrix A;
   A.assignMatrix(false, numcolsconv, numrowsconv,
                  start[numrowsconv], elem, ind,
                  start, length);
   if(A.getNumCols() != numcolsconv || A.getNumRows() != numrowsconv){
     std::cout<<"Error in row number"<<std::endl;
   }
   assert(A.getNumElements() == nnz);
   // Objective function
   double * obj = new double[numcolsconv];
   CoinFillN(obj,numcolsconv,0.);

   // some instances have no (or null) objective function, check it here
   if (p -> nObjs () > 0)
     p -> fillObjCoeff (obj);

   // Finally, load interface si with the initial LP relaxation
   si.loadProblem (A, colLower, colUpper, obj, rowLower, rowUpper);

   delete [] rowLower;
   delete [] rowUpper;
   delete [] colLower;
   delete [] colUpper;
   delete [] obj;

   for (int i=0; i<numcolsconv; i++)
     if ((p -> Var (i) -> Multiplicity () > 0) &&
	 (p -> Var (i) -> isDefinedInteger ()))
       si.setInteger (i);

   //si.writeMpsNative("toto",NULL,NULL,1);
   //si.writeLp ("toto");
   app_ -> enableWarmStart();

   // restored check. With "TOO FEW DEGREES OF FREEDOM" exception, x_sol() is null
   if (problem () -> x_sol ()) {
     setColSolution (problem () -> x_sol     ());
     setRowPrice    (problem () -> duals_sol ());
   }
}


/** To set some application specific defaults. */
void CouenneInterface::setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options){
  Options->SetStringValue("bonmin.algorithm", "B-Couenne", true, true);
  Options->SetIntegerValue("bonmin.filmint_ecp_cuts", 1, true, true);
}

