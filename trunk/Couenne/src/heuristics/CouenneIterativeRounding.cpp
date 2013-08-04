// (C) Copyright International Business Machines Corporation 2010
// and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Giacomo Nannicini, Tepper School of Business, Carnegie Mellon University
//
// Date : 07/25/2010

#include "CouenneIterativeRounding.hpp"
#include "BonTMINLP2Quad.hpp"
#include "BonTMINLPLinObj.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#endif

#include "CouenneRecordBestSol.hpp"

#define MILPTIME 5
#define CBCMILPTIME 20

namespace Couenne{

  CouenneIterativeRounding::CouenneIterativeRounding():
    CbcHeuristic(),
    nlp_(NULL), cinlp_(NULL), milp_(NULL),
    maxRoundingIter_(10),
    maxFirPoints_(5), maxTime_(60), maxTimeFirstCall_(60), numInitialRows_(0),
    numSol_(-1), colLower_(NULL), colUpper_(NULL), 
    colLowerNlp_(NULL), colUpperNlp_(NULL),
    omega_(0.2), baseLbRhs_(15),
    couenne_(NULL)
  {
    setHeuristicName("CouenneIterativeRounding");
  }
 
  CouenneIterativeRounding::CouenneIterativeRounding(Bonmin::OsiTMINLPInterface* nlp, 
						     OsiSolverInterface* cinlp,
						     CouenneProblem* couenne,
						     Ipopt::SmartPtr<Ipopt::OptionsList> options):
    CbcHeuristic(), nlp_(NULL),
    cinlp_(NULL), milp_(NULL),
    numSol_(-1), colLower_(NULL), colUpper_(NULL), 
    colLowerNlp_(NULL), colUpperNlp_(NULL),
    couenne_(couenne)
  {
    // Initialize Dynamic NLP
    setNlp(nlp, cinlp);

    // Read options
    int irAggressiveness;
    options->GetIntegerValue("iterative_rounding_aggressiveness",
			     irAggressiveness, "couenne.");
    setAggressiveness(irAggressiveness);
    double maxTime, maxTimeInit;
    options->GetNumericValue("iterative_rounding_time", maxTime,
			     "couenne.");
    options->GetNumericValue("iterative_rounding_time_firstcall",
			     maxTimeInit, "couenne.");
    if (maxTime >= 0){
      setMaxTime(maxTime);
    }
    if (maxTimeInit >= 0){
      setMaxTimeFirstCall(maxTimeInit);
    }
    int irLbrhs;
    options->GetIntegerValue("iterative_rounding_base_lbrhs",
			     irLbrhs, "couenne.");
    setBaseLbRhs(irLbrhs);
    int numFirPoints;
    options->GetIntegerValue("iterative_rounding_num_fir_points",
			     numFirPoints, "couenne.");
    setMaxFirPoints(numFirPoints);
    double omega;
    options->GetNumericValue("iterative_rounding_omega",
			     omega, "couenne.");
    setOmega(omega);
  }
  
  CouenneIterativeRounding::CouenneIterativeRounding(const CouenneIterativeRounding & other):
    CbcHeuristic(other), nlp_(other.nlp_), 
    cinlp_(other.cinlp_), milp_(other.milp_), 
    maxRoundingIter_(other.maxRoundingIter_),
    maxFirPoints_(other.maxFirPoints_),
    maxTime_(other.maxTime_),
    maxTimeFirstCall_(other.maxTimeFirstCall_),
    numInitialRows_(other.numInitialRows_),
    numSol_(other.numSol_),
    omega_(other.omega_), baseLbRhs_(other.baseLbRhs_),
    couenne_(other.couenne_)
  {
    if(nlp_ != NULL){
      nlp_ = dynamic_cast<Bonmin::OsiTMINLPInterface*>(other.nlp_->clone());
    }
    if(milp_ != NULL)
#ifdef COIN_HAS_CPX
      milp_ =  dynamic_cast<OsiCpxSolverInterface*>(other.milp_->clone());
#else
      milp_ =  dynamic_cast<OsiClpSolverInterface*>(other.milp_->clone());
#endif
    if (other.colLower_ != NULL){
      if (colLower_ != NULL)
	delete colLower_;
      colLower_ = new double[milp_->getNumCols()];
      CoinCopyN (other.colLower_, milp_->getNumCols(), colLower_);
    }
    if (other.colUpper_ != NULL){
      if (colUpper_ != NULL)
	delete colUpper_;
      colUpper_ = new double[milp_->getNumCols()];
      CoinCopyN (other.colUpper_, milp_->getNumCols(), colUpper_);
    }
    if (other.colLowerNlp_ != NULL){
      if (colLowerNlp_ != NULL)
	delete colLowerNlp_;
      colLowerNlp_ = new double[nlp_->getNumCols()];
      CoinCopyN (other.colLowerNlp_, nlp_->getNumCols(), colLowerNlp_);
    }
    if (other.colUpperNlp_ != NULL){
      if (colUpperNlp_ != NULL)
	delete colUpperNlp_;
      colUpperNlp_ = new double[nlp_->getNumCols()];
      CoinCopyN (other.colUpperNlp_, nlp_->getNumCols(), colLowerNlp_);
    }
  }
  
  CbcHeuristic * 
  CouenneIterativeRounding::clone() const{
    return new CouenneIterativeRounding(*this);
  }
  
  CouenneIterativeRounding &
  CouenneIterativeRounding::operator=(const CouenneIterativeRounding & rhs){
    if(this != &rhs){
      CbcHeuristic::operator=(rhs);
      if(nlp_)
        delete nlp_;
      
      if(rhs.nlp_ != NULL){
	nlp_ = dynamic_cast<Bonmin::OsiTMINLPInterface*>(rhs.nlp_->clone());
      }
      cinlp_ = rhs.cinlp_;
      maxRoundingIter_ = rhs.maxRoundingIter_;
      maxFirPoints_ = rhs.maxFirPoints_;
      maxTime_ = rhs.maxTime_;
      maxTimeFirstCall_ = rhs.maxTimeFirstCall_;
      numSol_ = rhs.numSol_;
      numInitialRows_ = rhs.numInitialRows_;
      omega_ = rhs.omega_;
      baseLbRhs_ = rhs.baseLbRhs_;
      couenne_ = rhs.couenne_;
      if (rhs.colLower_ != NULL){
	if (colLower_ != NULL)
	  delete colLower_;
	colLower_ = new double[milp_->getNumCols()];
	CoinCopyN (rhs.colLower_, milp_->getNumCols(), colLower_);
      }
      if (rhs.colUpper_ != NULL){
	if (colUpper_ != NULL)
	  delete colUpper_;
	colUpper_ = new double[milp_->getNumCols()];
	CoinCopyN (rhs.colUpper_, milp_->getNumCols(), colLower_);
      }
      if (rhs.colLowerNlp_ != NULL){
	if (colLowerNlp_ != NULL)
	  delete colLowerNlp_;
	colLowerNlp_ = new double[nlp_->getNumCols()];
	CoinCopyN (rhs.colLowerNlp_, nlp_->getNumCols(), colLowerNlp_);
      }
      if (rhs.colUpperNlp_ != NULL){
	if (colUpperNlp_ != NULL)
	  delete colUpperNlp_;
	colUpperNlp_ = new double[nlp_->getNumCols()];
	CoinCopyN (rhs.colUpperNlp_, nlp_->getNumCols(), colLowerNlp_);
      }
    }
    return *this;
  }
  
  CouenneIterativeRounding::~CouenneIterativeRounding(){
    delete nlp_;
    nlp_ = NULL;
    if (colLower_)
      delete[] colLower_;
    if (colUpper_)
      delete[] colUpper_;
    if (colLowerNlp_)
      delete[] colLowerNlp_;
    if (colUpperNlp_)
      delete[] colUpperNlp_;
    if (milp_)
      delete milp_;
    milp_ = NULL;
  }
  
  void
  CouenneIterativeRounding::setNlp(Bonmin::OsiTMINLPInterface* nlp, 
				   OsiSolverInterface * cinlp){
    // Create a dynamic NLP (i.e. one that can be used to add/remove
    // linear inequalities) from the initial one
    if(nlp_ != NULL)
      delete nlp_;
    nlp_ = new Bonmin::OsiTMINLPInterface;
    Ipopt::SmartPtr<Bonmin::TMINLP> tminlp = nlp->model();
    if (tminlp->hasLinearObjective()){
      Ipopt::SmartPtr<Bonmin::TMINLPLinObj> linObj =
	new Bonmin::TMINLPLinObj;
      linObj->setTminlp(GetRawPtr(tminlp));
      tminlp = GetRawPtr(linObj);
    }
    nlp_->initialize(nlp->regOptions(), nlp->options(), nlp->solver()->journalist(), "bonmin.", tminlp);
    nlp_->use(new Bonmin::TMINLP2TNLPQuadCuts(tminlp));
    cinlp_ = cinlp;
  }

  void
  CouenneIterativeRounding::setMilp(){
    if(milp_ != NULL)
      delete milp_;

    // Save the LP relaxation of the root node

    OsiSolverInterface * milp = model_->solver();
    int n = milp->getNumCols();

#ifdef COIN_HAS_CPX
    milp_ = new OsiCpxSolverInterface();
    milp_->loadProblem(*(milp->getMatrixByRow()), milp->getColLower(), 
		       milp->getColUpper(), milp->getObjCoefficients(),
		       milp->getRowLower(), milp->getRowUpper());
    for (int i = 0; i < n; ++i){
      if (milp->isInteger(i))
	milp_->setInteger(i); 
    }
#else
    milp_ =  dynamic_cast<OsiClpSolverInterface*>(milp->clone());
#endif

    colLower_ = new double[n];
    colUpper_ = new double[n];
    memcpy(colLower_, milp->getColLower(), n*sizeof(double));
    memcpy(colUpper_, milp->getColUpper(), n*sizeof(double));

    int nNlp = cinlp_->getNumCols();

    colLowerNlp_ = new double[nNlp];
    colUpperNlp_ = new double[nNlp];
    memcpy(colLowerNlp_, cinlp_->getColLower(), nNlp*sizeof(double));
    memcpy(colUpperNlp_, cinlp_->getColUpper(), nNlp*sizeof(double));

    numIntegers_ = 0;
    for (int i = 0; i < nNlp; ++i){
      if (cinlp_->isInteger(i)){
	numIntegers_++;
      }
    }

    // workaround for inconsistent bounds:
    // sometimes using getColLower and getColUpper we get *slightly*
    // inconsistent bounds, which make everything crash afterwards
    double swap;
    for (int i = 0; i < n; ++i){
      if (colUpper_[i] < colLower_[i]){
	swap = colUpper_[i];
	colUpper_[i] = colLower_[i];
	colLower_[i] = swap;
      }
    }

    numInitialRows_ = milp_->getNumRows();
    
    // this array contains zeroes; we use it to set empty elements
    double * tmpArray = new double[n];
    CoinFillN(tmpArray, n, 0.0);
    milp_->setObjective(tmpArray);
    milp_->setObjSense(1);

    // now there is no objective function
    // create new variables with objective coefficient 1
    for (int i = 0; i < n; ++i){
      milp_->addCol(0, NULL, NULL, 0.0, COIN_DBL_MAX, 1.0);
    }

    milp_->setHintParam(OsiDoDualInResolve,true,OsiHintDo);
    milp_->setHintParam(OsiDoPresolveInResolve,true,OsiHintDo);
    milp_->setHintParam(OsiDoReducePrint,true,OsiHintDo);
    milp_->setDblParam(OsiPrimalTolerance, COUENNE_EPS_INT);
    milp_->messageHandler()->setLogLevel(0);
    milp_->setDblParam(OsiDualTolerance, COUENNE_EPS_INT);

#ifdef COIN_HAS_CPX
    // Set options if we have Cplex
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_HIDDENFEAS);
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_FPHEUR, 2);
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_NODESEL, 0);
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_LBHEUR, 1);
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_RINSHEUR, 99);
    CPXsetintparam(milp_->getEnvironmentPtr(), CPX_PARAM_NODELIM, 50);
    CPXsetdblparam(milp_->getEnvironmentPtr(), CPX_PARAM_CUTSFACTOR, 1.0);
    CPXsetdblparam(milp_->getEnvironmentPtr(), CPX_PARAM_TILIM, MILPTIME);
#else
    // Set up heuristics if we use Cbc (not recommended)
    heuristics_ = new CbcHeuristic*[1];
    numHeuristics_ = 1;
    CbcHeuristicFPump * feaspump = new CbcHeuristicFPump();
    feaspump->setMaximumRetries(2);
    feaspump->setMaximumPasses(100);
    feaspump->setAccumulate(3);
    heuristics_[0] = feaspump;
#endif

    delete[] tmpArray;


  }
  
  void CouenneIterativeRounding::setAggressiveness(int value){
    switch (value){
    case 0:
      setMaxRoundingIter(5);
      setMaxTimeFirstCall(300);
      setMaxFirPoints(5);
      setMaxTime(60);
      break;
    case 1:
      setMaxRoundingIter(10);
      setMaxTimeFirstCall(300);
      setMaxFirPoints(5);
      setMaxTime(120);
      break;
    case 2:
      setMaxRoundingIter(20);
      setMaxTimeFirstCall(1000);
      setMaxFirPoints(5);
      setMaxTime(300);
      break;
    default:
      std::cerr << "CouenneIterativeRounding::setAggressiveness() : unknown value!\n" << std::endl;
    }
  }

  int
  CouenneIterativeRounding::solution(double & objectiveValue, double* newSolution){
    if (milp_ == NULL){
      // On the very first call, we set the MILP
      setMilp();
      return 0;
    }

    if ((model_->numberIntegers() == 0) || 
	(numSol_ == model_->getSolutionCount())){
      // There are no integers, or we already tried to improve current
      // solution. Exit immediately.
      return 0;
    }

    numSol_ = model_->getSolutionCount();

    std::cout << "Launching IterativeRounding with parameters:" << std::endl;
    std::cout << "Max rounding iter: " << maxRoundingIter_ << std::endl;
    std::cout << "Max feas point: " << maxFirPoints_ << std::endl;
    std::cout << "Base lbrhs: " << baseLbRhs_ << std::endl;
    std::cout << "Omega: " << omega_ << std::endl;
    std::cout << "Max time firstcall: " << maxTimeFirstCall_ << std::endl;

    // write down starting time and apply heuristics
    startTime_ = CoinCpuTime();
    endTime_ = ((numSol_ == 0) ? maxTimeFirstCall_ : maxTime_);

    const double* bestKnownSolution = model_->bestSolution();
    bool found = false;
    bool hasSolution = true;
    bool improved = true;
    if (numSol_ == 0){
      // No solution available; start with our feasibility heuristic
      hasSolution = feasibilityIR(objectiveValue, newSolution);
      if (hasSolution){
	bestKnownSolution = newSolution;
	found = hasSolution;
      }
    }
    if (!hasSolution){
      // Still do not have a starting solution, we cannot run improvementIR
      return found;
    }
    while (improved && (CoinCpuTime()-startTime_) < (endTime_ - MILPTIME)){
      // Keep improving current solution as long as it works
      improved = false;
      improved = improvementIR(objectiveValue, newSolution, bestKnownSolution);
      if (improved){
	bestKnownSolution = newSolution;
      }
      found = (found || improved);
    }
    if (found){
      // make sure we do not run improvementIR again on the same solution
      numSol_++;
    }

    return found;
  }

  int
  CouenneIterativeRounding::feasibilityIR(double& objectiveValue, 
					    double* newSolution){

    std::cout << "starting feasibility IR" << std::endl;

    OsiSolverInterface * solver = model_->solver();

    OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
    Bonmin::BabInfo * babInfo = dynamic_cast<Bonmin::BabInfo *> (auxInfo);

    if(babInfo){
      babInfo->setHasNlpSolution(false);
      if(babInfo->infeasibleNode()){
	return 0;
      }
    }
   
    int n = couenne_->nVars();
    int nNlp = cinlp_->getNumCols();

    // get best solution and count the number of binary variables
    int numIntAtBound = 0;

    OsiRowCut cut;
    OsiRowCut lbcut1;
    OsiRowCut lbcut2;
    OsiCuts lbcuts;
    double obj;

    bool boundsChanged = false;
    std::vector<int> previousBranches;

    // Try to find a feasible solution

    // apply NLP solver to obtain point xprime
    nlp_->resolve();
    if (nlp_->isProvenPrimalInfeasible() || nlp_->isProvenDualInfeasible() ||
	nlp_->isAbandoned()){
      nlp_->resolveForRobustness(3);
    }
    
    if (nlp_->isProvenPrimalInfeasible() || nlp_->isProvenDualInfeasible() ||
	nlp_->isAbandoned()){
      obj = COIN_DBL_MAX;
    }
    else{
      obj = nlp_->getObjValue();
    }
    
    if (obj == COIN_DBL_MAX){std::cout << " IR: no feasible solution found " << std::endl;
      std::cout << " IR: elapsed time " << CoinCpuTime()-startTime_ << std::endl;
      return 0;
    }

    double* xprime = new double [n];
    CoinCopyN (nlp_->getColSolution(), nlp_->getNumCols(), xprime);
    couenne_ -> getAuxs (xprime);

    // now prepare the MILP that we will solve

    // the objective is : minimize sum_{j = 1}^n w_j
    int constrInd[2];
    double constrElem[2];

    // main loop: solve a sequence of MILPs and NLPs
    bool foundSolution = false;
    double* tmpSolution = new double[n];
    OsiCuts revlb_all;
    OsiSolverInterface* curr_milp = milp_;

    int outerLoop = maxFirPoints_;
    for (int h = 0; h < outerLoop && 
	   ((CoinCpuTime()-startTime_) < (endTime_ - MILPTIME)); ++h){

      OsiCuts cs;
      cs.insert(lbcut1);
      cs.insert(lbcut2);
      // write constraints:
      // xbar[j] - x[j] <= w_j
      // x[j] - xbar[j] <= w_j
      for (int i = 0; i < n; i++){
	constrInd[0] = i;
	constrInd[1] = i+n;
	constrElem[0] = -1;
	constrElem[1] = -1;
	cut.mutableRow().setVector(2, constrInd, constrElem);
	cut.setLb(-COIN_DBL_MAX);
	cut.setUb(-xprime[i]);
	cs.insert(cut);
	constrElem[0] = +1;
	constrElem[1] = -1;
	cut.mutableRow().setVector(2, constrInd, constrElem);
	cut.setLb(-COIN_DBL_MAX);
	cut.setUb(xprime[i]);
	cs.insert(cut);
      }
      curr_milp->applyCuts(cs);

      for (int k = 0; k < maxRoundingIter_ &&
	     ((CoinCpuTime()-startTime_) < (endTime_ - MILPTIME)); ++k){
	// Solve MILP to obtain an integral feasible point near xprime
	curr_milp->restoreBaseModel(numInitialRows_+cs.sizeRowCuts());
	curr_milp->applyCuts(revlb_all);
	bool solFound = solveMilp(curr_milp, 
				  endTime_-(CoinCpuTime()-startTime_)-2);
	if (!solFound && !boundsChanged){
	  std::cout << " MILP cannot be solved, terminating LB " << std::endl;
	  // we cannot solve the MILP and bounds were not changed; exit
	  curr_milp->restoreBaseModel(numInitialRows_);
	  delete[] xprime;
	  delete[] tmpSolution;
	  if (boundsChanged){
	    curr_milp->setColLower(colLower_);
	    curr_milp->setColUpper(colUpper_);
	  }
	  std::cout << " IR: elapsed time " << CoinCpuTime()-startTime_ << std::endl;
	  return foundSolution;
	}
	else if (!solFound && boundsChanged){
	  // the MILP is infeasible but bounds were changed;
	  // restore original bounds and branch on a random variable
	  // of the solution at the previous iteration (so we do not cycle)
	  curr_milp->setColLower(colLower_);
	  curr_milp->setColUpper(colUpper_);
	  branchToCut(tmpSolution, curr_milp, previousBranches);
	  continue;
	}

	const double * xtilde = curr_milp->getColSolution();

	// now fix integer variables and solve NLP
	// xtilde has more columns than the NLP, so this should be ok
	CoinCopyN (xtilde, n, tmpSolution);
	for (int i = 0; i < nNlp; ++i){
	  if (model_->isInteger(i)){
	    tmpSolution[i] = floor(tmpSolution[i]+0.5);
	    cinlp_->setColLower(i, tmpSolution[i]);
	    cinlp_->setColUpper(i, tmpSolution[i]);
	  }
	}
	cinlp_->setColSolution(tmpSolution);

	cinlp_->messageHandler()->setLogLevel(1);
	cinlp_->resolve();
	obj = ((cinlp_->isProvenOptimal()) ? 
	       cinlp_->getObjValue():COIN_DBL_MAX);
	memcpy(tmpSolution, cinlp_->getColSolution(), 
	       nNlp*sizeof(double));
	// check solution of the NLP;
	// if we have a new incumbent we are done, otherwise we iterate

	bool isChecked = false;

	isChecked = couenne_ -> checkNLP0 (tmpSolution, obj, true,
					   false,
					   true,
					   false);

// #ifdef FM_CHECKNLP2
// 	isChecked = couenne_->checkNLP2(tmpSolution, 0, false, // do not care about obj
// 					true, // stopAtFirstViol
// 					false, // checkALL
// 					couenne_->getFeasTol());
// 	if(isChecked) {
// 	  obj = couenne_->getRecordBestSol()->getModSolVal();
// 	}
// #else /* not FM_CHECKNLP2 */
// 	isChecked = couenne_->checkNLP(tmpSolution, obj, true);
// #endif  /* not FM_CHECKNLP2 */
	
	if (cinlp_->isProvenOptimal () &&
	    isChecked &&
	    (obj < couenne_->getCutOff())) {
	  
#ifdef FM_CHECKNLP2
#ifdef FM_TRACE_OPTSOL
	  couenne_->getRecordBestSol()->update();
	  CoinCopyN (couenne_->getRecordBestSol()->getSol(), n, tmpSolution);
	  obj = couenne_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
	  CoinCopyN (couenne_->getRecordBestSol()->getModSol(n), n, tmpSolution);
#endif /* not FM_TRACE_OPTSOL */
#else /* not FM_CHECKNLP2 */
	  
	  //Get correct values for all auxiliary variables
	  couenne_ -> getAuxs (tmpSolution);
	  
#ifdef FM_TRACE_OPTSOL
	  couenne_->getRecordBestSol()->update(tmpSolution, n,
					       obj, couenne_->getFeasTol());
	  CoinCopyN (couenne_->getRecordBestSol()->getSol(), n, tmpSolution);
	  obj = couenne_->getRecordBestSol()->getVal();
#endif /* FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */
	  
	  if (babInfo){
	    babInfo->setNlpSolution (tmpSolution, n, obj);
	    babInfo->setHasNlpSolution (true);
	  }
	  
	  std::cout << "Final Nlp solution with objective " << obj << " :" << std::endl;
	  
	  if (obj < objectiveValue - COUENNE_EPS) { // found better solution?
	    std::cout << "New incumbent found" << std::endl;
	    const CouNumber 
	      *lb = solver -> getColLower (),
	      *ub = solver -> getColUpper ();
	    
	    // check bounds once more after getAux. This avoids false
	    // asserts in CbcModel.cpp:8305 on integerTolerance violated
	    for (int i=0; i < n; ++i, ++lb, ++ub) {
	      CouNumber &t = tmpSolution [i];
	      if      (t < *lb) t = *lb;
	      else if (t > *ub) t = *ub;
	    }

	    couenne_ -> setCutOff (obj);
	    foundSolution = true;
	    objectiveValue = obj;
	    CoinCopyN (tmpSolution, n, newSolution);
	    cinlp_->setColLower(colLowerNlp_);
	    cinlp_->setColUpper(colUpperNlp_);
	  
	    curr_milp->restoreBaseModel(numInitialRows_);
	    delete[] xprime;
	    delete[] tmpSolution;
	    if (boundsChanged){
	      curr_milp->setColLower(colLower_);
	      curr_milp->setColUpper(colUpper_);
	    }
	    //delete curr_milp;
	    double elapsed = CoinCpuTime()-startTime_;
	    std::cout << "IR: Heuristic: " << objectiveValue << std::endl;
	    std::cout << "IR: Heuristic time: " << elapsed << std::endl;
	    return foundSolution;
	  }
	}
	cinlp_->setColLower(colLowerNlp_);
	cinlp_->setColUpper(colUpperNlp_);
	double avgBound;
	numIntAtBound = computeIntAtBound(xtilde, avgBound);
	
	if (numIntAtBound >= 50 || 
	    numIntAtBound >= ((numIntegers_*0.1 > 5) ? numIntegers_*0.1 : 5)){
	  // write reverse local branching constraint
	  // this avoids finding the same xtilde again
	  avgBound = floor(avgBound + 0.5);
	  std::cout << "Using reverse LB with rhs " << avgBound << std::endl;
	  writeLB(cut, xtilde, 'G', avgBound);
	  revlb_all.insert(cut);
	}
	else{
	  // Branch on a random variable to obtain a different integer point
	  branchToCut(xtilde, curr_milp, previousBranches);
	  boundsChanged = true;
	}
      }
      if (h <= outerLoop -2){
	// Compute a NLP-feasible point by solving a log-barrier problem
	// with a given minimum value for the log-barrier parameter mu
	Bonmin::OsiTMINLPInterface * nlp = dynamic_cast<Bonmin::OsiTMINLPInterface *> (nlp_);
	Ipopt::SmartPtr< Ipopt::OptionsList > opt = nlp->solver()->options();
	nlp->messageHandler()->setLogLevel(10);
	double mu_target;
	double kappa_d;
	opt->GetNumericValue("mu_target", mu_target, "ipopt.");
	opt->SetNumericValue("mu_target", omega_*(h+1), "ipopt.");
	opt->GetNumericValue("kappa_d", kappa_d, "ipopt.");
	opt->SetNumericValue("kappa_d", 0.0, "ipopt.");
	nlp_->resolve();
	if (nlp_->isProvenPrimalInfeasible() || 
	    nlp_->isProvenDualInfeasible() ||
	    nlp_->isAbandoned()){
	  nlp_->resolveForRobustness(3);
	}
	opt->SetNumericValue("mu_target", mu_target, "ipopt.");
	opt->SetNumericValue("kappa_d", kappa_d, "ipopt.");

	if (!nlp->isProvenPrimalInfeasible() &&
	    !nlp->isProvenDualInfeasible() && 
	    !nlp->isAbandoned()){
	  CoinCopyN (nlp_->getColSolution(), nlp_->getNumCols(), xprime);
	  couenne_ -> getAuxs (xprime);
	  curr_milp->restoreBaseModel(numInitialRows_);
	  if (boundsChanged){
	    curr_milp->setColLower(colLower_);
	    curr_milp->setColUpper(colUpper_);
	  }
	}
      }
    }
    
    curr_milp->restoreBaseModel(numInitialRows_);
    delete[] xprime;
    delete[] tmpSolution;
    if (boundsChanged){
      curr_milp->setColLower(colLower_);
      curr_milp->setColUpper(colUpper_);
    }
    double elapsed = CoinCpuTime()-startTime_;
    std::cout << "IR: Heuristic: " << COUENNE_INFINITY << std::endl;
    std::cout << "IR: Heuristic time: " << elapsed << std::endl;
    return foundSolution;
  }

  int
  CouenneIterativeRounding::improvementIR(double& objectiveValue,
					    double* newSolution,
					    const double* solution){
    std::cout << "starting Improvement IR" << std::endl;

    OsiSolverInterface * solver = model_->solver();

    OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
    Bonmin::BabInfo * babInfo = dynamic_cast<Bonmin::BabInfo *> (auxInfo);

    if(babInfo){
      babInfo->setHasNlpSolution(false);
      if(babInfo->infeasibleNode()){
	return 0;
      }
    }


    int n = couenne_->nVars();
    int nNlp = cinlp_->getNumCols();

    // count the number of binary variables
    int numIntAtBound = 0;

    double lbrhs = CoinMin(baseLbRhs_, CoinMax(1,numIntegers_/2) );
    OsiRowCut cut;
    OsiRowCut lbcut1;
    OsiRowCut lbcut2;
    OsiCuts lbcuts;
    int currentIndex = 0;
    double obj;

    bool boundsChanged = false;
    std::vector<int> previousBranches;

    double avgBound = 0.0;
    numIntAtBound = computeIntAtBound(solution, avgBound);

    if (numIntAtBound >= 50 ||
	numIntAtBound >= ((numIntegers_*0.1 > 5) ? numIntegers_*0.1 : 5)){
      // write local branching constraint
      writeLB(lbcut1, solution, 'L', lbrhs + floor(avgBound - 0.5));
      lbcuts.insert(lbcut1);
      writeLB(lbcut2, solution, 'G', 1);
      lbcuts.insert(lbcut2);

      // We cannot add a row because OsiTMINLP will not let us do that!
      // Therefore we add a cut -- does the same thing.
      // We must delete the cut afterwards!
      nlp_->applyCuts(lbcuts);
    }
    else{
      // we cannot use a local branching constraint, so we simply
      // reduce the size of the box bounds, to have a smaller NLP and
      // stay somehow "close" to the current point
      for (int i = 0; i < nlp_->getNumCols(); ++i){
	if (nlp_->isInteger(i)){
	  nlp_->setColLower(i, colLowerNlp_[i]+(solution[i]-colLower_[i])*0.5);
	  nlp_->setColUpper(i, colUpperNlp_[i]+(solution[i]-colUpper_[i])*0.5);
	}
      }
    }
    // apply NLP solver to obtain point xprime
    nlp_->setColSolution(solution);
    nlp_->resolve();
    if (nlp_->isProvenPrimalInfeasible() || nlp_->isProvenDualInfeasible() ||
	nlp_->isAbandoned()){
      nlp_->resolveForRobustness(3);
    }

    if (nlp_->isProvenPrimalInfeasible() || nlp_->isProvenDualInfeasible() ||
	nlp_->isAbandoned()){
      obj = COIN_DBL_MAX;
    }
    else{
      obj = nlp_->getObjValue();
    }

    // Restore NLP
    if (numIntAtBound > 0){
      currentIndex = nlp_->getNumRows()-1;
      nlp_->deleteRows(1, &currentIndex);
      currentIndex = nlp_->getNumRows()-1;
      nlp_->deleteRows(1, &currentIndex);
    }
    else{
      nlp_->setColLower(colLowerNlp_);
      nlp_->setColUpper(colUpperNlp_);
    }
      
    if (obj == COIN_DBL_MAX || obj >= objectiveValue - COUENNE_EPS){
      std::cout << " IR: no improvement possible " << std::endl;
      std::cout << " IR: elapsed time " << CoinCpuTime()-startTime_ << std::endl;
      return 0;
    }

    double* xprime = new double [n];
    CoinCopyN (nlp_->getColSolution(), nlp_->getNumCols(), xprime);
    couenne_ -> getAuxs (xprime);

    // now prepare the MILP that we will solve

    // the objective is : minimize sum_{j = 1}^n w_j
    int constrInd[2];
    double constrElem[2];

    // main loop: solve a sequence of MILPs and NLPs
    bool foundSolution = false;
    double* tmpSolution = new double[n];
    OsiCuts revlb_all;
    OsiSolverInterface* curr_milp = milp_;
    
    OsiCuts cs;
    cs.insert(lbcut1);
    cs.insert(lbcut2);
    // write constraints:
    // xbar[j] - x[j] <= w_j
    // x[j] - xbar[j] <= w_j
    for (int i = 0; i < n; i++){
      constrInd[0] = i;
      constrInd[1] = i+n;
      constrElem[0] = -1;
      constrElem[1] = -1;
      cut.mutableRow().setVector(2, constrInd, constrElem);
      cut.setLb(-COIN_DBL_MAX);
      cut.setUb(-xprime[i]);
      cs.insert(cut);
      constrElem[0] = +1;
      constrElem[1] = -1;
      cut.mutableRow().setVector(2, constrInd, constrElem);
      cut.setLb(-COIN_DBL_MAX);
      cut.setUb(xprime[i]);
      cs.insert(cut);
    }
    curr_milp->applyCuts(cs);

    for (int k = 0; k < maxRoundingIter_ &&
	   ((CoinCpuTime()-startTime_) < (endTime_ - MILPTIME)); ++k){
      // Solve MILP to obtain an integral feasible point near xprime
      curr_milp->restoreBaseModel(numInitialRows_+cs.sizeRowCuts());
      curr_milp->applyCuts(revlb_all);
      bool solFound = solveMilp(curr_milp, 
				endTime_-(CoinCpuTime()-startTime_)-2);
      if (!solFound && !boundsChanged){
	std::cout << " MILP cannot be solved, terminating LB " << std::endl;
	// we cannot solve the MILP and bounds were not changed; exit
	curr_milp->restoreBaseModel(numInitialRows_);
	delete[] xprime;
	delete[] tmpSolution;
	if (boundsChanged){
	  curr_milp->setColLower(colLower_);
	  curr_milp->setColUpper(colUpper_);
	}
	std::cout << " IR: elapsed time " << CoinCpuTime()-startTime_ << std::endl;
	return foundSolution;
      }
      else if (!solFound && boundsChanged){
	// the MILP is infeasible but bounds were changed;
	// restore original bounds and branch on a random variable
	// of the solution at the previous iteration (so we do not cycle)
	curr_milp->setColLower(colLower_);
	curr_milp->setColUpper(colUpper_);
	branchToCut(tmpSolution, curr_milp, previousBranches);
	continue;
      }
      const double * xtilde = curr_milp->getColSolution();

      // now fix integer variables and solve NLP
      // xtilde has more columns than the NLP, so this should be ok
      CoinCopyN (xtilde, n, tmpSolution);
      for (int i = 0; i < nNlp; ++i){
	if (model_->isInteger(i)){
	  tmpSolution[i] = floor(tmpSolution[i]+0.5);
	  cinlp_->setColLower(i, tmpSolution[i]);
	  cinlp_->setColUpper(i, tmpSolution[i]);
	}
      }
      cinlp_->setColSolution(tmpSolution);

      cinlp_->messageHandler()->setLogLevel(1);
      cinlp_->resolve();
      obj = ((cinlp_->isProvenOptimal()) ? cinlp_->getObjValue():COIN_DBL_MAX);
      memcpy(tmpSolution, cinlp_->getColSolution(), 
	     nNlp*sizeof(double));
      // check solution of the NLP;
      // if we have a new incumbent we are done, otherwise we iterate

      bool isChecked = false;

      isChecked = couenne_ -> checkNLP0 (tmpSolution, obj, true,
					 false,
					 true,
					 false);

// #ifdef FM_CHECKNLP2
//       isChecked = couenne_->checkNLP2(tmpSolution, 0, false, // do not care about obj
// 				      true, // stopAtFirstViol
// 				      false, // checkALL
// 				      couenne_->getFeasTol());
//       if(isChecked) {
// 	obj = couenne_->getRecordBestSol()->getModSolVal();
//       }
// #else /* not FM_CHECKNLP2 */
//       isChecked = couenne_->checkNLP(tmpSolution, obj, true);
// #endif  /* not FM_CHECKNLP2 */
      
      if (cinlp_->isProvenOptimal () &&
	  isChecked &&
	  (obj < couenne_->getCutOff())) {
	
#ifdef FM_CHECKNLP2
#ifdef FM_TRACE_OPTSOL
	couenne_->getRecordBestSol()->update();
	CoinCopyN (couenne_->getRecordBestSol()->getSol(), n, tmpSolution);
	obj = couenne_->getRecordBestSol()->getVal();
#else /* not FM_TRACE_OPTSOL */
	CoinCopyN (couenne_->getRecordBestSol()->getModSol(n), n, tmpSolution);
#endif /* not FM_TRACE_OPTSOL */
#else /* not FM_CHECKNLP2 */
	
	//Get correct values for all auxiliary variables
	couenne_ -> getAuxs (tmpSolution);
	
#ifdef FM_TRACE_OPTSOL
	couenne_->getRecordBestSol()->update(tmpSolution, n,
					     obj, couenne_->getFeasTol());
	CoinCopyN (couenne_->getRecordBestSol()->getSol(), n, tmpSolution);
	obj = couenne_->getRecordBestSol()->getVal();
#endif /* FM_TRACE_OPTSOL */
#endif /* not FM_CHECKNLP2 */
	
	if (babInfo){
	  babInfo->setNlpSolution (tmpSolution, n, obj);
	  babInfo->setHasNlpSolution (true);
	}

	std::cout << "Final Nlp solution with objective " << obj << " :" << std::endl;

	if (obj < objectiveValue - COUENNE_EPS) { // found better solution?
	  std::cout << "New incumbent found" << std::endl;
	  const CouNumber 
	    *lb = solver -> getColLower (),
	    *ub = solver -> getColUpper ();

	  // check bounds once more after getAux. This avoids false
	  // asserts in CbcModel.cpp:8305 on integerTolerance violated
	  for (int i=0; i < n; ++i, ++lb, ++ub) {
	    CouNumber &t = tmpSolution [i];
	    if      (t < *lb) t = *lb;
	    else if (t > *ub) t = *ub;
	  }

	  couenne_ -> setCutOff (obj);
	  foundSolution = true;
	  objectiveValue = obj;
	  CoinCopyN (tmpSolution, n, newSolution);
	  cinlp_->setColLower(colLowerNlp_);
	  cinlp_->setColUpper(colUpperNlp_);
	  
	  curr_milp->restoreBaseModel(numInitialRows_);
	  delete[] xprime;
	  delete[] tmpSolution;
	  if (boundsChanged){
	    curr_milp->setColLower(colLower_);
	    curr_milp->setColUpper(colUpper_);
	  }
	  //delete curr_milp;
	  double elapsed = CoinCpuTime()-startTime_;
	  std::cout << "IR: Heuristic: " << objectiveValue << std::endl;
	  std::cout << "IR: Heuristic time: " << elapsed << std::endl;
	  return foundSolution;
	}
      }
      cinlp_->setColLower(colLowerNlp_);
      cinlp_->setColUpper(colUpperNlp_);
      numIntAtBound = computeIntAtBound(xtilde, avgBound);
	
      if (numIntAtBound >= 50 || 
	  numIntAtBound >= ((numIntegers_*0.1 > 5) ? numIntegers_*0.1 : 5)){
	// write reverse local branching constraint
	// this avoids finding the same xtilde again
	avgBound = floor(avgBound + 0.5);
	std::cout << "Using reverse LB with rhs " << avgBound << std::endl;
	writeLB(cut, xtilde, 'G', avgBound);
	revlb_all.insert(cut);
      }
      else{
	// Branch (randomly) so that we can obtain a new integer point
	branchToCut(xtilde, curr_milp, previousBranches);
	boundsChanged = true;
      }
    }
    
    curr_milp->restoreBaseModel(numInitialRows_);
    delete[] xprime;
    delete[] tmpSolution;
    if (boundsChanged){
      curr_milp->setColLower(colLower_);
      curr_milp->setColUpper(colUpper_);
    }
    double elapsed = CoinCpuTime()-startTime_;
    std::cout << "IR: Heuristic: " << COUENNE_INFINITY << std::endl;
    std::cout << "IR: Heuristic time: " << elapsed << std::endl;
    return foundSolution;
  }

  int 
  CouenneIterativeRounding::computeIntAtBound(const double* x){
    int numIntAtBound = 0;
    for (int i = 0; i < nlp_->getNumCols(); ++i){
      if (nlp_->isInteger(i) && (areEqual(x[i], colLower_[i]) || 
				 areEqual(x[i], colUpper_[i]))){
	numIntAtBound++;
      }
    }
    return numIntAtBound;
  }

  int 
  CouenneIterativeRounding::computeIntAtBound(const double* x, 
						double& avgBoundSize){
    int numIntAtBound = 0;
    avgBoundSize = 0;
    for (int i = 0; i < nlp_->getNumCols(); ++i){
      if (nlp_->isInteger(i) && (areEqual(x[i], colLower_[i]) || 
				 areEqual(x[i], colUpper_[i]))){
	numIntAtBound++;
	avgBoundSize += colUpper_[i] - colLower_[i];
      }
    }
    avgBoundSize /= numIntAtBound;
    return numIntAtBound;
  }

  void
  CouenneIterativeRounding::writeLB(OsiRowCut& cut, const double* x, 
				    char sense, double rhs){
    cut.mutableRow().clear();
    for (int i = 0; i < nlp_->getNumCols(); ++i){
      if (nlp_->isInteger(i)){
	if (areEqual(x[i], colUpper_[i])){
	  cut.mutableRow().insert(i, -1);
	  rhs -= x[i];
	}
	else if (areEqual(x[i], colLower_[i])){
	  cut.mutableRow().insert(i, 1);
	  rhs += x[i];
	}
      }
    }
    if (sense == 'L'){
      cut.setLb(-COIN_DBL_MAX);
      cut.setUb(rhs);
    }
    else if (sense == 'G'){
      cut.setLb(rhs);
      cut.setUb(COIN_DBL_MAX);
    }
    else{
      std::cerr << "### ERROR: wrong sense of LB constraint" << std::endl;
      exit(1);
    }
  }

  bool
  CouenneIterativeRounding::solveMilp(OsiSolverInterface* milp,
				      double maxTime){
    double start = CoinCpuTime();
#ifdef COIN_HAS_CPX
    OsiCpxSolverInterface * solver = dynamic_cast<OsiCpxSolverInterface*>(milp);
    CPXENVptr env = solver->getEnvironmentPtr();
    CPXLPptr cpxlp = solver->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
    int status = 0;
    bool solFound = false;
    bool infeasible = false;
    while (!solFound && !infeasible && ((CoinCpuTime() - start) < maxTime)){
      solver->branchAndBound();
      status = CPXgetstat(env, cpxlp);
      solFound = ((status == CPXMIP_NODE_LIM_FEAS) ||
		  (status == CPXMIP_TIME_LIM_FEAS) ||
		  (status == CPXMIP_MEM_LIM_FEAS) ||
		  (status == CPXMIP_ABORT_FEAS) ||
		  (status == CPXMIP_OPTIMAL) ||
		  (status == CPXMIP_OPTIMAL_TOL) ||		     
		  (status == CPXMIP_ABORT_FEAS) ||
		  (status == CPXMIP_FAIL_FEAS) ||
		  (status == CPXMIP_FAIL_FEAS_NO_TREE) ||
		  (status == CPXMIP_FEASIBLE));
      infeasible = ((status == CPXMIP_INForUNBD) ||
		    (status == CPXMIP_OPTIMAL_INFEAS) ||
		    (status == CPXMIP_INFEASIBLE));
    }
    if (solFound){
      return true;
    }
    return false;
#else
    CbcModel cbcModel(*milp);
    for (int i = 0; i < numHeuristics_; ++i){
      cbcModel.addHeuristic(heuristics_[i]);
    }
    cbcModel.setMaximumSeconds(CBCMILPTIME);
    //CglPreProcess * prep = cbcModel.preProcess(0,2,5);
    while ((cbcModel.getSolutionCount() == 0) && 
	   (!cbcModel.isProvenInfeasible()) &&
	   (!cbcModel.isProvenDualInfeasible()) &&
	   (!cbcModel.isAbandoned()) &&
	   ((CoinCpuTime() - start) < maxTime)){
      cbcModel.branchAndBound();
    }
    //cbcModel.postProcess(prep);
    milp = cbcModel.solver();
    if (cbcModel.getSolutionCount() > 0){
      return true;
    }
    return false;
#endif
  }

  int
  CouenneIterativeRounding::branchToCut(const double* x, 
					OsiSolverInterface* solver,
					std::vector<int>& previousBranches){
    int branch;
    bool found = false;
    while (!found){
      branch = rand()%numIntegers_;
      found = true;
      for (unsigned int i = 0; i < previousBranches.size(); ++i){
	if (branch == previousBranches[i]){
	  found = false;
	  break;
	}
      }
      if (found){
	previousBranches.push_back(branch);
      }
      else{
	continue;
      }
      for (int i = 0; i < nlp_->getNumCols(); ++i){
	if (model_->isInteger(i)){
	  if (branch == 0){
	    branch = i;
	    break;
	  }
	  else{
	    branch--;
	  }
	}
      }
    }
    double sample = rand();
    if (sample <= ((x[branch]-colLower_[branch])/(colUpper_[branch]-colLower_[branch]))){
      solver->setColUpper(branch, x[branch]-1);
    }
    else{
      solver->setColLower(branch, x[branch]+1);
    }
    return branch;
  }

  void 
  CouenneIterativeRounding::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions){
    roptions -> AddStringOption2
      ("iterative_rounding_heuristic",
       "Do we use the Iterative Rounding heuristic",
       "no",
       "no","",
       "yes","",
       "If enabled, a heuristic based on Iterative Rounding is used "
       "to find feasible solutions for the problem. "
       "The heuristic may take some time, but usually finds good solutions. "
       "Recommended if you want good upper bounds and have Cplex. "
       "Not recommended if you do not have Cplex");

    roptions -> AddNumberOption
      ("iterative_rounding_time",
       "Specify the maximum time allowed for the Iterative Rounding heuristic",
       -1, "Maximum CPU time employed by the Iterative Rounding heuristic; "
       "if no solution found in this time, failure is reported. "
       "This overrides the CPU time set by Aggressiveness if positive.");

    roptions -> AddNumberOption
      ("iterative_rounding_time_firstcall",
       "Specify the maximum time allowed for the Iterative Rounding heuristic "
       "when no feasible solution is known",
       -1, "Maximum CPU time employed by the Iterative Rounding heuristic "
       "when no solution is known; if no solution found in this time, "
       "failure is reported."
       "This overrides the CPU time set by Aggressiveness if  posive.");
    
    roptions -> AddBoundedIntegerOption 
      ("iterative_rounding_aggressiveness",
       "Aggressiveness of the Iterative Rounding heuristic",
       0, 2, 1,
       "Set the aggressiveness of the heuristic; i.e., how many iterations "
       "should be run, and with which parameters. The maximum time can be "
       "overridden by setting the _time and _time_firstcall options. "
       "0 = non aggressive, 1 = standard (default), 2 = aggressive.");

    roptions -> AddLowerBoundedIntegerOption 
      ("iterative_rounding_num_fir_points",
       "Max number of points rounded at the beginning of Iterative Rounding",
       1, 5,
       "Number of different points (obtained solving a log-barrier problem) "
       "that the heuristic will try to round at most, during its execution "
       "at the root node (i.e. the F-IR heuristic). Default 5.");

    roptions -> AddBoundedNumberOption 
      ("iterative_rounding_omega",
       "Omega parameter of the Iterative Rounding heuristic",
       0, true, 1, true, 0.2,
       "Set the omega parameter of the heuristic, which represents a "
       "multiplicative factor for the minimum log-barrier parameter "
       "of the NLP which is solved to obtain feasible points. This "
       "corresponds to $\\omega'$ in the paper. Default 0.2.");

    roptions -> AddLowerBoundedIntegerOption 
      ("iterative_rounding_base_lbrhs",
       "Base rhs of the local branching constraint for Iterative Rounding",
       0, 15,
       "Base rhs for the local branching constraint that defines a "
       "neighbourhood of the local incumbent. The base rhs is modified by "
       "the algorithm according to variable bounds. This corresponds to "
       "k' in the paper. Default 15.");

  }

}/** Ends namespace Couenne.*/
