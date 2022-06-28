/*
 *
 * Name:    CouenneBab.cpp
 * Author:  Pietro Belotti
 * Purpose: B&B object  
 * Created: 2012-01-25
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "BonCbc.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonBabInfos.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareActual.hpp"
#include "CbcCompareObjective.hpp"
#include "CbcCompareEstimate.hpp"

#include "BonExitCodes.hpp"

#include "BonChooseVariable.hpp"
#include "BonGuessHeuristic.hpp"

#include "BonDiver.hpp"
#include "BonLinearCutsGenerator.hpp"
#include "BonTMINLPLinObj.hpp"

#include "CouenneBab.hpp"
#include "CouenneProblem.hpp"
#include "CouenneRecordBestSol.hpp"

// sets cutoff a bit above real one, to avoid single-point feasible sets
#define CUTOFF_TOL 1e-6

#define SIGNAL
#ifdef SIGNAL

#include "CoinSignal.hpp"

// Code to enable user interruption
static CbcModel * currentBranchModel = NULL; // pointer to the main b&b
extern BONMINLIB_EXPORT Bonmin::OACutGenerator2 *currentOA;   // pointer to the OA generator
extern BONMINLIB_EXPORT CbcModel                *OAModel;     // pointer to the submip if using Cbc

extern "C" {

  static bool interruptedOnce = false;

  static void couenne_signal_handler (int whichSignal) {

    if (interruptedOnce) {
      std::cerr<<"[BREAK]"<<std::endl;
      //abort ();
      exit (-1);
    } else 
      std::cerr<<"Ctrl+C detected, stopping Couenne..." << std::endl;

    if (currentBranchModel) currentBranchModel -> sayEventHappened();           // stop at next node
    if (OAModel)            OAModel   -> sayEventHappened();                      // stop at next node
    if (currentOA)          currentOA -> parameter ().maxLocalSearchTime_ = 0.; // stop OA

    interruptedOnce = true;
  }
}
#endif

using namespace Couenne;
using namespace Bonmin;

/** Constructor.*/
CouenneBab::CouenneBab (): Bab (), problem_(NULL) {}

/** Destructor.*/
CouenneBab::~CouenneBab () {}

void CouenneBab::setProblem (CouenneProblem *p)
{problem_ = p;}


/** Perform a branch-and-bound on given setup.*/
void CouenneBab::branchAndBound (Bonmin::BabSetupBase & s) {

  double remaining_time = s.getDoubleParameter(Bonmin::BabSetupBase::MaxTime) + CoinCpuTime();

  /* Put a link to this into solver.*/
  OsiBabSolver *  babInfo = dynamic_cast<OsiBabSolver *>(s.continuousSolver()->getAuxiliaryInfo());
  assert(babInfo);
  Bonmin::BabInfo *  bonBabInfoPtr = dynamic_cast<Bonmin::BabInfo*>(babInfo);

  if (bonBabInfoPtr == NULL) { //Replace with a Bonmin::babInfo
    bonBabInfoPtr = new Bonmin::BabInfo(*babInfo);
    s.continuousSolver()->setAuxiliaryInfo(bonBabInfoPtr);
    delete bonBabInfoPtr;
    bonBabInfoPtr = dynamic_cast<Bonmin::BabInfo*>(s.continuousSolver()->getAuxiliaryInfo());
  }

  bonBabInfoPtr->setBabPtr(this);

  s.nonlinearSolver()->solver()->setup_global_time_limit(s.getDoubleParameter(Bonmin::BabSetupBase::MaxTime));
  OsiSolverInterface * solver = s.continuousSolver()->clone();
  delete modelHandler_;
  modelHandler_ = s.continuousSolver()->messageHandler()->clone();
  model_.passInMessageHandler(modelHandler_);
  model_.assignSolver(solver, true);

  //  s.continuousSolver() = model_.solver();
  //   if(s.continuousSolver()->objects()!=NULL){
  //     model_.addObjects(s.continuousSolver()->numberObjects(),s.continuousSolver()->objects());
  //   }

  int specOpt = s.getIntParameter(Bonmin::BabSetupBase::SpecialOption);
  if (specOpt) {
    model_.setSpecialOptions(specOpt);
    if (specOpt==16) {
      Bonmin::CbcNlpStrategy strat(s.getIntParameter(Bonmin::BabSetupBase::MaxFailures), 
				   s.getIntParameter(Bonmin::BabSetupBase::MaxInfeasible), 
				   s.getIntParameter(Bonmin::BabSetupBase::FailureBehavior));
      model_.setStrategy(strat);
    }
  }

  model_.setMaximumCutPasses(s.getIntParameter(Bonmin::BabSetupBase::NumCutPasses));
  model_.setMaximumCutPassesAtRoot(s.getIntParameter(Bonmin::BabSetupBase::NumCutPassesAtRoot));

  //Setup cutting plane methods
  for (Bonmin::BabSetupBase::CuttingMethods::iterator i = s.cutGenerators().begin() ;
       i != s.cutGenerators().end() ; i++) {

    Bonmin::OaDecompositionBase * oa = dynamic_cast<Bonmin::OaDecompositionBase *>(i->cgl);
    if (oa && oa->reassignLpsolver())
      oa->assignLpInterface(model_.solver());
    model_.addCutGenerator(i->cgl,i->frequency,i->id.c_str(), i->normal,
			   i->atSolution);
    if(i->always){
      model_.cutGenerators()[model_.numberCutGenerators()-1]
	->setMustCallAgain(true);
    }
  }

  for (Bonmin::BabSetupBase::HeuristicMethods::iterator i = s.heuristics().begin() ;
       i != s.heuristics().end() ; i++) {
    CbcHeuristic * heu = i->heuristic;
    heu->setModel(&model_);
    model_.addHeuristic(heu, i->id.c_str());
  }

  //need to record solver logLevel here
  int logLevel = s.continuousSolver()->messageHandler()->logLevel();

  //Set true branch-and-bound parameters
  model_.setLogLevel(s.getIntParameter(Bonmin::BabSetupBase::BabLogLevel));

  // Put back solver logLevel
  model_.solver()->messageHandler()->setLogLevel(logLevel);

  model_.setPrintFrequency(s.getIntParameter(Bonmin::BabSetupBase::BabLogInterval));

  bool ChangedObject = false;
  //Pass over user set branching priorities to Cbc
  if (s.continuousSolver()->objects()==NULL) {
    //assert (s.branchingMethod() == NULL);
    const OsiTMINLPInterface * nlpSolver = s.nonlinearSolver();
    //set priorities, prefered directions...
    const int * priorities = nlpSolver->getPriorities();
    const double * upPsCosts = nlpSolver->getUpPsCosts();
    const double * downPsCosts = nlpSolver->getDownPsCosts();
    const int * directions = nlpSolver->getBranchingDirections();
    bool hasPseudo = (upPsCosts!=NULL);
    model_.findIntegers(true,hasPseudo);
    OsiObject ** simpleIntegerObjects = model_.objects();
    int numberObjects = model_.numberObjects();
    if (priorities != NULL || directions != NULL || hasPseudo) {
      ChangedObject = true;
      for (int i = 0 ; i < numberObjects ; i++) {
	CbcObject * object = dynamic_cast<CbcObject *>
	  (simpleIntegerObjects[i]);
	int iCol = object->columnNumber();
	if (priorities)
	  object->setPriority(priorities[iCol]);
	if (directions)
	  object->setPreferredWay(directions[iCol]);
	if (upPsCosts) {
	  CbcSimpleIntegerPseudoCost * pscObject =
	    dynamic_cast<CbcSimpleIntegerPseudoCost*> (object);
	  pscObject->setUpPseudoCost(upPsCosts[iCol]);
	  pscObject->setDownPseudoCost(downPsCosts[iCol]);
	}
      }
    }

#if 1
    // Now pass user set Sos constraints (code inspired from CoinSolve.cpp)
    const TMINLP::SosInfo * sos = s.nonlinearSolver()->model()->sosConstraints();

    if (!s.getIntParameter(Bonmin::BabSetupBase::DisableSos) && sos && sos->num > 0) {

      // we have some sos constraints

        const OsiTMINLPInterface * nlpSolver = s.nonlinearSolver();
        const int & numSos = sos->num;
	(*nlpSolver->messageHandler())<<"Adding "<<sos->num<<" sos constraints."
				      <<CoinMessageEol;

        CbcObject ** objects = new CbcObject*[numSos];
        const int * starts = sos->starts;
        const int * indices = sos->indices;
        const char * types = sos->types;
        const double * weights = sos->weights;
        //verify if model has user set priorities
        bool hasPriorities = false;
        const int * varPriorities = nlpSolver->getPriorities();
        int numberObjects = model_.numberObjects();
        if (varPriorities)
	  {
	    for (int i = 0 ; i < numberObjects ; i++) {
	      if (varPriorities[i]) {
		hasPriorities = true;
		break;
	      }
	    }
	  }
        const int * sosPriorities = sos->priorities;
        if (sosPriorities)
	  {
	    for (int i = 0 ; i < numSos ; i++) {
	      if (sosPriorities[i]) {
		hasPriorities = true;
		break;
	      }
	    }
	  }
        for (int i = 0 ; i < numSos ; i++)
	  {
	    int start = starts[i];
	    int length = starts[i + 1] - start;
#ifdef DO_IT_NWAY
	    printf("setting nway object\n"),
	      objects[i] = new CbcNWay(&model_, length, &indices[start],
				       i);
	    objects[i]->setPriority(1);
#else
	    objects[i] = new CbcSOS(&model_, length, &indices[start],
				    &weights[start], i, types[i]);
	    objects[i]->setPriority(10);
#endif
	    if (hasPriorities && sosPriorities && sosPriorities[i]) {
	      objects[i]->setPriority(sosPriorities[i]);
	    }
	  }
        model_.addObjects (numSos, objects);
        for (int i = 0 ; i < numSos ; i++)
          delete objects[i];
        delete [] objects;
      }
#endif
    //If Setup contains more objects add them to Cbc
    if (s.objects().size()) {
      CbcObject ** objects = new CbcObject *[s.objects().size()];
      for (unsigned int i = 0 ; i < s.objects().size() ; i++) {
	objects[i] = dynamic_cast<CbcObject *> (s.objects()[i]);
	assert(objects[i]);
	objects[i]->setModel(&model_);
      }
      model_.addObjects ((int) s.objects().size(), objects);
      delete [] objects;
    }

    replaceIntegers(model_.objects(), model_.numberObjects());

  } else { // Pass in objects to Cbc

    // Redundant definition of default branching (as Default == User)
    assert (s.branchingMethod() != NULL);

    // Add nonlinear and integer objects (need to add OsiSOS)
    model_.addObjects (s.continuousSolver () -> numberObjects (), s.continuousSolver () -> objects ());

    // Now model_ has only CouenneObjects and SOS objects

    // for (int i=0; i<nco; i++) 
    //   if (!(dynamic_cast <CbcSimpleInteger *> (s.continuousSolver () -> objects () [i])))
    // 	model_ . objects () [nRealObj++] = s.continuousSolver () -> objects () [i] -> clone ();

    CbcBranchDefaultDecision branch;
    s.branchingMethod()->setSolver(model_.solver());
    BonChooseVariable * strong2 = dynamic_cast<BonChooseVariable *>(s.branchingMethod());
    if (strong2)
      strong2->setCbcModel(&model_);
    branch.setChooseMethod(*s.branchingMethod());

    model_.setBranchingMethod(&branch);
    // prevent duplicating object when copying in CbcModel.cpp
    model_.solver()->deleteObjects();
  }

  model_.setDblParam(CbcModel::CbcCutoffIncrement, s.getDoubleParameter(Bonmin::BabSetupBase::CutoffDecr));

  model_.setCutoff(s.getDoubleParameter(Bonmin::BabSetupBase::Cutoff) + CUTOFF_TOL);

  model_.setDblParam(CbcModel::CbcAllowableGap, s.getDoubleParameter(Bonmin::BabSetupBase::AllowableGap));
  model_.setDblParam(CbcModel::CbcAllowableFractionGap, s.getDoubleParameter(Bonmin::BabSetupBase::AllowableFractionGap));

  // Definition of node selection strategy

  if (s.nodeComparisonMethod()==Bonmin::BabSetupBase::bestBound) {
    CbcCompareObjective compare;
    model_.setNodeComparison(compare);
  }
  else if (s.nodeComparisonMethod()==Bonmin::BabSetupBase::DFS) {
    CbcCompareDepth compare;
    model_.setNodeComparison(compare);
  }
  else if (s.nodeComparisonMethod()==Bonmin::BabSetupBase::BFS) {
    CbcCompareDefault compare;
    compare.setWeight(0.0);
    model_.setNodeComparison(compare);
  }
  else if (s.nodeComparisonMethod()==Bonmin::BabSetupBase::dynamic) {
    CbcCompareDefault compare;
    model_.setNodeComparison(compare);
  }
  else if (s.nodeComparisonMethod()==Bonmin::BabSetupBase::bestGuess) {
    // Right now, this is a mess.  We need a separation of the
    // pseudo costs from the ChooseVariable method
    CbcCompareEstimate compare;
    model_.setNodeComparison(compare);
    GuessHeuristic * guessHeu = new GuessHeuristic(model_);
    model_.addHeuristic(guessHeu);
    delete guessHeu;
  }

  if (s.treeTraversalMethod() == Bonmin::BabSetupBase::HeapOnly) {
    //Do nothing this is the default of Cbc.
  }
  else if (s.treeTraversalMethod() == Bonmin::BabSetupBase::DiveFromBest) {
    CbcDiver treeTraversal;
    treeTraversal.initialize(s);
    model_.passInTreeHandler(treeTraversal);
  }
  else if (s.treeTraversalMethod() == Bonmin::BabSetupBase::ProbedDive) {
    CbcProbedDiver treeTraversal;
    treeTraversal.initialize(s);
    model_.passInTreeHandler(treeTraversal);
  }
  else if (s.treeTraversalMethod() == Bonmin::BabSetupBase::DfsDiveFromBest) {
    CbcDfsDiver treeTraversal;
    treeTraversal.initialize(s);
    model_.passInTreeHandler(treeTraversal);
  }
  else if (s.treeTraversalMethod() == Bonmin::BabSetupBase::DfsDiveDynamic) {
    CbcDfsDiver treeTraversal;
    treeTraversal.initialize(s);
    model_.passInTreeHandler(treeTraversal);

    DiverCompare compare;
    compare.setComparisonDive(*model_.nodeComparison());
    compare.setComparisonBound(CbcCompareObjective());
    CbcDfsDiver * dfs = dynamic_cast<CbcDfsDiver *> (model_.tree());
    assert(dfs);
    compare.setDiver(dfs);
    model_.setNodeComparison(compare);
  }

  model_.setNumberStrong(s.getIntParameter(Bonmin::BabSetupBase::NumberStrong));
  model_.setNumberBeforeTrust(s.getIntParameter(Bonmin::BabSetupBase::MinReliability));
  model_.setNumberPenalties(8);

  model_.setDblParam(CbcModel::CbcMaximumSeconds, s.getDoubleParameter(Bonmin::BabSetupBase::MaxTime));

  model_.setMaximumNodes(s.getIntParameter(Bonmin::BabSetupBase::MaxNodes));

  model_.setMaximumNumberIterations(s.getIntParameter(Bonmin::BabSetupBase::MaxIterations));

  model_.setMaximumSolutions(s.getIntParameter(Bonmin::BabSetupBase::MaxSolutions));

  model_.setIntegerTolerance(s.getDoubleParameter(Bonmin::BabSetupBase::IntTol));

  //Get objects from model_ if it is not null means there are some sos constraints or non-integer branching object
  // pass them to cut generators.
  OsiObject ** objects = model_.objects();

  if (specOpt!=16 && objects) {

    int numberObjects = model_.numberObjects();
    if (objects_ != NULL) {
      for (int i = 0 ; i < nObjects_; i++)
  	delete objects_[i];
    }
    delete [] objects_;
    objects_ = new OsiObject*[numberObjects];
    nObjects_ = numberObjects;
    for (int i = 0 ; i < numberObjects; i++) {
      OsiObject * obj = objects[i];
      CbcSimpleInteger * intObj = dynamic_cast<CbcSimpleInteger *> (obj);
      if (intObj) {
  	objects_[i] = intObj->osiObject();
      }
      else {
  	CbcSOS * sosObj = dynamic_cast<CbcSOS *>(obj);
  	if (sosObj) objects_[i] = sosObj->osiObject(model_.solver());
  	else {//Maybe an unsupported CbcObject
  	  CbcObject * cbcObj = dynamic_cast<CbcObject *>(obj);
  	  if (cbcObj) {
  	    std::cerr<<"Unsupported CbcObject appears in the code"<<std::endl;
  	    throw UNSUPPORTED_CBC_OBJECT;
  	  }
  	  else {//It has to be an OsiObject.
  	    objects_[i]=obj->clone();
  	  }
  	}
      }
    }
    CbcCutGenerator ** gen = model_.cutGenerators();
    int numGen = model_.numberCutGenerators();
    for (int i = 0 ; i < numGen ; i++) {
      Bonmin::OaDecompositionBase * oa = dynamic_cast<Bonmin::OaDecompositionBase * >(gen[i]->generator());
      // if (oa)
      // 	printf ("\n\n\nat least one OADecompBase\n\n\n");
      if (oa) // pass objects
  	oa->setObjects(objects_,nObjects_);
    }
  }

  // if (objects_) {

  //   for (int i = 0 ; i < nObjects_; i++)
  //     delete objects_ [i];

  //   delete [] objects_;
  // }

  // OsiObject ** objects = model_.objects();
  // int numObjects = model_.numberObjects();

  // nObjects_ = 0;
  // objects_ = new OsiObject* [numObjects];

  // for (int i=0; i < numObjects; ++i)
  //   if (objects [i])
  //     objects_ [nObjects_++] = objects [i] -> clone ();

  try {

    //Get the time and start.

    {
      OsiTMINLPInterface * tmpOsi = NULL;
      if(s.nonlinearSolver() == s.continuousSolver()){
        tmpOsi = dynamic_cast<OsiTMINLPInterface *> (model_.solver());
        tmpOsi->forceSolverOutput(s.getIntParameter(Bonmin::BabSetupBase::RootLogLevel)); 
      }

      model_.initialSolve();

      if(tmpOsi != NULL){
        tmpOsi->setSolverOutputToDefault(); 
      }
    }

    int ival;

    s.options()->GetEnumValue("enable_dynamic_nlp", ival, "bonmin.");

    if(s.nonlinearSolver() == s.continuousSolver() && ival) {

	if(!model_.solver()->isProvenOptimal() ){//Something went wrong check if objective is linear and alternate model
	  // can be solved
	  OsiTMINLPInterface * tmpOsi = dynamic_cast<OsiTMINLPInterface *> (model_.solver());
	  TMINLPLinObj * tmp_tminlp = dynamic_cast<TMINLPLinObj *> (tmpOsi->model());
	  tmpOsi->setModel(tmp_tminlp->tminlp());
	  model_.initialSolve();
	} 
	else {
	  LinearCutsGenerator cgl;
	  cgl.initialize(s); 
	  OsiCuts cuts;
	  cgl.generateCuts(*model_.solver(), cuts);
	  std::vector<const OsiRowCut *> mycuts(cuts.sizeRowCuts());
	  for(int i = 0 ; i < cuts.sizeRowCuts() ; i++){
	    mycuts[i] = cuts.rowCutPtr(i);
	  }
	  model_. solver () -> applyRowCuts ((int) mycuts.size(), (const OsiRowCut **) &mycuts[0]);
	}

	//Added by Claudia
	OsiTMINLPInterface * nlpSolver = dynamic_cast<OsiTMINLPInterface *>(model_.solver());
	if(nlpSolver && nlpSolver->getNewCutoffDecr()!=COIN_DBL_MAX)
          model_.setDblParam(CbcModel::CbcCutoffIncrement, nlpSolver->getNewCutoffDecr());

	model_.solver()->resolve();
      }

    // for Couenne
    model_.passInSolverCharacteristics (bonBabInfoPtr);

    continuousRelaxation_ =model_.solver()->getObjValue();
    if (specOpt==16)//Set warm start point for Ipopt
      {
#if 1
	const double * colsol = model_.solver()->getColSolution();
	const double * duals = model_.solver()->getRowPrice();

	OsiTMINLPInterface * tnlpSolver = dynamic_cast<OsiTMINLPInterface *>(model_.solver());
	// Primal dual point is not copied if one (supposedly a better one) has already been put into the solver.
	if(tnlpSolver->problem()->has_x_init() != 2){
	  model_.solver()->setColSolution(colsol);
	  model_.solver()->setRowPrice(duals);
	}
#else
	OsiTMINLPInterface * tnlpSolver = dynamic_cast<OsiTMINLPInterface *>(model_.solver());
	CoinWarmStart * warm = tnlpSolver->solver()->getWarmStart(tnlpSolver->problem());
	tnlpSolver->solver()->setWarmStart(warm, tnlpSolver->problem());
	delete warm;
#endif

#if 0 // Sometimes primal dual point is problematic in the context of Cut-and-branch
	model_.solver()->resolve();
	if(!model_.solver()->isProvenOptimal())
	  model_.solver()->setColSolution(NULL);
#endif 
      }

#ifdef SIGNAL
    CoinSighandler_t saveSignal = SIG_DFL;
    // signal handler
    saveSignal = signal (SIGINT,couenne_signal_handler);
    currentBranchModel = &model_;
#endif


    // to get node parent info in Cbc, pass parameter 3.
    //model_.branchAndBound(3);
    remaining_time -= CoinCpuTime();
    model_.setDblParam(CbcModel::CbcMaximumSeconds, remaining_time);
    if(remaining_time > 0.)
      model_.branchAndBound();
  }

  catch(TNLPSolver::UnsolvedError *E){
    s.nonlinearSolver()->model()->finalize_solution
      (TMINLP::MINLP_ERROR, 0, NULL, DBL_MAX);
    throw E;
  }

  numNodes_ = model_.getNodeCount();
  bestObj_ = model_.getObjValue();
  bestBound_ = model_.getBestPossibleObjValue();
  mipIterationCount_ = model_.getIterationCount();

  bool hasFailed = false;
  if (specOpt==16)//Did we continue branching on a failure
    {
      CbcNlpStrategy * nlpStrategy = dynamic_cast<CbcNlpStrategy *>(model_.strategy());
      if (nlpStrategy)
        hasFailed = nlpStrategy->hasFailed();
      else
        throw -1;
    }
  else
    hasFailed = s.nonlinearSolver()->hasContinuedOnAFailure();

  // Output summarizing cut generators (taken from CbcSolver.cpp)
  // ToDo put into proper print level

  int numberGenerators = model_.numberCutGenerators();
  for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CbcCutGenerator * generator = model_.cutGenerator(iGenerator);
    //CglStored * stored = dynamic_cast<CglStored*>(generator->generator());
    if (true&&!(generator->numberCutsInTotal() || generator->numberColumnCuts()))
      continue;
    if(modelHandler_->logLevel() >= 1) {
      *modelHandler_ << generator->cutGeneratorName()
		     << "was tried" << generator->numberTimesEntered()
		     << "times and created" << generator->numberCutsInTotal()+generator->numberColumnCuts()
		     << "cuts of which" << generator->numberCutsActive()
		     << "were active after adding rounds of cuts";
      // if (generator->timing()) {
      // 	char timebuf[20];
      // 	sprintf(timebuf, "(%.3fs)", generator->timeInCutGenerator());
      // 	*modelHandler_ << timebuf << CoinMessageEol;
      // }
      // else {
      // 	*modelHandler_ << CoinMessageEol;
      // }
    }
  }

  TMINLP::SolverReturn status = TMINLP::MINLP_ERROR;

  if (model_.numberObjects()==0) {
    if (bestSolution_)
      delete [] bestSolution_;
    OsiSolverInterface * solver = 
      (s.nonlinearSolver() == s.continuousSolver())? 
      model_.solver() : s.nonlinearSolver();
    bestSolution_ = new double[solver->getNumCols()];
    CoinCopyN(solver->getColSolution(), solver->getNumCols(),
	      bestSolution_);
    bestObj_ = bestBound_ = solver->getObjValue();
  }

  if (bonBabInfoPtr->bestSolution2().size() > 0) {
    assert((int) bonBabInfoPtr->bestSolution2().size() == s.nonlinearSolver()->getNumCols());
    if (bestSolution_)
      delete [] bestSolution_;
    bestSolution_ = new double[s.nonlinearSolver()->getNumCols()];
    std::copy(bonBabInfoPtr->bestSolution2().begin(), bonBabInfoPtr->bestSolution2().end(),
	      bestSolution_);
    bestObj_ = (bonBabInfoPtr->bestObj2());
    (*s.nonlinearSolver()->messageHandler())<<"\nReal objective function: "
                                            <<bestObj_<<CoinMessageEol;
  }
  else if (model_.bestSolution()) {
    if (bestSolution_)
      delete [] bestSolution_;
    bestSolution_ = new double[s.nonlinearSolver()->getNumCols()];
    CoinCopyN(model_.bestSolution(), s.nonlinearSolver()->getNumCols(), bestSolution_);
  }
  if(remaining_time <= 0.){
    status = TMINLP::LIMIT_EXCEEDED;
    if (bestSolution_) {
      mipStatus_ = Feasible;
    }
    else {
      mipStatus_ = NoSolutionKnown;
    }
  }
  else if (model_.status() == 0) {
    if(model_.isContinuousUnbounded()){
      status = TMINLP::CONTINUOUS_UNBOUNDED;
      mipStatus_ = UnboundedOrInfeasible;
    }
    else
      if (bestSolution_) {
        status = TMINLP::SUCCESS;
        mipStatus_ = FeasibleOptimal;
      }
      else {
        status = TMINLP::INFEASIBLE;
        mipStatus_ = ProvenInfeasible;
      }
  }
  else if (model_.status() == 1 || model_.status() == 5) {
#if (BONMIN_VERSION_MAJOR > 1) || (BONMIN_VERSION_MINOR > 6)
    status = model_.status() == 1 ? TMINLP::LIMIT_EXCEEDED : TMINLP::USER_INTERRUPT;
#else
    status = TMINLP::LIMIT_EXCEEDED;
#endif
    if (bestSolution_) {
      mipStatus_ = Feasible;
    }
    else {
      mipStatus_ = NoSolutionKnown;
    }
  }
  else if (model_.status()==2) {
    status = TMINLP::MINLP_ERROR;
  }

  // Which solution should we use? false if RBS's, true if Cbc's
  bool use_RBS_Cbc = 
    !problem_ ||
    !(problem_ -> getRecordBestSol ()) ||
    !(problem_ -> getRecordBestSol () -> getHasSol()) ||
    (((fabs (bestObj_) < COUENNE_INFINITY / 1e4) && 
      (problem_ -> getRecordBestSol () -> getVal () > bestObj_)));

  /* if we do not pass the cbc solution and problem_ -> getRecordBestSol () -> getHasSol() is true, then there should be a solution vector in problem_ -> getRecordBestSol () */
  assert(use_RBS_Cbc || problem_ -> getRecordBestSol () -> getSol() != NULL);

  s.nonlinearSolver () -> model () -> finalize_solution 
    (status,
     s.nonlinearSolver () -> getNumCols (),
     use_RBS_Cbc ? bestSolution_ : problem_ -> getRecordBestSol () -> getSol (),
     use_RBS_Cbc ? bestObj_      : problem_ -> getRecordBestSol () -> getVal ());
}

const double * CouenneBab::bestSolution() const {
  if(!problem_ ||
     !(problem_ -> getRecordBestSol ()) ||
     !(problem_ -> getRecordBestSol () -> getHasSol()) ||
     (((fabs (bestObj_) < COUENNE_INFINITY / 1e4) &&
       (problem_ -> getRecordBestSol () -> getVal () > bestObj_))))
    return bestSolution_;
  return problem_ -> getRecordBestSol () -> getSol ();
}

double CouenneBab::bestObj() const {
  if(!problem_ ||
     !(problem_ -> getRecordBestSol ()) ||
     !(problem_ -> getRecordBestSol () -> getHasSol()) ||
     (((fabs (bestObj_) < COUENNE_INFINITY / 1e4) &&
       (problem_ -> getRecordBestSol () -> getVal () > bestObj_))))
    return bestObj_;
  return problem_ -> getRecordBestSol () -> getVal ();
}
