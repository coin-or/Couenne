/* $Id$ */
// (C) Copyright International Business Machines Corporation 2007 
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
// Pietro Belotti, Lehigh University
//
// Date : 04/09/2007

#include "CouenneCutGenerator.hpp"

#include "BonCouenneInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"
#include "CbcCutGenerator.hpp"
//#include "CbcBranchActual.hpp"
#include "BonAuxInfos.hpp"
#include "CoinHelperFunctions.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonNlpHeuristic.hpp"
#include "CouenneRecordBestSol.hpp"

using namespace Ipopt;
using namespace Couenne;

NlpSolveHeuristic::NlpSolveHeuristic():
  CbcHeuristic(),
  nlp_(NULL),
  hasCloned_(false),
  maxNlpInf_(maxNlpInf_0),
  numberSolvePerLevel_(-1),
  couenne_(NULL){
  setHeuristicName("NlpSolveHeuristic");
}
  
NlpSolveHeuristic::NlpSolveHeuristic(CbcModel & model, Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp, CouenneProblem * couenne):
  CbcHeuristic(model), nlp_(&nlp), hasCloned_(cloneNlp),maxNlpInf_(maxNlpInf_0),
  numberSolvePerLevel_(-1),
  couenne_(couenne){
  setHeuristicName("NlpSolveHeuristic");
  if(cloneNlp)
    nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (nlp.clone());
  }
  
NlpSolveHeuristic::NlpSolveHeuristic(const NlpSolveHeuristic & other):
  CbcHeuristic(other), nlp_(other.nlp_), 
  hasCloned_(other.hasCloned_),
  maxNlpInf_(other.maxNlpInf_),
  numberSolvePerLevel_(other.numberSolvePerLevel_),
  couenne_(other.couenne_){
  if(hasCloned_ && nlp_ != NULL)
    nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (other.nlp_->clone());
}
  
CbcHeuristic * 
NlpSolveHeuristic::clone() const{
  return new NlpSolveHeuristic(*this);
}
  
NlpSolveHeuristic &
NlpSolveHeuristic::operator=(const NlpSolveHeuristic & rhs){
  if(this != &rhs){
    CbcHeuristic::operator=(rhs);
    if(hasCloned_ && nlp_)
      delete nlp_;
      
    hasCloned_ = rhs.hasCloned_;
    if(nlp_ != NULL){
      if(hasCloned_)
	nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (rhs.nlp_->clone());
      else
	nlp_ = rhs.nlp_;
    }
  }
  maxNlpInf_ = rhs.maxNlpInf_;
  numberSolvePerLevel_ = rhs.numberSolvePerLevel_;
  couenne_ = rhs.couenne_;
  return *this;
}
  
NlpSolveHeuristic::~NlpSolveHeuristic(){
  if(hasCloned_)
    delete nlp_;
  nlp_ = NULL;
}
  
void
NlpSolveHeuristic::setNlp (Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp){
  if(hasCloned_ && nlp_ != NULL)
    delete nlp_;
  hasCloned_ = cloneNlp;
  if(cloneNlp)
    nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (nlp.clone());
  else
    nlp_ = &nlp;
}
  
void
NlpSolveHeuristic::setCouenneProblem(CouenneProblem * couenne)
{couenne_ = couenne;}


int
NlpSolveHeuristic::solution (double & objectiveValue, double * newSolution) {

  int noSolution = 1, maxTime = 2;

  // do heuristic the usual way, but if for any reason (time is up, no
  // better solution found) there is no improvement, get the best
  // solution from the GlobalCutOff object in the pointer to the
  // CouenneProblem and return it instead.
  //
  // Although this should be handled by Cbc, very often this doesn't
  // happen.

  //  int nodeDepth = -1;

  const int depth = (model_ -> currentNode ()) ? model_ -> currentNode () -> depth () : 0;

  if (depth <= 0)
    couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "NLP Heuristic: "); fflush (stdout);

  try {

  if (CoinCpuTime () > couenne_ -> getMaxCpuTime ())
    throw maxTime;

  OsiSolverInterface * solver = model_ -> solver();

  OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (auxInfo);

  if(babInfo){
    babInfo->setHasNlpSolution(false);
    if(babInfo->infeasibleNode()){
      throw noSolution;
    }
  }

  // if too deep in the BB tree, only run NLP heuristic if
  // feasibility is low
  bool too_deep = false;

  // check depth
  if (numberSolvePerLevel_ > -1) {

    if (numberSolvePerLevel_ == 0) 
      throw maxTime;

    //if (CoinDrand48 () > pow (2., numberSolvePerLevel_ - depth))
    if (CoinDrand48 () > 1. / CoinMax 
	(1., (double) ((depth - numberSolvePerLevel_) * (depth - numberSolvePerLevel_))))
      too_deep = true;
  }

  if (too_deep)
    throw maxTime;

  double *lower = new double [couenne_ -> nVars ()];
  double *upper = new double [couenne_ -> nVars ()];

  CoinFillN (lower, couenne_ -> nVars (), -COUENNE_INFINITY);
  CoinFillN (upper, couenne_ -> nVars (),  COUENNE_INFINITY);

  CoinCopyN (solver->getColLower(), nlp_ -> getNumCols (), lower);
  CoinCopyN (solver->getColUpper(), nlp_ -> getNumCols (), upper);

  /*printf ("-- int candidate, before: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
    printf ("[%g %g] ", lower [i], upper [i]);
    printf ("\n");*/

  const double * solution = solver->getColSolution();
  OsiBranchingInformation info (solver, true);
  const int & numberObjects = model_->numberObjects();
  OsiObject ** objects = model_->objects();
  double maxInfeasibility = 0;

  bool haveRoundedIntVars = false;

  for (int i = 0 ; i < numberObjects ; i++) {

    CouenneObject * couObj = dynamic_cast <CouenneObject *> (objects [i]);

    if (couObj) {
      if (too_deep) { // only test infeasibility if BB level is high
	int dummy;
	double infeas;
	maxInfeasibility = CoinMax ( maxInfeasibility, infeas = couObj->infeasibility(&info, dummy));

	if (maxInfeasibility > maxNlpInf_){
	  delete [] lower;
	  delete [] upper;
	  throw noSolution;
	}
      }
    } else {

      OsiSimpleInteger * intObj = dynamic_cast<OsiSimpleInteger *>(objects[i]);

      if (intObj) {
	const int & i = intObj -> columnNumber ();
	// Round the variable in the solver
	double value = solution [i];
	if (value < lower[i])
	  value = lower[i];
	else if (value > upper[i])
	  value = upper[i];

	double rounded = floor (value + 0.5);

	if (fabs (value - rounded) > COUENNE_EPS) {
	  haveRoundedIntVars = true;
	  //value = rounded;
	}

	// fix bounds anyway, if a better candidate is not found
	// below at least we have an integer point
	//lower[i] = upper[i] = value;
      }
      else{

	// Probably a SOS object -- do not stop here
	//throw CoinError("Bonmin::NlpSolveHeuristic","solution",
	//"Unknown object.");
      }
    }
  }

  // if here, it means the infeasibility is not too high. Generate a
  // better integer point as there are rounded integer variables

  bool skipOnInfeasibility = false;

  double *Y = new double [couenne_ -> nVars ()];
  CoinFillN (Y, couenne_ -> nVars (), 0.);
  CoinCopyN (solution, nlp_ -> getNumCols (), Y);

  /*printf ("-- int candidate, upon call: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
    if (couenne_ -> Var (i) -> isInteger ())
    printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
    else printf ("%g ", Y [i]);
    printf ("\n");*/

  if (haveRoundedIntVars) // create "good" integer candidate for Ipopt
    skipOnInfeasibility = (couenne_ -> getIntegerCandidate (solution, Y, lower, upper) < 0);

  /*printf ("-- int candidate, after: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
    if (couenne_ -> Var (i) -> isInteger ())
    printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
    else printf ("%g ", Y [i]);
    printf ("\n");*/

  bool foundSolution = false;

  if (haveRoundedIntVars && skipOnInfeasibility) 
    // no integer initial point could be found, make up some random rounding

    for (int i = couenne_ -> nOrigVars (); i--;) 

      if (couenne_ -> Var (i) -> isDefinedInteger ())
	lower [i] = upper [i] = Y [i] = 
	  (CoinDrand48 () < 0.5) ? 
	  floor (Y [i] + COUENNE_EPS) : 
	  ceil  (Y [i] - COUENNE_EPS);

      else if (lower [i] > upper [i]) { 

	// sanity check (should avoid problems in ex1263 with
	// couenne.opt.obbt)

	double swap = lower [i];
	lower [i] = upper [i];
	upper [i] = swap;
      }

  {
    //	printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);

    /*printf ("int candidate: ");
      for (int i=0; i<couenne_ -> nOrig (); i++) 
      if (couenne_ -> Var (i) -> isInteger ())
      printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
      else printf ("%g ", Y [i]);
      printf ("\n");*/

    // Now set column bounds and solve the NLP with starting point
    double * saveColLower = CoinCopyOfArray (nlp_ -> getColLower (), nlp_ -> getNumCols ());
    double * saveColUpper = CoinCopyOfArray (nlp_ -> getColUpper (), nlp_ -> getNumCols ());

    for (int i = nlp_ -> getNumCols (); i--;) {

      if (lower [i] > upper [i]) {
	double swap = lower [i];
	lower [i] = upper [i];
	upper [i] = swap;
      }

      if      (Y [i] < lower [i]) Y [i] = lower [i];
      else if (Y [i] > upper [i]) Y [i] = upper [i];
    }

    nlp_ -> setColLower    (lower);
    nlp_ -> setColUpper    (upper);
    nlp_ -> setColSolution (Y);

    // apply NLP solver /////////////////////////////////
    try {
      nlp_ -> options () -> SetNumericValue ("max_cpu_time", CoinMax (0., couenne_ -> getMaxCpuTime () - CoinCpuTime ()));
      nlp_ -> initialSolve ();
    }
    catch (Bonmin::TNLPSolver::UnsolvedError *E) {}

    double obj = (nlp_ -> isProvenOptimal()) ? nlp_ -> getObjValue (): COIN_DBL_MAX;

    if (nlp_ -> isProvenOptimal () &&
	couenne_ -> checkNLP (nlp_ -> getColSolution (), obj, true) && // true for recomputing obj
	(obj < couenne_ -> getCutOff ())) {

      // store solution in Aux info

      const int nVars = solver->getNumCols();
      double* tmpSolution = new double [nVars];
      CoinCopyN (nlp_ -> getColSolution(), nlp_ -> getNumCols(), tmpSolution);

      //Get correct values for all auxiliary variables
      CouenneInterface * couenne = dynamic_cast <CouenneInterface *> (nlp_);

      if (couenne)
	couenne_ -> getAuxs (tmpSolution);

#ifdef FM_CHECKNLP2
      if(!couenne_->checkNLP2(tmpSolution, 
			      0, false, // do not care about obj
			      true, // stopAtFirstViol 
			      false, // checkAll
			      couenne_->getFeasTol())) {
#ifdef FM_USE_REL_VIOL_CONS
	printf("NlpSolveHeuristic::solution(): ### ERROR: checkNLP(): returns true,  checkNLP2() returns false\n");
	exit(1);
#endif
      }
      obj = couenne_->getRecordBestSol()->getModSolVal(); 
      couenne_->getRecordBestSol()->update();
#else
      couenne_->getRecordBestSol()->update(tmpSolution, nVars,
					   obj, couenne_->getFeasTol());
#endif

      if (babInfo){
	babInfo->setNlpSolution (tmpSolution, nVars, obj);
	babInfo->setHasNlpSolution (true);
      }

      if (obj < objectiveValue) { // found better solution?

	const CouNumber 
	  *lb = solver -> getColLower (),
	  *ub = solver -> getColUpper ();

	// check bounds once more after getAux. This avoids false
	// asserts in CbcModel.cpp:8305 on integerTolerance violated
	for (int i=0; i < nVars; i++, lb++, ub++) {

	  CouNumber &t = tmpSolution [i];
	  if      (t < *lb) t = *lb;
	  else if (t > *ub) t = *ub;
	}
	  
	//printf ("new cutoff %g from BonNlpHeuristic\n", obj);
	couenne_ -> setCutOff (obj);
	foundSolution = true;
	objectiveValue = obj;
	CoinCopyN (tmpSolution, nVars, newSolution);
      }
      delete [] tmpSolution;
    }

    nlp_->setColLower (saveColLower);
    nlp_->setColUpper (saveColUpper);

    delete [] saveColLower;
    delete [] saveColUpper;
  }

  delete [] Y;

  delete [] lower;
  delete [] upper;

  if (depth <= 0) {

    if (foundSolution) couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "solution found, obj. %g\n", objectiveValue);
    else               couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "no solution.\n");
  }

  return foundSolution;

  }
  catch (int &e) {

    // no solution available? Use the one from the global cutoff

    if ((couenne_ -> getCutOff () < objectiveValue) &&
	couenne_ -> getCutOffSol ()) {

      objectiveValue = couenne_ -> getCutOff    ();
      CoinCopyN       (couenne_ -> getCutOffSol (), couenne_ -> nVars (), newSolution);

      if (depth <= 0)
	couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "solution found, obj. %g\n", objectiveValue);

      return 1;

    } else {

      if (depth <= 0 && e==noSolution)
	couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "no solution.\n", objectiveValue);

      return 0;
    }
  }
}


/// initialize options
void NlpSolveHeuristic::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2
    ("local_optimization_heuristic",
     "Search for local solutions of MINLPs",
     "yes",
     "no","",
     "yes","",
     "If enabled, a heuristic based on Ipopt is used to find feasible solutions for the problem. "
     "It is highly recommended that this option is left enabled, as it would be difficult to find feasible solutions otherwise.");

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_local_optimization_per_level",
     "Specify the logarithm of the number of local optimizations to perform" 
     " on average for each level of given depth of the tree.",
     -1,
     2, "Solve as many nlp's at the nodes for each level of the tree. "
     "Nodes are randomly selected. If for a "
     "given level there are less nodes than this number nlp are solved for every nodes. "
     "For example if parameter is 8, nlp's are solved for all node until level 8, " 
     "then for half the node at level 9, 1/4 at level 10.... "
     "Value -1 specify to perform at all nodes.");
}
