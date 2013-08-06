/* $Id$
 *
 * Name:    CouenneNlpRoundOne.hpp
 * Author:  Pietro Belotti
 * Purpose: Implementation of a rounding heuristic that simply calls a
 *          NLP solver on a (possibly fractional) solution to obtain
 *          an integer feasible one
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneExprVar.hpp"

#include "BonCouenneInterface.hpp"
#include "CouenneProblem.hpp"
#include "BonAuxInfos.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "CouenneRecordBestSol.hpp"

#include "CouenneNlpRoundOne.hpp"

using namespace Ipopt;
using namespace Couenne;

CouenneNlpRoundOne::CouenneNlpRoundOne ():

  CbcHeuristic (),

  nlp_                 (NULL),
  hasCloned_           (false),
  couenne_             (NULL),
  numberSolvePerLevel_ (-1) {

  setHeuristicName     ("ClassifierNlp");
}
  
CouenneNlpRoundOne::CouenneNlpRoundOne (CbcModel &model, Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp, CouenneProblem *couenne):

  CbcHeuristic (model), 

  nlp_                 (&nlp), 
  hasCloned_           (cloneNlp),
  couenne_             (couenne),
  numberSolvePerLevel_ (-1) {

  setHeuristicName ("ClassifierNlp");

  if (cloneNlp)
    nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (nlp.clone ());
}
  
CouenneNlpRoundOne::CouenneNlpRoundOne (const CouenneNlpRoundOne &other):

  CbcHeuristic         (other), 

  nlp_                 (other.nlp_),
  hasCloned_           (other.hasCloned_),
  couenne_             (other.couenne_), 
  numberSolvePerLevel_ (other.numberSolvePerLevel_) {

  if (hasCloned_ && nlp_ != NULL)
    nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (other.nlp_ -> clone());
}
 
CbcHeuristic *CouenneNlpRoundOne::clone () const {
  return new CouenneNlpRoundOne (*this);
}
  
CouenneNlpRoundOne &CouenneNlpRoundOne::operator= (const CouenneNlpRoundOne & rhs) {

  if (this != &rhs){

    CbcHeuristic::operator= (rhs);

    if (hasCloned_ && nlp_)
      delete nlp_;

    hasCloned_ = rhs.hasCloned_;
    numberSolvePerLevel_ = rhs. numberSolvePerLevel_;

    if (nlp_ != NULL) {

      if(hasCloned_) nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (rhs.nlp_ -> clone ());
      else           nlp_ = rhs.nlp_;
    }
  }

  couenne_ = rhs.couenne_;
  return *this;
}
  
CouenneNlpRoundOne::~CouenneNlpRoundOne(){
  if (hasCloned_)
    delete nlp_;
  nlp_ = NULL;
}
  
void CouenneNlpRoundOne::setNlp (Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp) {

  if (hasCloned_ && nlp_ != NULL)
    delete nlp_;

  hasCloned_ = cloneNlp;

  if (cloneNlp) nlp_ = dynamic_cast <Bonmin::OsiTMINLPInterface *> (nlp.clone ());
  else          nlp_ = &nlp;
}
  
void CouenneNlpRoundOne::setCouenneProblem (CouenneProblem * couenne)
{couenne_ = couenne;}


int CouenneNlpRoundOne::solution (double &objectiveValue, double *newSolution) {

  // Simply run an NLP solver starting from the current solution. The
  // only effect should be that all binary variables that are still
  // fractional should be set to one

  int noSolution = 1, maxTime = 2;

  // Do heuristic the usual way, but if for any reason (time is up, no
  // better solution found) there is no improvement, get the best
  // solution from the GlobalCutOff object in the pointer to the
  // CouenneProblem and return it instead.
  //
  // Although this should be handled by Cbc, very often this doesn't
  // happen.

  const int depth = (model_ -> currentNode ()) ? model_ -> currentNode () -> depth () : 0;

  if (depth <= 0)
    couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "Classifier heuristic: "); fflush (stdout);

  try {

    if (CoinCpuTime () > couenne_ -> getMaxCpuTime ())
      throw maxTime;

    OsiSolverInterface * solver = model_ -> solver();

    OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
    Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (auxInfo);

    if (babInfo){
      babInfo -> setHasNlpSolution (false);
      if (babInfo->infeasibleNode()){
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
	  (1., (double) ((depth - numberSolvePerLevel_) * 
			 (depth - numberSolvePerLevel_))))
	too_deep = true;
    }

    if (too_deep)
      throw maxTime;

    double *lower = new double [couenne_ -> nVars ()];
    double *upper = new double [couenne_ -> nVars ()];

    CoinFillN (lower, couenne_ -> nVars (), -COUENNE_INFINITY);
    CoinFillN (upper, couenne_ -> nVars (),  COUENNE_INFINITY);

    CoinCopyN (solver -> getColLower (), nlp_ -> getNumCols (), lower);
    CoinCopyN (solver -> getColUpper (), nlp_ -> getNumCols (), upper);

    const double * solution = solver -> getColSolution();
    OsiBranchingInformation info (solver, true);

    double *Y = new double [couenne_ -> nVars ()];
    CoinFillN (Y, couenne_ -> nVars (), 0.);
    CoinCopyN (solution, nlp_ -> getNumCols (), Y);

    bool foundSolution = false;

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

    /*
    printf ("---------------------------------------------------------------------------------------------------\nInitial solution:  ");
    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("%c", 
		Y [i] <= 1e-7   ? '0' : 
		Y [i] <= 1e-4   ? '.' :
		Y [i] <= 1-1e-4 ? '-' : 
		Y [i] <= 1-1e-7 ? '^' : '1');
    printf ("\n");
    */

    // apply NLP solver /////////////////////////////////
    try {
      nlp_ -> options () -> SetNumericValue ("max_cpu_time", CoinMax (0., couenne_ -> getMaxCpuTime () - CoinCpuTime ()));
      //nlp_ -> messageHandler () -> setLogLevel (8);
      nlp_ -> initialSolve ();
    }

    catch (Bonmin::TNLPSolver::UnsolvedError *E) {}

    /*
    const double *sol = nlp_ -> getColSolution ();

    printf ("-- Solution change:");
    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("%c", 
		(Y [i] - sol [i] <= -1e-4)   ? '^' : 
		(Y [i] - sol [i] >=  1e-4)   ? 'V' : ' ');
    printf ("\n");

    printf ("-- Solution:       ");
    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("%c", 
		(sol   [i]       <=1e-4)   ? '0' : 
		(sol   [i]       >=1-1e-4) ? '1' : '-');
    printf ("\n");

    printf ("-- lower    bounds:");
    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("%c", 
		(lower [i]       <=1e-4)   ? '0' : 
		(lower [i]       >=1-1e-4) ? '1' : '-');
    printf ("\n");

    printf ("-- upper    bounds:");
    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("%c", 
		(upper [i]       <=  1e-4) ? '0' : 
		(upper [i]       >=1-1e-4) ? '1' : '-');
    printf ("\n----------------------------------------------------------------------------------------\n");
    */

    double obj = (nlp_ -> isProvenOptimal()) ? nlp_ -> getObjValue (): COIN_DBL_MAX;

    for (int i=0; i<couenne_ -> nVars (); ++i) 
      if (couenne_ -> Var (i) -> isInteger ())
	if (!(::isInteger (Y [i])))
	  Y [i] = ceil (Y [i]);

    if (nlp_ -> isProvenOptimal () &&
	couenne_ -> checkNLP (nlp_ -> getColSolution (), obj, true) && // true for recomputing obj
	(obj < couenne_ -> getCutOff ())) {

      //printf ("====================================== improved solution: %g\n", obj);

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
	printf("CouenneNlpRoundOne::solution(): ### ERROR: checkNLP(): returns true,  checkNLP2() returns false\n");
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

    nlp_ -> setColLower (saveColLower);
    nlp_ -> setColUpper (saveColUpper);

    delete [] saveColLower;
    delete [] saveColUpper;
    //}

    delete [] Y;

    delete [] lower;
    delete [] upper;

    if (couenne_ -> Jnlst () -> ProduceOutput (J_ERROR, J_COUENNE) && (depth <= 0)) {

      if (foundSolution) printf ("Solution found, obj. %g\n", objectiveValue);
      else               printf ("No solution.\n");
    }

    return foundSolution;
  }

  catch (int &e) {

    // forget about using the global cutoff. That has to trickle up to
    // Cbc some other way

    bool output = couenne_ -> Jnlst () -> ProduceOutput (J_ERROR, J_COUENNE) && (depth <= 0);

    if      (e==noSolution) {if (output) printf ("No solution.\n");                            return 0;}
    else if (e==maxTime)    {if (output) printf ("Time limit reached.\n");                     return 0;}
    else                    {if (output) printf ("Solution found, obj. %g\n", objectiveValue); return 1;}

    // // no solution available? Use the one from the global cutoff
    // if ((couenne_ -> getCutOff () < objectiveValue) &&
    // 	couenne_ -> getCutOffSol ()) {
    //   objectiveValue = couenne_ -> getCutOff    ();
    //   CoinCopyN       (couenne_ -> getCutOffSol (), couenne_ -> nVars (), newSolution);
    //   if (depth <= 0)
    // 	couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "solution found, obj. %g\n", objectiveValue);
    //   return 1;
    // } else {
    //   if (depth <= 0 && e==noSolution)
    // 	couenne_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "no solution.\n", objectiveValue);
    //   return 0;
    // }
  }
}


/// initialize options
void CouenneNlpRoundOne::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2
    ("classifier_heuristic",
     "Search for local solutions of classification MINLPs by runnning an NLP",
     "yes",
     "no","",
     "yes","",
     "If enabled, a heuristic based on Ipopt is used to find feasible solutions for the problem. "
     "It is highly recommended that this option is left enabled, as it would be difficult to find feasible solutions otherwise.");

  roptions -> AddLowerBoundedIntegerOption
    ("heur_classifier_level",
     "Specify the logarithm of the number of classifier heuristic calls to perform" 
     " on average for each level of given depth of the tree.",
     -1,
     2, "Solve as many nlp's at the nodes for each level of the tree. "
     "Nodes are randomly selected. If for a "
     "given level there are less nodes than this number nlp are solved for every nodes. "
     "For example if parameter is 8, nlp's are solved for all node until level 8, " 
     "then for half the node at level 9, 1/4 at level 10.... "
     "Value -1 specify to perform at all nodes.");
}
