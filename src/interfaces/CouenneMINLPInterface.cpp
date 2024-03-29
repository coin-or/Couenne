/*
 *
 * Name:    CouenneMINLPInterface
 * Author:  Pietro Belotti
 * Purpose: MINLP interface with computation of gradient and Jacobian
 *          through CouenneExpressions
 *
 * (C) Pietro Belotti 2010
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneMINLPInterface.hpp"

using namespace Couenne;

/// sets the initial solution for the NLP solver
void CouenneMINLPInterface::setInitSol (const CouNumber *sol) {

}

/// solves and returns the optimal objective function and the
/// solution
CouNumber CouenneMINLPInterface::solve (CouNumber *solution) {

  return 0.;
}


#if 0
//#include "BonminConfig.h"

//#include "BonOsiTMINLPInterface.hpp"
#include "CoinTime.hpp"
#include <climits>
#include <string>
#include <sstream>
#include "BonAuxInfos.hpp"

#include "Ipopt/BonIpoptSolver.hpp"
#include "Ipopt/BonIpoptWarmStart.hpp"

#ifdef COUENNE_HAS_FILTERSQP
#include "Filter/BonFilterSolver.hpp"
#include "Filter/BonFilterWarmStart.hpp"
//#include "Filter/BonBqpdWarmStart.hpp"
#endif

#include "OsiBranchingObject.hpp"
#include "BonStrongBranchingSolver.hpp"

//Macros to try and make messages definition less heavy
#include "BonMsgUtils.hpp"

using namespace Ipopt;
using namespace Couenne;

///Register options
static void
register_general_options
(SmartPtr<RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("nlp interface option", RegisteredOptions::BonminCategory);
  roptions->AddStringOption3("nlp_solver",
                             "Choice of the solver for local optima of continuous nlp's",
                             "Ipopt",
                             "Ipopt", "Interior Point OPTimizer (https://projects.coin-or.org/Ipopt)",
                             "filterSQP", "Sequential quadratic programming trust region "
			     "algorithm (http://www-unix.mcs.anl.gov/~leyffer/solvers.html)",
                             "all", "run all available solvers at each node",
                             "Note that option will work only if the specified solver has been installed. Ipopt will usualy be installed with Bonmin by default. For FilterSQP please see http://www-unix.mcs.anl.gov/~leyffer/solvers.html on how to obtain it and https://projects.coin-or.org/Bonmin/wiki/HintTricks on how to configure Bonmin to use it.");
  roptions->setOptionExtraInfo("nlp_solver",15);
  roptions->AddBoundedIntegerOption("nlp_log_level",
                                    "specify NLP solver interface log level (independent from ipopt print_level).",
				    0,2,1,
                                    "Set the level of output of the CouenneMINLPInterface : "
                                    "0 - none, 1 - normal, 2 - verbose"
				    );
  roptions->setOptionExtraInfo("nlp_log_level",15);

  roptions->AddStringOption2("file_solution",
			     "Write a file bonmin.sol with the solution",
			     "no",
			     "yes","",
			     "no","","");

  roptions->AddStringOption3("warm_start",
			     "Select the warm start method",
			     "none",
			     "none","No warm start",
			     "optimum","Warm start with direct parent optimum",
			     "interior_point","Warm start with an interior point of direct parent",
			     "This will affect the function getWarmStart(), and as a consequence the warm starting in the various algorithms.");
  roptions->setOptionExtraInfo("warm_start",8);

  roptions->SetRegisteringCategory("Nlp solution robustness", RegisteredOptions::BonminCategory);

  roptions->AddLowerBoundedNumberOption("max_random_point_radius",
					"Set max value r for coordinate of a random point.",
					0.,1,1e5,
					"When picking a random point coordinate i will be in the interval [min(max(l,-r),u-r), max(min(u,r),l+r)] "
					"(where l is the lower bound for the variable and u is its upper bound)");
  roptions->setOptionExtraInfo("max_random_point_radius",8);

  roptions->AddStringOption3("random_point_type","method to choose a random starting point",
			     "Jon",
			     "Jon", "Choose random point uniformly between the bounds",
			     "Andreas", "perturb the starting point of the problem within a prescribed interval",
			     "Claudia", "perturb the starting point using the perturbation radius suffix information",
			     "");
  roptions->setOptionExtraInfo("random_point_type",8);

  roptions->AddLowerBoundedNumberOption("random_point_perturbation_interval",
					"Amount by which starting point is perturbed when choosing to pick random point by perturbating starting point",
					0.,true, 1.,
					"");
  roptions->setOptionExtraInfo("random_point_perturbation_interval",8);
					

  roptions->AddLowerBoundedIntegerOption
    ("num_iterations_suspect",
     "Number of iterations over which a node is considered \"suspect\" (for debugging purposes only, see detailed documentation).",
     -1,-1,
     "When the number of iterations to solve a node is above this number, the subproblem at this"
     " node is considered to be suspect and it will be outputed in a file (set to -1 to deactivate this).");
  roptions->setOptionExtraInfo("num_iterations_suspect",15);



  roptions->AddLowerBoundedIntegerOption("num_retry_unsolved_random_point",
					 "Number $k$ of times that the algorithm will try to resolve an unsolved NLP with a random starting point "
					 "(we call unsolved an NLP for which Ipopt is not "
					 "able to guarantee optimality within the specified tolerances).",
					 0,0,
					 "When Ipopt fails to solve a continuous NLP sub-problem, if $k > 0$, the algorithm will "
					 "try again to solve the failed NLP with $k$ new randomly chosen starting points "
					 " or until the problem is solved with success.");
  roptions->setOptionExtraInfo("num_retry_unsolved_random_point",15);


  roptions->SetRegisteringCategory("Options for non-convex problems", RegisteredOptions::BonminCategory);


  roptions->AddLowerBoundedIntegerOption("num_resolve_at_root",
					 "Number $k$ of tries to resolve the root node with different starting points.",
					 0,0,
					 "The algorithm will solve the root node with $k$ random starting points"
					 " and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_root",8);

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_node",
					 "Number $k$ of tries to resolve a node (other than the root) of the tree with different starting point.",
					 0,0,
					 "The algorithm will solve all the nodes with $k$ different random starting points "
					 "and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_node",8);

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_infeasibles",
					 "Number $k$ of tries to resolve an infeasible node (other than the root) of the tree with different starting point.",
					 0,0,
					 "The algorithm will solve all the infeasible nodes with $k$ different random starting points "
					 "and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_infeasibles",8);



}

static void register_OA_options
(SmartPtr<RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("Outer Approximation cuts generation", RegisteredOptions::BonminCategory);

  roptions->AddStringOption2("disjunctive_cut_type",
			     "Determine if and what kind of disjunctive cuts should be computed.",
			     "none",
			     "none", "No disjunctive cuts.",
			     "most-fractional", "If discrete variables present, compute disjunction for most-fractional variable");
  roptions->setOptionExtraInfo("disjunctive_cut_type",7);

  roptions->AddStringOption2("oa_cuts_scope","Specify if OA cuts added are to be set globally or locally valid",
                             "global",
			     "local","Cuts are treated as globally valid",
			     "global", "Cuts are treated as locally valid",
			     "");
  roptions->setOptionExtraInfo("oa_cuts_scope",7);

  roptions->AddStringOption2("add_only_violated_oa","Do we add all OA cuts or only the ones violated by current point?",
			     "no",
			     "no","Add all cuts",
			     "yes","Add only violated Cuts","");
  roptions->setOptionExtraInfo("add_only_violated_oa",7);

  roptions->AddStringOption4("cut_strengthening_type",
                             "Determines if and what kind of cut strengthening should be performed.",
                             "none",
                             "none", "No strengthening of cuts.",
                             "sglobal", "Strengthen global cuts.",
                             "uglobal-slocal", "Unstrengthened global and strengthened local cuts",
                             "sglobal-slocal", "Strengthened global and strengthened local cuts",
                             "");
  roptions->setOptionExtraInfo("cut_strengthening_type",7);

  roptions->AddLowerBoundedNumberOption("tiny_element","Value for tiny element in OA cut",
					-0.,0,1e-08,
					"We will remove \"cleanly\" (by relaxing cut) an element lower"
					" than this.");
  roptions->setOptionExtraInfo("tiny_element",7);

  roptions->AddLowerBoundedNumberOption("very_tiny_element","Value for very tiny element in OA cut",
					-0.,0,1e-17,
					"Algorithm will take the risk of neglecting an element lower"
					" than this.");
  roptions->setOptionExtraInfo("very_tiny_element",7);

  roptions->AddLowerBoundedIntegerOption("oa_cuts_log_level",
                                         "level of log when generating OA cuts.",
                                         0, 0,
                                         "0: outputs nothing,\n"
                                         "1: when a cut is generated, its violation and index of row from which it originates,\n"
                                         "2: always output violation of the cut.\n"
                                         "3: output generated cuts incidence vectors.");
  roptions->setOptionExtraInfo("oa_cuts_log_level",7);

}


///Register options
void
CouenneMINLPInterface::registerOptions
(SmartPtr<RegisteredOptions> roptions)
{
  // We try to the options - if those have been registered
  // already, we catch the exception and don't need to do it again
  try {
    register_general_options(roptions);
    register_OA_options(roptions);
#ifdef COUENNE_HAS_FILTERSQP
    FilterSolver::RegisterOptions(roptions);
#endif
#ifdef COUENNE_HAS_IPOPT
    IpoptSolver::RegisterOptions(roptions);
#endif
  }
  catch(RegisteredOptions::OPTION_ALREADY_REGISTERED) {
    // skipping
  }


}

/** To set some application specific defaults. */
void
CouenneMINLPInterface::setAppDefaultOptions(SmartPtr<OptionsList> Options)
{}

CouenneMINLPInterface::Messages::Messages
():CoinMessages((int)OSITMINLPINTERFACE_DUMMY_END)
{
  strcpy(source_ ,"NLP");
  ADD_MSG(SOLUTION_FOUND, std_m, 2,
          "After %d tries found a solution of %g (previous best %g).");

  ADD_MSG(INFEASIBLE_SOLUTION_FOUND, std_m, 2,
          "After %d tries found an solution of %g infeasible problem.");

  ADD_MSG(UNSOLVED_PROBLEM_FOUND, std_m, 2,
          "After %d tries found an solution of %g unsolved problem.");
  ADD_MSG(WARN_SUCCESS_WS, warn_m, 2,
	  "Problem not solved with warm start but solved without");

  ADD_MSG(WARNING_RESOLVING, warn_m,2,
	  "Trying to resolve NLP with different starting point (%d attempts).");
  ADD_MSG(WARN_SUCCESS_RANDOM, warn_m, 1,
	  "Problem initially not solved but solved with a random starting point (success on %d attempt)");
  ADD_MSG(WARN_CONTINUING_ON_FAILURE, warn_m, 1,
	  "Warning : continuing branching, while there are unrecovered failures at nodes");

  ADD_MSG(SUSPECT_PROBLEM, std_m, 2, "NLP number %d is suspect (see bounds and start file)");
  ADD_MSG(IPOPT_SUMMARY, std_m, 2, "Ipopt return (for %s): status %2d, iter count %4d, time %g");
  ADD_MSG(BETTER_SOL, std_m, 2, "Solution of value %g found on %d'th attempt");

  ADD_MSG(LOG_HEAD, std_m, 1,
          "\n          "
          "    Num      Status      Obj             It       time");
  ADD_MSG(LOG_FIRST_LINE, std_m, 1,
          "    %-8d %-11s %-23.16g %-8d %-3g");
  ADD_MSG(LOG_LINE, std_m, 1, " %c  r%-7d %-11s %-23.16g %-8d %-3g");

  ADD_MSG(ALTERNATE_OBJECTIVE, std_m, 1, "Objective value recomputed with alternate objective: %g.");

  ADD_MSG(WARN_RESOLVE_BEFORE_INITIAL_SOLVE, warn_m, 1,
	  "resolve called before any call to initialSol  can not use warm starts.");
  ADD_MSG(ERROR_NO_TNLPSOLVER, warn_m, 1,"Can not parse options when no IpApplication has been created");
  ADD_MSG(WARNING_NON_CONVEX_OA, warn_m, 1,
          "OA on non-convex constraint is very experimental.");
  ADD_MSG(SOLVER_DISAGREE_STATUS, warn_m, 1, "%s says problem %s, %s says %s.");
  ADD_MSG(SOLVER_DISAGREE_VALUE, warn_m, 1, "%s gives objective %.16g, %s gives %.16g.");

}


void
CouenneMINLPInterface::OaMessageHandler::print(OsiRowCut &row){
  FILE * fp = filePointer();
  const int &n = row.row().getNumElements();
  fprintf(fp,"Row cut has %d elements. Lower bound: %g, upper bound %g.\n", n, row.lb(), row.ub());
  const int * idx = row.row().getIndices();
  const double * val = row.row().getElements();
  for(int i = 0 ; i < n ; i++){
    fprintf(fp,"%g, x%d",val[i], idx[i]);
    if(i && i % 7 == 0)
      fprintf(fp,"\n");
  }
}

CouenneMINLPInterface::OaMessages::OaMessages(): CoinMessages((int) OA_MESSAGES_DUMMY_END){
  strcpy(source_,"OaCg");
  ADD_MSG(VIOLATED_OA_CUT_GENERATED, std_m, 1,"Row %d, cut violation is %g: Outer approximation cut generated.");

  ADD_MSG(CUT_NOT_VIOLATED_ENOUGH, std_m, 2,"Row %d, cut violation is %g: Outer approximation cut not generated.");

  ADD_MSG(OA_CUT_GENERATED, std_m, 1,"Row %d: Outer approximation cut not generated.");
}
bool CouenneMINLPInterface::hasPrintedOptions=0;

////////////////////////////////////////////////////////////////////
// Constructors and desctructors                                  //
////////////////////////////////////////////////////////////////////
/// Default Constructor
CouenneMINLPInterface::CouenneMINLPInterface():
  OsiSolverInterface(),
  tminlp_(NULL),
  problem_(NULL),
  problem_to_optimize_(NULL),
  feasibility_mode_(false),
  app_(NULL),
  debug_apps_(),
  testOthers_(false),
  warmstart_(NULL),
  rowsense_(NULL),
  rhs_(NULL),
  rowrange_(NULL),
  reducedCosts_(NULL),
  OsiDualObjectiveLimit_(1e200),
  hasVarNamesFile_(true),
  nCallOptimizeTNLP_(0),
  totalNlpSolveTime_(0),
  totalIterations_(0),
  maxRandomRadius_(1e08),
  randomGenerationType_(0),
  max_perturbation_(COIN_DBL_MAX),
  pushValue_(1e-02),
  numRetryInitial_(-1),
  numRetryResolve_(-1),
  numRetryInfeasibles_(-1),
  numRetryUnsolved_(1),
  messages_(),
  pretendFailIsInfeasible_(0),
  hasContinuedAfterNlpFailure_(false),
  numIterationSuspect_(-1),
  hasBeenOptimized_(false),
  obj_(NULL),
  feasibilityProblem_(NULL),
  jRow_(NULL),
  jCol_(NULL),
  jValues_(NULL),
  nnz_jac(0),
  constTypes_(NULL),
  nNonLinear_(0),
  tiny_(1e-8),
  veryTiny_(1e-20),
  infty_(1e100),
  exposeWarmStart_(false),
  firstSolve_(true),
  cutStrengthener_(NULL),
  oaMessages_(),
  oaHandler_(NULL)
{
  oaHandler_ = new OaMessageHandler;
  oaHandler_->setLogLevel(0);
}

void
CouenneMINLPInterface::initialize(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
				  Ipopt::SmartPtr<Ipopt::OptionsList> options,
				  Ipopt::SmartPtr<Ipopt::Journalist> journalist,
				  Ipopt::SmartPtr<TMINLP> tminlp){
  if(!IsValid(app_))
    createApplication(roptions, options, journalist);
  setModel(tminlp);
}

void CouenneMINLPInterface::setSolver(Ipopt::SmartPtr<TNLPSolver> app){
  app_ = app;
}


void
CouenneMINLPInterface::createApplication(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
					 Ipopt::SmartPtr<Ipopt::OptionsList> options,
					 Ipopt::SmartPtr<Ipopt::Journalist> journalist
					 )
{
  assert(!IsValid(app_));
  int ival;
  options->GetEnumValue("nlp_solver", ival, "bonmin.");
  Solver s = (Solver) ival;
  if(s == EFilterSQP){
    testOthers_ = false;;
#ifdef COUENNE_HAS_FILTERSQP
    app_ = new Bonmin::FilterSolver(roptions, options, journalist);
#else
    throw SimpleError("createApplication",
		      "Bonmin not configured to run with FilterSQP.");
#endif
  }
  else if(s == EIpopt){
    testOthers_ = false;
#ifdef COUENNE_HAS_IPOPT
    app_ = new IpoptSolver(roptions, options, journalist);
#else
    throw SimpleError("createApplication",
		      "Bonmin not configured to run with Ipopt.");
#endif
  }
  else if(s == EAll){
#ifdef COUENNE_HAS_FILTERSQP
    app_ = new Bonmin::FilterSolver(roptions, options, journalist);
#else
    throw SimpleError("createApplication",
		      "Bonmin not configured to run with Ipopt.");
#endif
#ifdef COUENNE_HAS_IPOPT
    debug_apps_.push_back(new IpoptSolver(roptions, options, journalist));
#endif
    testOthers_ = true;
  }
  if (!app_->Initialize("")) {
    throw CoinError("Error during initialization of app_","createApplication", "CouenneMINLPInterface");
  }
  for(std::list<Ipopt::SmartPtr<TNLPSolver> >::iterator i = debug_apps_.begin() ;
      i != debug_apps_.end() ; i++){
    (*i)->Initialize("");
  }
  extractInterfaceParams();

}

/** Facilitator to allocate a tminlp and an application. */
void
CouenneMINLPInterface::setModel(SmartPtr<TMINLP> tminlp)
{
  assert(IsValid(tminlp));
  tminlp_ = tminlp;
  problem_ = new TMINLP2TNLP(tminlp_);
  feasibilityProblem_ = new TNLP2FPNLP
    (SmartPtr<TNLP>(GetRawPtr(problem_)));
  if(feasibility_mode_){
    problem_to_optimize_ = GetRawPtr(feasibilityProblem_);
  }
  else {
    problem_to_optimize_ = GetRawPtr(problem_);
  }
}



void
CouenneMINLPInterface::readOptionFile(const std::string & fileName)
{
  if(IsValid(app_)){
    std::ifstream is;
    if (fileName != "") {
      try {
	is.open(fileName.c_str());
      }
      catch(std::bad_alloc) {
	throw CoinError("Not enough memory to open option file.\n", "readOptionFile", "CouenneMINLPInterface");
      }
    }
    options()->ReadFromStream(*app_->journalist(), is);
    extractInterfaceParams();
  }
}

/// Copy constructor
CouenneMINLPInterface::CouenneMINLPInterface (const CouenneMINLPInterface &source):
  OsiSolverInterface(source),
  tminlp_(source.tminlp_),
  problem_(NULL),
  problem_to_optimize_(NULL),
  feasibility_mode_(source.feasibility_mode_),
  rowsense_(NULL),
  rhs_(NULL),
  rowrange_(NULL),
  reducedCosts_(NULL),
  OsiDualObjectiveLimit_(source.OsiDualObjectiveLimit_),
  hasVarNamesFile_(source.hasVarNamesFile_),
  nCallOptimizeTNLP_(0),
  totalNlpSolveTime_(0),
  totalIterations_(0),
  maxRandomRadius_(source.maxRandomRadius_),
  randomGenerationType_(source.randomGenerationType_),
  max_perturbation_(source.max_perturbation_),
  pushValue_(source.pushValue_),
  numRetryInitial_(source.numRetryInitial_),
  numRetryResolve_(source.numRetryResolve_),
  numRetryInfeasibles_(source.numRetryInfeasibles_),
  numRetryUnsolved_(source.numRetryUnsolved_),
  messages_(),
  pretendFailIsInfeasible_(source.pretendFailIsInfeasible_),
  hasContinuedAfterNlpFailure_(source.hasContinuedAfterNlpFailure_),
  numIterationSuspect_(source.numIterationSuspect_),
  hasBeenOptimized_(source.hasBeenOptimized_),
  obj_(NULL),
  feasibilityProblem_(NULL),
  jRow_(NULL),
  jCol_(NULL),
  jValues_(NULL),
  nnz_jac(source.nnz_jac),
  constTypes_(NULL),
  //    constTypesNum_(NULL),
  nNonLinear_(0),
  tiny_(source.tiny_),
  veryTiny_(source.veryTiny_),
  infty_(source.infty_),
  exposeWarmStart_(source.exposeWarmStart_),
  firstSolve_(true),
  cutStrengthener_(source.cutStrengthener_),
  oaMessages_(),
  oaHandler_(NULL),
  strong_branching_solver_(source.strong_branching_solver_)
{
  if(defaultHandler()){
    messageHandler()->setLogLevel(source.messageHandler()->logLevel());
  }
  //Pass in message handler
  if(source.messageHandler())
    passInMessageHandler(source.messageHandler());
  // Copy options from old application
  if(IsValid(source.tminlp_)) {
    problem_ = source.problem_->clone();
    feasibilityProblem_ = new TNLP2FPNLP
      (SmartPtr<TNLP>(GetRawPtr(problem_)), source.feasibilityProblem_);
    if(feasibility_mode_)
      problem_to_optimize_ = GetRawPtr(feasibilityProblem_);
    else
      problem_to_optimize_ = GetRawPtr(problem_);
    pretendFailIsInfeasible_ = source.pretendFailIsInfeasible_;

    setAuxiliaryInfo(source.getAuxiliaryInfo());
    // Copy options from old application
    app_ = source.app_->clone();
    for(std::list<Ipopt::SmartPtr<TNLPSolver> >::const_iterator i = source.debug_apps_.begin();
        i != source.debug_apps_.end() ; i++){
      debug_apps_.push_back((*i)->clone());
    }
    testOthers_ = source.testOthers_;
  }
  else {
    throw SimpleError("Don't know how to copy an empty IpoptInterface.",
		      "copy constructor");
  }

  warmstart_ = source.warmstart_ ? source.warmstart_->clone() : NULL;

  if(source.obj_) {
    obj_ = new double[source.getNumCols()];
    CoinCopyN(source.obj_, source.getNumCols(), obj_);
  }


  oaHandler_ = new OaMessageHandler(*source.oaHandler_);;

}

OsiSolverInterface *
CouenneMINLPInterface::clone(bool copyData ) const
{
  if(copyData)
    return new CouenneMINLPInterface(*this);
  else return new CouenneMINLPInterface;
}

/// Assignment operator
CouenneMINLPInterface & CouenneMINLPInterface::operator=(const CouenneMINLPInterface& rhs)
{
  if(this!= &rhs) {
    OsiSolverInterface::operator=(rhs);
    OsiDualObjectiveLimit_ = rhs.OsiDualObjectiveLimit_;
    nCallOptimizeTNLP_ = rhs.nCallOptimizeTNLP_;
    totalNlpSolveTime_ = rhs.nCallOptimizeTNLP_;
    totalIterations_ = rhs.totalIterations_;
    maxRandomRadius_ = rhs.maxRandomRadius_;
    hasVarNamesFile_ = rhs.hasVarNamesFile_;
    pushValue_ = rhs.pushValue_;

    delete warmstart_;
    warmstart_ = NULL;

    if(IsValid(rhs.tminlp_)) {

      tminlp_ = rhs.tminlp_;
      problem_ = new TMINLP2TNLP(tminlp_);
      problem_to_optimize_ = GetRawPtr(problem_);
      app_ = rhs.app_->clone();

      warmstart_ = rhs.warmstart_ ? rhs.warmstart_->clone() : NULL;

      feasibilityProblem_ = new TNLP2FPNLP
	(SmartPtr<TNLP>(GetRawPtr(problem_)));
      nnz_jac = rhs.nnz_jac;

      if(constTypes_ != NULL) {
        delete [] constTypes_;
        constTypes_ = NULL;
      }
      if(rhs.constTypes_ != NULL) {
        constTypes_ = new TNLP::LinearityType[getNumRows()];
        CoinCopyN(rhs.constTypes_, getNumRows(), constTypes_);
      }
      /*
	if(constTypesNum_ != NULL) {
        delete [] constTypesNum_;
        constTypesNum_ = NULL;
	}
	if(rhs.constTypesNum_ != NULL) {
        constTypesNum_ = new int[getNumRows()];
        CoinCopyN(rhs.constTypesNum_, getNumRows(), constTypesNum_);
	}
      */
      if(rhs.jValues_!=NULL && rhs.jRow_ != NULL && rhs.jCol_ != NULL && nnz_jac>0) {
        jValues_ = new double [nnz_jac];
        jCol_    = new Index [nnz_jac];
        jRow_    = new Index [nnz_jac];
        CoinCopyN(rhs.jValues_ , nnz_jac,jValues_ );
        CoinCopyN(rhs.jCol_    , nnz_jac,jCol_    );
        CoinCopyN(rhs.jRow_    , nnz_jac,jRow_    );
      }
      else if(nnz_jac > 0) {
        throw CoinError("Arrays for storing jacobian are inconsistant.",
			"copy constructor",
			"IpoptOAInterface");
      }
      tiny_ = rhs.tiny_;
      veryTiny_ = rhs.veryTiny_;
      infty_ = rhs.infty_;
      exposeWarmStart_ = rhs.exposeWarmStart_;

    }
    else {
      tminlp_ =NULL;
      problem_ = NULL;
      app_ = NULL;
      feasibilityProblem_ = NULL;
    }


    if(obj_) {
      delete [] obj_;
      obj_ = NULL;
    }
    if(rhs.obj_) {
      obj_ = new double[rhs.getNumCols()];
      CoinCopyN(rhs.obj_, rhs.getNumCols(), obj_);
    }

    hasVarNamesFile_ = rhs.hasVarNamesFile_;

    nCallOptimizeTNLP_ = rhs.nCallOptimizeTNLP_;
    totalNlpSolveTime_ = rhs.totalNlpSolveTime_;
    totalIterations_ = rhs.totalIterations_;
    maxRandomRadius_ = rhs.maxRandomRadius_;
    pushValue_ = rhs.pushValue_;
    numRetryInitial_ = rhs.numRetryInitial_;
    numRetryResolve_ = rhs.numRetryResolve_;
    numRetryInfeasibles_ = rhs.numRetryInfeasibles_;
    numRetryUnsolved_ = rhs.numRetryUnsolved_;
    pretendFailIsInfeasible_ = rhs.pretendFailIsInfeasible_;
    numIterationSuspect_ = rhs.numIterationSuspect_;

    hasBeenOptimized_ = rhs.hasBeenOptimized_;
    cutStrengthener_ = rhs.cutStrengthener_;

    delete oaHandler_;
    oaHandler_ = new OaMessageHandler(*rhs.oaHandler_);
    strong_branching_solver_ = rhs.strong_branching_solver_;

    freeCachedData();
  }
  return *this;
}

const SmartPtr<OptionsList> CouenneMINLPInterface::options() const
{
  if(!IsValid(app_)) {
    messageHandler()->message(ERROR_NO_TNLPSOLVER, messages_)<<CoinMessageEol;
    return NULL;
  }
  else
    return app_->options();
}

SmartPtr<OptionsList> CouenneMINLPInterface::options()
{
  if(!IsValid(app_)) {
    messageHandler()->message(ERROR_NO_TNLPSOLVER, messages_)<<CoinMessageEol;
    return NULL;
  }
  else
    return app_->options();
}

/// Destructor
CouenneMINLPInterface::~CouenneMINLPInterface ()
{
  freeCachedData();
  delete [] jRow_;
  delete [] jCol_;
  delete [] jValues_;
  delete [] constTypes_;
  delete [] obj_;
  delete oaHandler_;
  delete warmstart_;
}

void
CouenneMINLPInterface::freeCachedColRim()
{
  if(reducedCosts_!=NULL) {
    delete []  reducedCosts_;
    reducedCosts_ = NULL;
  }

}

void
CouenneMINLPInterface::freeCachedRowRim()
{
  if(rowsense_!=NULL) {
    delete []  rowsense_;
    rowsense_ = NULL;
  }
  if(rhs_!=NULL) {
    delete []  rhs_;
    rhs_ = NULL;
  }
  if(rowrange_!=NULL) {
    delete []  rowrange_;
    rowrange_ = NULL;
  }
  //   if(dualsol_!=NULL)
  //     {
  //       delete []  dualsol_; dualsol_ = NULL;
  //     }
}

void
CouenneMINLPInterface::freeCachedData()
{
  freeCachedColRim();
  freeCachedRowRim();
}

const char * CouenneMINLPInterface::OPT_SYMB="OPT";
const char * CouenneMINLPInterface::FAILED_SYMB="FAILED";
const char * CouenneMINLPInterface::UNBOUND_SYMB="UNBOUNDED";
const char * CouenneMINLPInterface::INFEAS_SYMB="INFEAS";
///////////////////////////////////////////////////////////////////
// WarmStart Information                                                                           //
///////////////////////////////////////////////////////////////////


void
CouenneMINLPInterface::resolveForCost(int numsolve, bool keepWarmStart)
{
  // This method assumes that a problem has just been solved and we try for a
  // different solution. So disregard (in fact, clear out) any warmstart info
  // we might have, and acquire a new one before returning.
  delete warmstart_;
  warmstart_ = NULL;

  Coin::SmartPtr<SimpleReferencedPtr<CoinWarmStart> > ws_backup = NULL;
  if(!exposeWarmStart_ && keepWarmStart){
    //if warm start is not exposed, we need to store the current starting point to
    //restore it at the end of the method
    ws_backup = make_referenced(app_->getUsedWarmStart(problem_));
  }
  /** Save current optimum. */
  vector<double> point(getNumCols()*3+ getNumRows());
  double bestBound = getObjValue();
  CoinCopyN(getColSolution(),
	    getNumCols(), point());
  CoinCopyN(getRowPrice(),
	    2*getNumCols()+ getNumRows(),
	    point() + getNumCols());

  if(isProvenOptimal())
    messageHandler()->message(SOLUTION_FOUND,
			      messages_)
      <<1<<getObjValue()<<bestBound
      <<CoinMessageEol;
  else
    messageHandler()->message(INFEASIBLE_SOLUTION_FOUND,
			      messages_)
      <<1
      <<CoinMessageEol;
  for(int f = 0; f < numsolve ; f++) {
    messageHandler()->message(WARNING_RESOLVING,
			      messages_)
      <<f+1<< CoinMessageEol ;
    randomStartingPoint();
    solveAndCheckErrors(0,0,"resolveForCost");


    char c=' ';
    //Is solution better than previous
    if(isProvenOptimal() &&
       getObjValue()<bestBound) {
      c='*';
      messageHandler()->message(BETTER_SOL, messages_)<<getObjValue()<<f+1<< CoinMessageEol;
      CoinCopyN(getColSolution(),
		getNumCols(), point());
      CoinCopyN(getRowPrice(),
		2*getNumCols()+ getNumRows(),
		point() + getNumCols());
      bestBound = getObjValue();
    }

    messageHandler()->message(LOG_LINE, messages_)
      <<c<<f+1<<statusAsString()<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    if(isProvenOptimal())
      messageHandler()->message(SOLUTION_FOUND,
				messages_)
	<<f+2<<getObjValue()<<bestBound
	<<CoinMessageEol;
    else if(!isAbandoned())
      messageHandler()->message(UNSOLVED_PROBLEM_FOUND,
				messages_)
	<<f+2
	<<CoinMessageEol;
    else
      messageHandler()->message(INFEASIBLE_SOLUTION_FOUND,
				messages_)
	<<f+2
	<<CoinMessageEol;
  }
  setColSolution(point());
  setRowPrice(point() + getNumCols());
  app_->enableWarmStart();

  optimizationStatus_ = app_->ReOptimizeTNLP(GetRawPtr(problem_to_optimize_));
  hasBeenOptimized_ = true;

  if(!exposeWarmStart_ && keepWarmStart) {
    app_->setWarmStart(ws_backup->ptr(), problem_);
  }
}

void
CouenneMINLPInterface::resolveForRobustness(int numsolve)
{
  // This method assumes that a problem has just been solved and we try for a
  // different solution. So disregard (in fact, clear out) any warmstart info
  // we might have, and acquire a new one before returning.
  delete warmstart_;
  warmstart_ = NULL;


  CoinWarmStart * ws_backup = NULL;
  if(!exposeWarmStart_){
    //if warm start is not exposed, we need to store the current starting point to
    //restore it at the end of the method
    ws_backup = app_->getUsedWarmStart(problem_);
  }
  //std::cerr<<"Resolving the problem for robustness"<<std::endl;
  //First remove warm start point and resolve
  app_->disableWarmStart();
  problem()->resetStartingPoint();
  messageHandler()->message(WARNING_RESOLVING,
			    messages_)
    <<1<< CoinMessageEol ;
  solveAndCheckErrors(0,0,"resolveForRobustness");


  char c='*';
  if(isAbandoned()) {
    c=' ';
  }
  messageHandler()->message(LOG_LINE, messages_)
    <<c<<1<<statusAsString()<<getObjValue()<<app_->IterationCount()<<
    app_->CPUTime()<<CoinMessageEol;


  if(!isAbandoned()) {
    messageHandler()->message(WARN_SUCCESS_WS, messages_) << CoinMessageEol ;
    // re-enable warmstart and get it
    app_->enableWarmStart();
    if (!exposeWarmStart_) {
      app_->setWarmStart(ws_backup, problem_);
      delete ws_backup;
    }
    return; //we won go on
  }

  //still unsolved try again with different random starting points
  for(int f = 0; f < numsolve ; f++) {
    messageHandler()->message(WARNING_RESOLVING,
			      messages_)
      <<f+2<< CoinMessageEol ;

    randomStartingPoint();
    solveAndCheckErrors(0,0,"resolveForRobustness");


    messageHandler()->message(IPOPT_SUMMARY, messages_)
      <<"resolveForRobustness"<<optimizationStatus_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    char c='*';
    if(isAbandoned()) {
      c=' ';
    }
    messageHandler()->message(LOG_LINE, messages_)
      <<c<<f+2<<statusAsString()<<getObjValue()
      <<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    if(!isAbandoned()) {
      messageHandler()->message(WARN_SUCCESS_RANDOM, messages_)
	<< f+2 << CoinMessageEol ;
      // re-enable warmstart and get it
      app_->enableWarmStart();
      if (!exposeWarmStart_) {
        app_->setWarmStart(ws_backup, problem_);
        delete ws_backup;
      }

      return; //we have found a solution and continue
    }
  }


  if(!exposeWarmStart_){
    app_->setWarmStart(ws_backup, problem_);
    delete ws_backup;
  }
  if(pretendFailIsInfeasible_) {
    if(pretendFailIsInfeasible_ == 1) {
      messageHandler()->message(WARN_CONTINUING_ON_FAILURE,
				messages_)
	<<CoinMessageEol;
      hasContinuedAfterNlpFailure_ = 1;
    }
    return;
  }
  else {
    std::string probName;
    getStrParam(OsiProbName,probName);
    throw newUnsolvedError(app_->errorCode(), problem_,
                           probName);
  }
  // Do NOT get warmstart in other cases
}

////////////////////////////////////////////////////////////////////
// Problem information methods                                    //
////////////////////////////////////////////////////////////////////
/// Get number of columns
int CouenneMINLPInterface::getNumCols() const
{

  return problem_->num_variables();
}


/// Get number of rows
int
CouenneMINLPInterface::getNumRows() const
{
  return problem_->num_constraints();
}

const double *
CouenneMINLPInterface::getColLower() const
{
  return problem_->x_l();
}

const double *
CouenneMINLPInterface::getColUpper() const
{
  return problem_->x_u();
}

#if 1


///get name of variables
const OsiSolverInterface::OsiNameVec&
CouenneMINLPInterface::getVarNames() {
  return getColNames();
}
#endif


void CouenneMINLPInterface::extractSenseRhsAndRange() const
{
  assert(rowsense_==NULL&&rhs_==NULL&&rowrange_==NULL);
  int numrows = problem_->num_constraints();
  if(numrows == 0) return;
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  rowsense_ = new char [numrows];
  rhs_ = new double [numrows];
  rowrange_ = new double [numrows];
  for(int i = 0 ; i < numrows ; i++) {
    rowrange_[i]=0.;
    convertBoundToSense(rowLower[i], rowUpper[i], rowsense_[i], rhs_[i], rowrange_[i]);
  }
}

/** Get pointer to array[getNumRows()] of row constraint senses.
    <ul>
    <li>'L': <= constraint
    <li>'E': =  constraint
    <li>'G': >= constraint
    <li>'R': ranged constraint
    <li>'N': free constraint
    </ul>
*/
const char *
CouenneMINLPInterface::getRowSense() const
{
  if(rowsense_==NULL) {
    extractSenseRhsAndRange();
  }
  return rowsense_;
}

/** Get pointer to array[getNumRows()] of rows right-hand sides
    <ul>
    <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
    <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
    <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
    <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
    </ul>
*/
const double *
CouenneMINLPInterface::getRightHandSide() const
{
  if(rhs_==NULL) {
    extractSenseRhsAndRange();
  }
  return rhs_;
}

/** Get pointer to array[getNumRows()] of row ranges.
    <ul>
    <li> if rowsense()[i] == 'R' then
    rowrange()[i] == rowupper()[i] - rowlower()[i]
    <li> if rowsense()[i] != 'R' then
    rowrange()[i] is 0.0
    </ul>
*/
const double *
CouenneMINLPInterface::getRowRange() const
{
  if(rowrange_==NULL) {
    extractSenseRhsAndRange();
  }
  return rowrange_;
}

const double *
CouenneMINLPInterface::getRowLower() const
{
  return problem_->g_l();
}

const double *
CouenneMINLPInterface::getRowUpper() const
{
  return problem_->g_u();
}

/// Return true if column is continuous
bool
CouenneMINLPInterface::isContinuous(int colNumber) const
{
  return (problem_->var_types()[colNumber]==TMINLP::CONTINUOUS);
}

/// Return true if column is binary
bool
CouenneMINLPInterface::isBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==TMINLP::BINARY);
}

/** Return true if column is integer.
    Note: This function returns true if the the column
    is binary or a general integer.
*/
bool
CouenneMINLPInterface::isInteger(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==TMINLP::BINARY)||
	  (problem_->var_types()[colNumber]==TMINLP::INTEGER));
}

/// Return true if column is general integer
bool
CouenneMINLPInterface::isIntegerNonBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==TMINLP::INTEGER);
}
/// Return true if column is binary and not fixed at either bound
bool
CouenneMINLPInterface::isFreeBinary(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==TMINLP::BINARY)
	  &&((getColUpper()[colNumber]-getColLower()[colNumber]) > 1 - 1e-09));
}

/// Get solver's value for infinity
double
CouenneMINLPInterface::getInfinity() const
{
  return COIN_DBL_MAX;
}

/// Get pointer to array[getNumCols()] of primal solution vector
const double *
CouenneMINLPInterface::getColSolution() const
{
  if(hasBeenOptimized_)
    return problem_->x_sol();
  else
    return problem_->x_init();
}

/// Get pointer to array[getNumRows()] of dual prices
const double *
CouenneMINLPInterface::getRowPrice() const
{
  if(hasBeenOptimized_)
    return problem_->duals_sol();
  else
    return problem_->duals_init();
}

/// Get a pointer to array[getNumCols()] of reduced costs
const double *
CouenneMINLPInterface::getReducedCost() const
{
  (*handler_)<<"WARNING : trying to access reduced cost in Ipopt always retrun 0"<<CoinMessageEol;
  if(reducedCosts_==NULL) {
    reducedCosts_ = new double [getNumCols()];
    CoinFillN(reducedCosts_,getNumCols(),0.);
  }
  return reducedCosts_;
}

/** Get pointer to array[getNumRows()] of row activity levels (constraint
    matrix times the solution vector */
const double *
CouenneMINLPInterface::getRowActivity() const
{
  return problem_->g_sol();
}

/** Get how many iterations it took to solve the problem (whatever
    "iteration" mean to the solver.
*/
int
CouenneMINLPInterface::getIterationCount() const
{
  return app_->IterationCount();
}


/** Set a single column lower bound.
    Use -getInfinity() for -infinity. */
void
CouenneMINLPInterface::setColLower( int elementIndex, double elementValue )
{
  //  if(fabs(problem_->x_l()[elementIndex]-elementValue)>1e-06)
  problem_->SetVariableLowerBound(elementIndex,elementValue);
  hasBeenOptimized_ = false;
}

/** Set a single column upper bound.
    Use getInfinity() for infinity. */
void
CouenneMINLPInterface::setColUpper( int elementIndex, double elementValue )
{
  //  if(fabs(problem_->x_u()[elementIndex]-elementValue)>1e-06)
  problem_->SetVariableUpperBound(elementIndex,elementValue);
  hasBeenOptimized_ = false;
}

/** Set the lower bounds for all columns
    Use -getInfinity() for -infinity. */
void
CouenneMINLPInterface::setColLower( const double* array )
{
  problem_->SetVariablesLowerBounds(problem_->num_variables(),
				    array);
  hasBeenOptimized_ = false;
}

/** Set Set the upper bounds for all columns
    Use getInfinity() for infinity. */
void
CouenneMINLPInterface::setColUpper( const double* array )
{
  problem_->SetVariablesUpperBounds(problem_->num_variables(),
				    array);
  hasBeenOptimized_ = false;
}

/** Set a single row lower bound.
    Use -getInfinity() for -infinity. */
void
CouenneMINLPInterface::setRowLower( int elementIndex, double elementValue )
{
  throw SimpleError("Not implemented yet but should be if necessary.",
		    "setRowLower");
  hasBeenOptimized_ = false;
}

/** Set a single row upper bound.
    Use getInfinity() for infinity. */
void
CouenneMINLPInterface::setRowUpper( int elementIndex, double elementValue )
{
  throw SimpleError("Not implemented yet but should be if necessary.",
		    "setRowUpper");
  hasBeenOptimized_ = false;
}

/** Set the type of a single row */
void
CouenneMINLPInterface::setRowType(int index, char sense, double rightHandSide,
				  double range)
{
  throw SimpleError("Not implemented yet but should be if necessary.",
		    "setRowType");
  hasBeenOptimized_ = false;
}


/// Set the objective function sense.
/// (1 for min (default), -1 for max)
void
CouenneMINLPInterface::setObjSense(double s)
{
  throw SimpleError("Can not change objective sense of an Ipopt problem.",
		    "setObjSense");
  hasBeenOptimized_ = false;
}

/** Set the primal solution variable values

    colsol[getNumCols()] is an array of values for the primal variables.
    These values are copied to memory owned by the solver interface object
    or the solver.  They will be returned as the result of getColSolution()
    until changed by another call to setColSolution() or by a call to any
    solver routine.  Whether the solver makes use of the solution in any
    way is solver-dependent.
*/
void
CouenneMINLPInterface::setColSolution(const double *colsol)
{
  problem_->setxInit(getNumCols(), colsol);
  hasBeenOptimized_ = false;
}

/** Set dual solution variable values

    rowprice[getNumRows()] is an array of values for the dual
    variables. These values are copied to memory owned by the solver
    interface object or the solver.  They will be returned as the result of
    getRowPrice() until changed by another call to setRowPrice() or by a
    call to any solver routine.  Whether the solver makes use of the
    solution in any way is solver-dependent.
*/

void
CouenneMINLPInterface::setRowPrice(const double * rowprice)
{
  problem_->setDualsInit(getNumCols()*2 + getNumRows(), rowprice);
  hasBeenOptimized_ = false;
}

/*! \brief Get an empty warm start object

  This routine returns an empty CoinWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
*/
CoinWarmStart *
CouenneMINLPInterface::getEmptyWarmStart () const
{return app_->getEmptyWarmStart();}

/** Get warmstarting information */
CoinWarmStart*
CouenneMINLPInterface::getWarmStart() const
{
  if (exposeWarmStart_) {
    return internal_getWarmStart();;
  }
  else {
    return getEmptyWarmStart();
  }
}
/** Set warmstarting information. Return true/false depending on whether
    the warmstart information was accepted or not. */
bool
CouenneMINLPInterface::setWarmStart(const CoinWarmStart* ws)
{
  if (exposeWarmStart_) {
    return internal_setWarmStart(ws);
  }
  else {
    return true;
  }
}
/** Get warmstarting information */
CoinWarmStart*
CouenneMINLPInterface::internal_getWarmStart() const
{
  if (exposeWarmStart_ && warmstart_) {
    return warmstart_->clone();
  }
  else {
    return getEmptyWarmStart();
  }
}
/** Set warmstarting information. Return true/false depending on whether
    the warmstart information was accepted or not. */
bool
CouenneMINLPInterface::internal_setWarmStart(const CoinWarmStart* ws)
{
  delete warmstart_;
  warmstart_ = NULL;
  hasBeenOptimized_ = false;
  if (exposeWarmStart_) {
    if (ws == NULL) {
      return true;
    }
    if(app_->warmStartIsValid(ws)) {
      warmstart_ = ws->clone();
      return true;
    }
    // See if it is anything else than the CoinWarmStartBasis that all others
    // derive from
    const CoinWarmStartPrimalDual* pdws =
      dynamic_cast<const CoinWarmStartPrimalDual*>(ws);
    if (pdws) {
      // Must be an IpoptWarmStart, since the others do not derive from this.
      // Accept it.
      warmstart_ = new IpoptWarmStart(*pdws);
      return true;
    }
    return false;
  }
  else {
    return true;
  }
}

/** Set the index-th variable to be a continuous variable */
void
CouenneMINLPInterface::setContinuous(int index)
{
  problem_->SetVariableType(index, TMINLP::CONTINUOUS);
  hasBeenOptimized_ = false;
}
/** Set the index-th variable to be an integer variable */
void
CouenneMINLPInterface::setInteger(int index)
{
  problem_->SetVariableType(index, TMINLP::INTEGER);
  hasBeenOptimized_ = false;
}

/// Get objective function value (can't use default)
double
CouenneMINLPInterface::getObjValue() const
{
  return problem_->obj_value();
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
CouenneMINLPInterface::setIntParam(OsiIntParam key, int value)
{
  //  debugMessage("OsiCpxSolverInterface::setIntParam(%d, %d)\n", key, value);

  bool retval = false;
  switch (key) {
  case OsiMaxNumIteration:
    retval = false;
    break;
  case OsiMaxNumIterationHotStart:
    if( value >= 0 ) {
      retval = false;
    }
    else
      retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  default:
    retval = false;
    (*handler_)<< "Unhandled case in setIntParam\n"<<CoinMessageEol;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
CouenneMINLPInterface::setDblParam(OsiDblParam key, double value)
{
  //  debugMessage("CouenneMINLPInterface::setDblParam(%d, %g)\n", key, value);

  bool retval = false;
  switch (key) {
  case OsiDualObjectiveLimit:
    OsiDualObjectiveLimit_ = value;
    retval = true;
    break;
  case OsiPrimalObjectiveLimit:
    (*handler_)<<"Can not set primal objective limit parameter"<<CoinMessageEol;
    retval = false;
    break;
  case OsiDualTolerance:
    (*handler_)<<"Can not set dual tolerance parameter"<<CoinMessageEol;
    retval = false;
    break;
  case OsiPrimalTolerance:
    (*handler_)<<"Can not set primal tolerance parameter"<<CoinMessageEol;
    retval = false;
  case OsiObjOffset:
    retval = OsiSolverInterface::setDblParam(key,value);
    break;
  case OsiLastDblParam:
    retval = false;
    break;
  default:
    retval = false;
    (*handler_) << "Unhandled case in setDblParam"<<CoinMessageEol;
    break;
  }
  return retval;
}


//-----------------------------------------------------------------------------

bool
CouenneMINLPInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  //  debugMessage("CouenneMINLPInterface::setStrParam(%d, %s)\n", key, value.c_str());

  bool retval=false;
  switch (key) {
  case OsiProbName:
    OsiSolverInterface::setStrParam(key,value);
    return retval = true;
  case OsiSolverName:
    return false;
  case OsiLastStrParam:
    return false;
  }
  return false;
}

//-----------------------------------------------------------------------------

bool
CouenneMINLPInterface::getIntParam(OsiIntParam key, int& value) const
{
  //  debugMessage("CouenneMINLPInterface::getIntParam(%d)\n", key);

  value = -COIN_INT_MAX; // Give a dummy value
  bool retval = false;
  switch (key) {
  case OsiMaxNumIteration:
    retval = false;
    break;
  case OsiMaxNumIterationHotStart:
    retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  default:
    retval = false;
    (*handler_) << "Unhandled case in setIntParam"<<CoinMessageEol;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
CouenneMINLPInterface::getDblParam(OsiDblParam key, double& value) const
{
  //  debugMessage("CouenneMINLPInterface::getDblParam(%d)\n", key);

  bool retval = false;
  switch (key) {
  case OsiDualObjectiveLimit:
    value = OsiDualObjectiveLimit_;
    retval = true;
    break;
  case OsiPrimalObjectiveLimit:
    value = getInfinity();
    retval = true;
    break;
  case OsiDualTolerance:
    retval = false;
    break;
  case OsiPrimalTolerance:
    options()->GetNumericValue("tol", value,"");
    //value = 1e-07;
    retval = true;
    break;
  case OsiObjOffset:
    retval = OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiLastDblParam:
    retval = false;
    break;
  }
  return retval;
}


//-----------------------------------------------------------------------------

bool
CouenneMINLPInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  //  debugMessage("CouenneMINLPInterface::getStrParam(%d)\n", key);

  switch (key) {
  case OsiProbName:
    OsiSolverInterface::getStrParam(key, value);
    break;
  case OsiSolverName:
    value = "Ipopt";
    break;
  case OsiLastStrParam:
    return false;
  }

  return true;
}

void
CouenneMINLPInterface::randomStartingPoint()
{
  int numcols = getNumCols();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  double * sol = new double[numcols];
  const Number * x_init = problem_->x_init_user();
  const double* perturb_radius = NULL;
  if (randomGenerationType_ == perturb_suffix) {
    const TMINLP::PerturbInfo* pertubinfo = tminlp_->perturbInfo();
    if (pertubinfo) {
      perturb_radius = pertubinfo->GetPerturbationArray();
    }
    if (!perturb_radius) {
      throw SimpleError("Can't use perturb_radius if no radii are given.",
			"randomStartingPoint");
    }
  }
  for(int i = 0 ; i < numcols ; i++) {
    int randomGenerationType = randomGenerationType_;
    if(x_init[i] < colLower[i] || x_init[i] > colUpper[i])
      randomGenerationType = uniform;
    if(randomGenerationType_ == uniform){
      double lower = std::min(-maxRandomRadius_,colUpper[i] - maxRandomRadius_);
      lower = std::max(colLower[i], lower);
      double upper = std::max(maxRandomRadius_,colLower[i] + maxRandomRadius_);
      upper = std::min(colUpper[i],upper);
      lower = std::min(upper,lower);
      upper = std::max(upper, lower);
      double interval = upper - lower;
      sol[i] = CoinDrand48()*(interval) + lower;}
    else if (randomGenerationType_ == perturb){
      const double lower = std::max(x_init[i] - max_perturbation_, colLower[i]);
      const double upper = std::min(x_init[i] + max_perturbation_, colUpper[i]);
      const double interval = upper - lower;
      sol[i]  = lower + CoinDrand48()*(interval);
    }
    else if (randomGenerationType_ == perturb_suffix){
      const double radius = perturb_radius[i];
      const double lower = std::max(x_init[i] - radius*max_perturbation_, colLower[i]);
      const double upper = std::min(x_init[i] + radius*max_perturbation_, colUpper[i]);
      const double interval = upper - lower;
      sol[i]  = lower + CoinDrand48()*(interval);
    }
  }
  app_->disableWarmStart();
  setColSolution(sol);
  delete [] sol;
}



/** This methods initialiaze arrays for storing the jacobian */
int CouenneMINLPInterface::initializeJacobianArrays()
{
  Index n, m, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  tminlp_->get_nlp_info( n, m, nnz_jac, nnz_h_lag, index_style);

  if(jRow_ != NULL) delete jRow_;
  if(jCol_ != NULL) delete jCol_;
  if(jValues_ != NULL) delete jValues_;

  jRow_ = new Index[nnz_jac];
  jCol_ = new Index[nnz_jac];
  jValues_ = new Number[nnz_jac];
  tminlp_->eval_jac_g(n, NULL, 0, m, nnz_jac, jRow_, jCol_, NULL);
  if(index_style == Ipopt::TNLP::FORTRAN_STYLE)//put C-style
    {
      for(int i = 0 ; i < nnz_jac ; i++){
	jRow_[i]--;
	jCol_[i]--;
      }
    }

  if(constTypes_ != NULL) delete [] constTypes_;
  //  if(constTypesNum_ != NULL) delete [] constTypesNum_;

  constTypes_ = new TNLP::LinearityType[getNumRows()];
  tminlp_->get_constraints_linearity(getNumRows(), constTypes_);
  //  constTypesNum_ = new int[getNumRows()];
  for(int i = 0; i < getNumRows() ; i++) {
    if(constTypes_[i]==TNLP::NON_LINEAR) {
      //constTypesNum_[i] =
      nNonLinear_++;
    }
  }
  return nnz_jac;
}


double
CouenneMINLPInterface::getConstraintsViolation(const double *x, double &obj)
{
  int numcols = getNumCols();
  int numrows = getNumRows();
  double * g = new double[numrows];
  tminlp_->eval_g(numcols, x, 1, numrows, g);
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();


  double norm = 0;
  for(int i = 0; i< numrows ; i++) {
    if(!constTypes_ || constTypes_[i] == TNLP::NON_LINEAR) {
      double rowViolation = 0;
      if(rowLower[i] > -1e10)
	rowViolation = std::max(0.,rowLower[i] - g[i]);

      if(rowUpper[i] < 1e10);
      rowViolation = std::max(rowViolation, g[i] - rowUpper[i]);

      norm = rowViolation > norm ? rowViolation : norm;
    }
  }
  tminlp_->eval_f(numcols, x, 1, obj);
  delete [] g;
  return norm;
}

/** get infinity norm of constraint violation of a point and objective error*/
double
CouenneMINLPInterface::getNonLinearitiesViolation(const double *x, const double obj)
{
  double f;
  double norm = getConstraintsViolation(x, f);
  assert((f - obj) > -1e-08);
  norm =  (f - obj) > norm ? f - obj : norm;
  return norm;
}



//A procedure to try to remove small coefficients in OA cuts (or make it non small
static inline
bool cleanNnz(double &value, double colLower, double colUpper,
	      double rowLower, double rowUpper, double colsol,
	      double & lb, double &ub, double tiny, double veryTiny)
{
  if(fabs(value)>= tiny) return 1;

  if(fabs(value)<veryTiny) return 0;//Take the risk?

  //try and remove
  double infty = 1e20;
  bool colUpBounded = colUpper < 10000;
  bool colLoBounded = colLower > -10000;
  bool rowNotLoBounded =  rowLower <= - infty;
  bool rowNotUpBounded = rowUpper >= infty;
  bool pos =  value > 0;

  if(colLoBounded && pos && rowNotUpBounded) {
    lb += value * (colsol - colLower);
    return 0;
  }
  else
    if(colLoBounded && !pos && rowNotLoBounded) {
      ub += value * (colsol - colLower);
      return 0;
    }
    else
      if(colUpBounded && !pos && rowNotUpBounded) {
        lb += value * (colsol - colUpper);
        return 0;
      }
      else
        if(colUpBounded && pos && rowNotLoBounded) {
          ub += value * (colsol - colUpper);
          return 0;
        }
  //can not remove coefficient increase it to smallest non zero
  if(pos) value = tiny;
  else
    value = - tiny;
  return 1;
}

/** Get the outer approximation constraints at the point x.
 */
void
CouenneMINLPInterface::getOuterApproximation(OsiCuts &cs, const double * x, bool getObj, const double * x2, double theta, bool global)
{

  //printf ("this procedure has been called\n");

  int n,m, nnz_jac_g, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  problem_to_optimize_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
  if(jRow_ == NULL || jCol_ == NULL || jValues_ == NULL)
    initializeJacobianArrays();
  assert(jRow_ != NULL);
  assert(jCol_ != NULL);
  double * g = new double[m];

  problem_to_optimize_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, jValues_);
  problem_to_optimize_->eval_g(n,x,1,m,g);

  //double *xx = new double [n];

  // x ---> xx
  //
  //
  //




  problem_to_optimize_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, jValues_);
  problem_to_optimize_->eval_g(n,x,1,m,g);
  //As jacobian is stored by cols fill OsiCuts with cuts
  CoinPackedVector * cuts = new CoinPackedVector[nNonLinear_ + 1];
  double * lb = new double[nNonLinear_ + 1];
  double * ub = new double[nNonLinear_ + 1];

  int * row2cutIdx = new int[m];//store correspondance between index of row and index of cut (some cuts are not generated for rows because linear, or not binding). -1 if constraint does not generate a cut, otherwise index in cuts.
  int numCuts = 0;

  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice() + 2 * n;

  for (int i=0; i< n; i++)
    printf ("variable %d: [%g,%g]\n", i, rowLower [i], rowUpper [i]);

  double infty = getInfinity();
  double nlp_infty = infty_;

  for(int rowIdx = 0; rowIdx < m ; rowIdx++) {
    if(constTypes_[rowIdx] == TNLP::NON_LINEAR) {
#if 0
      if(fabs(duals[rowIdx]) == 0.)
	{
	  row2cutIdx[rowIdx] = -1;
#ifdef NDEBUG
	  (*handler_)<<"non binding constraint"<<CoinMessageEol;
#endif
	  continue;
	}
#endif
      row2cutIdx[rowIdx] = numCuts;
      if(rowLower[rowIdx] > - nlp_infty)
        lb[numCuts] = rowLower[rowIdx] - g[rowIdx];
      else
        lb[numCuts] = - infty;
      if(rowUpper[rowIdx] < nlp_infty)
        ub[numCuts] = rowUpper[rowIdx] - g[rowIdx];
      else
        ub[numCuts] = infty;
      if(rowLower[rowIdx] > -infty && rowUpper[rowIdx] < infty)
	{
	  if(duals[rowIdx] >= 0)// <= inequality
	    lb[numCuts] = - infty;
	  if(duals[rowIdx] <= 0)// >= inequality
	    ub[numCuts] = infty;
	}

      numCuts++;
    }
    else
      row2cutIdx[rowIdx] = -1;
  }


  for(int i = 0 ; i < nnz_jac_g ; i++) {
    const int &rowIdx = jRow_[i];
    const int & cutIdx = row2cutIdx[ rowIdx ];
    if(cutIdx != -1) {
      const int &colIdx = jCol_[i];
      //"clean" coefficient
      if(cleanNnz(jValues_[i],colLower[colIdx], colUpper[colIdx],
		  rowLower[rowIdx], rowUpper[rowIdx],
		  x[colIdx],
		  lb[cutIdx],
		  ub[cutIdx], tiny_, veryTiny_)) {
        cuts[cutIdx].insert(colIdx,jValues_[i]);
        if(lb[cutIdx] > - infty)
          lb[cutIdx] += jValues_[i] * x[colIdx];
        if(ub[cutIdx] < infty)
	  ub[cutIdx] += jValues_[i] * x[colIdx];
      }
    }
  }

  int * cut2rowIdx = NULL;
  if (IsValid(cutStrengthener_) || oaHandler_->logLevel() > 0) {
    cut2rowIdx = new int [numCuts];// Store correspondance between indices of cut and indices of rows. For each cut
    for(int rowIdx = 0 ; rowIdx < m ; rowIdx++){
      if(row2cutIdx[rowIdx] >= 0){
	cut2rowIdx[row2cutIdx[rowIdx]] = rowIdx;
      }
    }
  }

  for(int cutIdx = 0; cutIdx < numCuts; cutIdx++) {
    //Compute cut violation
    if(x2 != NULL) {
      double rhs = cuts[cutIdx].dotProduct(x2);
      double violation = 0.;
      violation = std::max(violation, rhs - ub[cutIdx]);
      violation = std::max(violation, lb[cutIdx] - rhs);
      if(violation < theta) {
        if(oaHandler_->logLevel() > 0)
          oaHandler_->message(CUT_NOT_VIOLATED_ENOUGH, oaMessages_)<<cut2rowIdx[cutIdx]<<violation<<CoinMessageEol;
        continue;}
      if(oaHandler_->logLevel() > 0)
	oaHandler_->message(VIOLATED_OA_CUT_GENERATED, oaMessages_)<<cut2rowIdx[cutIdx]<<violation<<CoinMessageEol;
    }
    else if (oaHandler_->logLevel() > 0)
      oaHandler_->message(OA_CUT_GENERATED, oaMessages_)<<cut2rowIdx[cutIdx]<<CoinMessageEol;
    OsiRowCut newCut;
    //    if(lb[i]>-1e20) assert (ub[i]>1e20);

    if (IsValid(cutStrengthener_)) {
      const int& rowIdx = cut2rowIdx[cutIdx];
      bool retval =
	cutStrengthener_->ComputeCuts(cs, GetRawPtr(tminlp_),
				      GetRawPtr(problem_), rowIdx,
				      cuts[cutIdx], lb[cutIdx], ub[cutIdx], g[rowIdx],
				      rowLower[rowIdx], rowUpper[rowIdx],
				      n, x, infty);
      if (!retval) {
	(*messageHandler()) << "error in cutStrengthener_->ComputeCuts\n";
	//exit(-2);
      }
    }
    if(global) {
      newCut.setGloballyValidAsInteger(1);
    }
    newCut.setEffectiveness(99.99e99);
    newCut.setLb(lb[cutIdx]);
    newCut.setUb(ub[cutIdx]);
    newCut.setRow(cuts[cutIdx]);
    //    CftValidator validator;
    //    validator(newCut);
    if(oaHandler_->logLevel()>2){
      oaHandler_->print(newCut);}
    cs.insert(newCut);
  }

  delete[] g;
  delete [] cuts;
  delete [] row2cutIdx;
  delete [] cut2rowIdx;

  if(getObj && ! problem_->hasLinearObjective()) { // Get the objective cuts
    double * obj = new double [n];
    problem_to_optimize_->eval_grad_f(n, x, 1,obj);
    double f;
    problem_to_optimize_->eval_f(n, x, 1, f);

    CoinPackedVector v;
    v.reserve(n);
    lb[nNonLinear_] = -f;
    ub[nNonLinear_] = -f;
    //double minCoeff = 1e50;
    for(int i = 0; i<n ; i++) {
      if(cleanNnz(obj[i],colLower[i], colUpper[i],
		  -getInfinity(), 0,
		  x[i],
		  lb[nNonLinear_],
		  ub[nNonLinear_],tiny_, 1e-15)) {
        //	      minCoeff = std::min(fabs(obj[i]), minCoeff);
        v.insert(i,obj[i]);
        lb[nNonLinear_] += obj[i] * x[i];
        ub[nNonLinear_] += obj[i] * x[i];
      }
    }
    v.insert(n,-1);
    //Compute cut violation
    bool genCut = true;
    if(x2 != NULL) {
      double rhs = v.dotProduct(x2);
      double violation = std::max(0., rhs - ub[nNonLinear_]);
      if(violation < theta) genCut = false;
    }
    if(genCut) {
      if (IsValid(cutStrengthener_)) {
	lb[nNonLinear_] = -infty;
	bool retval =
	  cutStrengthener_->ComputeCuts(cs, GetRawPtr(tminlp_),
					GetRawPtr(problem_), -1,
					v, lb[nNonLinear_], ub[nNonLinear_],
					ub[nNonLinear_], -infty, 0.,
					n, x, infty);
	if (!retval) {
	  (*handler_)<< "error in cutStrengthener_->ComputeCuts"<<CoinMessageEol;
	  //exit(-2);
	}
      }
      OsiRowCut newCut;
      if(global)
	newCut.setGloballyValidAsInteger(1);
      newCut.setEffectiveness(99.99e99);
      newCut.setRow(v);
      newCut.setLb(-COIN_DBL_MAX/*Infinity*/);
      newCut.setUb(ub[nNonLinear_]);
      //     CftValidator validator;
      //     validator(newCut);
      cs.insert(newCut);
    }
    delete [] obj;
  }

  delete []lb;
  delete[]ub;
}



/** Get the outer approximation of a single constraint at the point x.
 */
void
CouenneMINLPInterface::getConstraintOuterApproximation(OsiCuts &cs, int rowIdx,
						       const double * x,
						       const double * x2, bool global)
{
  double g;
  int * indices = new int[getNumCols()];
  double * values = new double[getNumCols()];
  int nnz;
  problem_->eval_grad_gi(getNumCols(), x, 1, rowIdx, nnz, indices, values);
  problem_->eval_gi(getNumCols(),x,1, rowIdx, g);

  CoinPackedVector cut;
  double lb;
  double ub;


  const double rowLower = getRowLower()[rowIdx];
  const double rowUpper = getRowUpper()[rowIdx];
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double dual = (getRowPrice() + 2 * getNumCols())[rowIdx];
  double infty = getInfinity();
  double nlp_infty = infty_;

  if(rowLower > - nlp_infty)
    lb = rowLower - g;
  else
    lb = - infty;
  if(rowUpper < nlp_infty)
    ub = rowUpper - g;
  else
    ub = infty;
  if(rowLower > -infty && rowUpper < infty)
    {
      if(dual >= 0)// <= inequality
	lb = - infty;
      if(dual <= 0)// >= inequality
	ub = infty;
    }

  for(int i = 0 ; i < nnz; i++) {
    const int &colIdx = indices[i];
    //"clean" coefficient
    if(cleanNnz(values[i],colLower[colIdx], colUpper[colIdx],
		rowLower, rowUpper,
		x[colIdx],
		lb,
		ub, tiny_, veryTiny_)) {
      cut.insert(colIdx,values[i]);
      if(lb > - infty)
	lb += values[i] * x[colIdx];
      if(ub < infty)
	ub += values[i] * x[colIdx];
    }
  }

  OsiRowCut newCut;

  if(global) {
    newCut.setGloballyValidAsInteger(1);
  }
  newCut.setEffectiveness(99.99e99);
  newCut.setLb(lb);
  newCut.setUb(ub);
  newCut.setRow(cut);
  cs.insert(newCut);

  delete [] indices;
  delete [] values;
}

void
CouenneMINLPInterface::switchToFeasibilityProblem(int n,const double * x_bar,const int *inds,
						  double a, double s, int L){
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_use_feasibility_pump_objective(true);
  feasibilityProblem_->set_dist2point_obj(n,(const Number *) x_bar,(const Index *) inds);
  feasibilityProblem_->setLambda(a);
  feasibilityProblem_->setSigma(s);
  feasibilityProblem_->setNorm(L);
  feasibilityProblem_->set_use_cutoff_constraint(false);
  feasibilityProblem_->set_use_local_branching_constraint(false);
  problem_to_optimize_ = GetRawPtr(feasibilityProblem_);
  feasibility_mode_ = true;
}

void
CouenneMINLPInterface::switchToFeasibilityProblem(int n,const double * x_bar,const int *inds,
						  double rhs_local_branching_constraint){
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_use_feasibility_pump_objective(false);
  feasibilityProblem_->set_dist2point_obj(n,(const Number *) x_bar,(const Index *) inds);
  feasibilityProblem_->set_use_cutoff_constraint(false);
  feasibilityProblem_->set_use_local_branching_constraint(true);
  feasibilityProblem_->set_rhs_local_branching_constraint(rhs_local_branching_constraint);
  problem_to_optimize_ = GetRawPtr(feasibilityProblem_);
  feasibility_mode_ = true;
}

void
CouenneMINLPInterface::switchToOriginalProblem(){
  problem_to_optimize_ = GetRawPtr(problem_);
  feasibility_mode_ = false;
}

double
CouenneMINLPInterface::solveFeasibilityProblem(int n,const double * x_bar,const int *inds,
					       double a, double s, int L)
{
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_use_feasibility_pump_objective(true);
  feasibilityProblem_->set_dist2point_obj(n,(const Number *) x_bar,(const Index *) inds);
  feasibilityProblem_->setLambda(a);
  feasibilityProblem_->setSigma(s);
  feasibilityProblem_->setNorm(L);
  feasibilityProblem_->set_use_cutoff_constraint(false);
  feasibilityProblem_->set_use_local_branching_constraint(false);
  nCallOptimizeTNLP_++;
  totalNlpSolveTime_-=CoinCpuTime();
  SmartPtr<TNLPSolver> app2 = app_->clone();
  app2->options()->SetIntegerValue("print_level", (Index) 0);
  optimizationStatus_ = app2->OptimizeTNLP(GetRawPtr(feasibilityProblem_));
  totalNlpSolveTime_+=CoinCpuTime();
  hasBeenOptimized_=true;
  return getObjValue();
}

double
CouenneMINLPInterface::solveFeasibilityProblem(int n,const double * x_bar,const int *inds,
					       int L, double cutoff)
{
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_use_feasibility_pump_objective(true);
  feasibilityProblem_->set_dist2point_obj(n,(const Number *) x_bar,(const Index *) inds);
  feasibilityProblem_->setLambda(1.0);
  feasibilityProblem_->setSigma(0.0);
  feasibilityProblem_->setNorm(L);
  feasibilityProblem_->set_use_cutoff_constraint(true);
  feasibilityProblem_->set_cutoff(cutoff);
  feasibilityProblem_->set_use_local_branching_constraint(false);
  nCallOptimizeTNLP_++;
  totalNlpSolveTime_-=CoinCpuTime();
  SmartPtr<TNLPSolver> app2 = app_->clone();
  app2->options()->SetIntegerValue("print_level", (Index) 0);
  optimizationStatus_ = app2->OptimizeTNLP(GetRawPtr(feasibilityProblem_));
  totalNlpSolveTime_+=CoinCpuTime();
  hasBeenOptimized_=true;
  return getObjValue();
}

double
CouenneMINLPInterface::getFeasibilityOuterApproximation(int n,const double * x_bar,const int *inds, OsiCuts &cs, bool addOnlyViolated, bool global)
{
  double ret_val = solveFeasibilityProblem(n, x_bar, inds, 1, 0, 2);
  getOuterApproximation(cs, getColSolution(), 0, (addOnlyViolated? x_bar:NULL)
			, global);
  return ret_val;
}


static bool WarnedForNonConvexOa=false;

void
CouenneMINLPInterface::extractLinearRelaxation(OsiSolverInterface &si,
					       const double * x, bool getObj)
{
  double * rowLow = NULL;
  double * rowUp = NULL;

  int n;
  int m;
  int nnz_jac_g;
  int nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  //Get problem information
  problem_to_optimize_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);

  //if not allocated allocate spaced for stroring jacobian
  if(jRow_ == NULL || jCol_ == NULL || jValues_ == NULL)
    initializeJacobianArrays();

  //get Jacobian
  problem_to_optimize_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, jValues_);


  double *g = new double[m];
  problem_to_optimize_->eval_g(n, x, 1, m, g);

  rowLow = new double[m];
  rowUp = new double[m];
  int * nonBindings = new int[m];//store non binding constraints (which are to be removed from OA)
  int numNonBindings = 0;
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice() + 2*n;
  assert(m==getNumRows());
  double infty = si.getInfinity();
  double nlp_infty = infty_;

  for(int i = 0 ; i < m ; i++) {
    if(constTypes_[i] == TNLP::NON_LINEAR) {
      //If constraint is range not binding prepare to remove it
      if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty && fabs(duals[i]) == 0.)
	{
	  nonBindings[numNonBindings++] = i;
	  continue;
	}
      else
        if(rowLower[i] > - nlp_infty){
          rowLow[i] = (rowLower[i] - g[i]) - 1e-07;
          if(! WarnedForNonConvexOa && rowUpper[i] < nlp_infty){
	    messageHandler()->message(WARNING_NON_CONVEX_OA, messages_)<<CoinMessageEol;
	    WarnedForNonConvexOa = true;
          }
        }
	else
	  rowLow[i] = - infty;
      if(rowUpper[i] < nlp_infty)
        rowUp[i] =  (rowUpper[i] - g[i]) + 1e-07;
      else
        rowUp[i] = infty;

      //If equality or ranged constraint only add one side by looking at sign of dual multiplier
      if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty)
	{
	  if(duals[i] >= 0.)// <= inequality
	    rowLow[i] = - infty;
	  if(duals[i] <= 0.)// >= inequality
	    rowUp[i] = infty;
	}
    }
    else {
      if(rowLower[i] > -nlp_infty){
	//   printf("Lower %g ", rowLower[i]);
	rowLow[i] = (rowLower[i] - g[i]);
      }
      else
        rowLow[i] = - infty;
      if(rowUpper[i] < nlp_infty){
	//   printf("Upper %g ", rowUpper[i]);
	rowUp[i] =  (rowUpper[i] - g[i]);
      }
      else
        rowUp[i] = infty;
    }
  }



  //Then convert everything to a CoinPackedMatrix
  //Go through values, clean coefficients and fix bounds
  for(int i = 0 ; i < nnz_jac_g ; i++) {
    if(constTypes_[jRow_[i]] != TNLP::LINEAR){//For linear just copy is fine.
      if(//For other clean tinys
	 cleanNnz(jValues_[i],colLower[jCol_[i]], colUpper[jCol_[i]],
		  rowLower[jRow_[i]], rowUpper[jRow_[i]],
		  x[jCol_[i]],
		  rowLow[jRow_[i]],
		  rowUp[jRow_[i]], tiny_, veryTiny_)) {
	rowLow[jRow_[i]] += jValues_[i] * x[jCol_ [i]];
	rowUp[jRow_[i]] += jValues_[i] *x[jCol_[i]];
      }
    }
    else {
      double value = jValues_[i] * getColSolution()[jCol_[i]];
      rowLow[jRow_[i]] += value;
      rowUp[jRow_[i]] += value;
    }
  }
  CoinPackedMatrix mat(true, jRow_, jCol_, jValues_, nnz_jac_g);
  mat.setDimensions(m,n); // In case matrix was empty, this should be enough

  //remove non-bindings equality constraints
  mat.deleteRows(numNonBindings, nonBindings);

  int numcols=getNumCols();
  double *obj = new double[numcols];
  for(int i = 0 ; i < numcols ; i++)
    obj[i] = 0.;


  si.loadProblem(mat, getColLower(), getColUpper(), obj, rowLow, rowUp);
  delete [] rowLow;
  delete [] rowUp;
  delete [] nonBindings;
  delete [] g;
  for(int i = 0 ; i < getNumCols() ; i++) {
    if(isInteger(i))
      si.setInteger(i);
  }
  if(getObj) {
    bool addObjVar = false;
    if(problem_->hasLinearObjective()){
      //Might be in trouble if objective has a constant part
      // for now just check that f(0) = 0.
      // If it is not adding a constant term does not seem supported by Osi
      // for now
      double zero;
      problem_to_optimize_->eval_f(n, obj, 1, zero);
      si.setDblParam(OsiObjOffset, -zero);
      //if(fabs(zero - 0) > 1e-10)
      //addObjVar = true;
      //else {
      //Copy the linear objective and don't create a dummy variable.
      problem_to_optimize_->eval_grad_f(n, x, 1,obj);
      si.setObjective(obj);
      //}
    }
    else {
      addObjVar = true;
    }

    if(addObjVar){
      //add variable alpha
      //(alpha should be empty in the matrix with a coefficient of -1 and unbounded)
      CoinPackedVector a;
      si.addCol(a,-si.getInfinity(), si.getInfinity(), 1.);

      // Now get the objective cuts
      // get the gradient, pack it and add the cut
      problem_to_optimize_->eval_grad_f(n, x, 1,obj);
      double ub;
      problem_to_optimize_->eval_f(n, x, 1, ub);
      ub*=-1;
      double lb = -1e300;
      CoinPackedVector objCut;
      CoinPackedVector * v = &objCut;
      v->reserve(n);
      for(int i = 0; i<n ; i++) {
	if(nnz_jac_g)
	  {
	    if(cleanNnz(obj[i],colLower[i], colUpper[i],
			-getInfinity(), 0,
			x[i],
			lb,
			ub, tiny_, veryTiny_)) {
	      v->insert(i,obj[i]);
	      lb += obj[i] * x[i];
	      ub += obj[i] * x[i];
	    }
	  }
	else //Unconstrained problem can not put clean coefficient
	  {
	    if(cleanNnz(obj[i],colLower[i], colUpper[i],
			-getInfinity(), 0,
			x[i],
			lb,
			ub, 1e-03, 1e-08)) {
	      v->insert(i,obj[i]);
	      lb += obj[i] * x[i];
	      ub += obj[i] * x[i];
	    }
	  }
      }
      v->insert(n,-1);
      si.addRow(objCut, lb, ub);
    }
  }
  //  si.writeMpsNative("OA.mps",NULL, NULL, 1);
  delete [] obj;
}

/** Add a collection of linear cuts to problem formulation.*/
void
CouenneMINLPInterface::applyRowCuts(int numberCuts, const OsiRowCut * cuts)
{
  if(numberCuts)
    freeCachedRowRim();
  const OsiRowCut ** cutsPtrs = new const OsiRowCut*[numberCuts];
  for(int i = 0 ; i < numberCuts ; i++)
    {
      cutsPtrs[i] = &cuts[i];
    }
  problem_->addCuts(numberCuts, cutsPtrs);
  delete [] cutsPtrs;
}

void
CouenneMINLPInterface::solveAndCheckErrors(bool warmStarted, bool throwOnFailure,
					   const char * whereFrom)
{
  totalNlpSolveTime_-=CoinCpuTime();
  if(warmStarted)
    optimizationStatus_ = app_->ReOptimizeTNLP(GetRawPtr(problem_to_optimize_));
  else
    optimizationStatus_ = app_->OptimizeTNLP(GetRawPtr(problem_to_optimize_));
  totalNlpSolveTime_+=CoinCpuTime();
  nCallOptimizeTNLP_++;
  hasBeenOptimized_ = true;


  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    //app_->Options()->SetIntegerValue("print_level",0, true, true);
    app_->options()->SetStringValue("print_user_options","no", false, true);
  }

  bool otherDisagree = false ;
#if 0
  if(optimizationStatus_ == TNLPSolver::notEnoughFreedom)//Too few degrees of freedom
    {
      (*messageHandler())<<"Too few degrees of freedom...."<<CoinMessageEol;
      int numrows = getNumRows();
      int numcols = getNumCols();

      const double * colLower = getColLower();
      const double * colUpper = getColUpper();


      const double * rowLower = getRowLower();
      const double * rowUpper = getRowUpper();

      int numberFixed = 0;
      for(int i = 0 ; i < numcols ; i++)
	{
	  if(colUpper[i] - colLower[i] <= INT_BIAS)
	    {
	      numberFixed++;
	    }
	}
      int numberEqualities = 0;
      for(int i = 0 ; i < numrows ; i++)
	{
	  if(rowUpper[i] - rowLower[i] <= INT_BIAS)
	    {
	      numberEqualities++;
	    }	
	}
      if(numcols - numberFixed > numberEqualities || numcols < numberEqualities)
	{
	  std::string probName;
	  getStrParam(OsiProbName, probName);
	  throw newUnsolvedError(app_->errorCode(), problem_, probName);
	}
      double * saveColLow = CoinCopyOfArray(getColLower(), getNumCols());
      double * saveColUp = CoinCopyOfArray(getColUpper(), getNumCols());

      for(int i = 0; i < numcols && numcols - numberFixed <= numberEqualities ; i++)
	{
	  if(colUpper[i] - colLower[i] <= INT_BIAS)
	    {
	      setColLower(i, saveColLow[i]-1e-06);
	      setColUpper(i, saveColLow[i]+1e-06);
	      numberFixed--;
	    }
	}
      solveAndCheckErrors(warmStarted, throwOnFailure, whereFrom);
      //restore
      for(int i = 0; i < numcols && numcols - numberFixed < numrows ; i++)
	{
	  problem_->SetVariableLowerBound(i,saveColLow[i]);
	  problem_->SetVariableUpperBound(i,saveColUp[i]);
	}
      delete [] saveColLow;
      delete [] saveColUp;
      return;
    }
  else
#endif
    if(!app_->isRecoverable(optimizationStatus_))//Solver failed and the error can not be recovered, throw it
      {
	std::string probName;
	getStrParam(OsiProbName, probName);
	throw newUnsolvedError(app_->errorCode(), problem_, probName);
      }
    else if(testOthers_ && !app_->isError(optimizationStatus_)){
      Ipopt::SmartPtr<TMINLP2TNLP> problem_copy = problem_->clone();
      //Try other solvers and see if they agree
      int f =1;
      for(std::list<Ipopt::SmartPtr<TNLPSolver> >::iterator i = debug_apps_.begin();
          i != debug_apps_.end() ; i++){
        TNLPSolver::ReturnStatus otherStatus = (*i)->OptimizeTNLP(GetRawPtr(problem_copy));
	messageHandler()->message(LOG_LINE, messages_)
	  <<'d'<<f++<<statusAsString(otherStatus)<<problem_copy->obj_value()
	  <<(*i)->IterationCount()<<(*i)->CPUTime()<<CoinMessageEol;
        if(!(*i)->isError(otherStatus)){
	  CoinRelFltEq eq(1e-05);
	  if(otherStatus != optimizationStatus_){
	    otherDisagree = true;
	    messageHandler()->message(SOLVER_DISAGREE_STATUS, messages_)
	      <<app_->solverName()<<statusAsString()
	      <<(*i)->solverName()<<statusAsString(otherStatus)<<CoinMessageEol;
	  }
	  else if(isProvenOptimal() && !eq(problem_->obj_value(),problem_copy->obj_value()))
	    {
	      otherDisagree = true;
	      messageHandler()->message(SOLVER_DISAGREE_VALUE, messages_)
		<<app_->solverName()<<problem_->obj_value()
		<<(*i)->solverName()<<problem_copy->obj_value()<<CoinMessageEol;
	    }
        }
      }
    }
  try{
    totalIterations_ += app_->IterationCount();
  }
  catch(SimpleError &E)
    {
      if (throwOnFailure)//something failed throw
	{
	  throw SimpleError("No statistics available from Ipopt",whereFrom);
	}
      else {
	return;
      }
    }
  if(problem_->hasUpperBoundingObjective()){//Check if solution is integer and recompute objective value using alternative objective function
    const double * sol = getColSolution();
    bool integerSol = true;
    double intTol = 1e-08;
    if(objects()){
      int nObjects = numberObjects();
      OsiObject ** object = objects();
      for(int i = 0 ; i< nObjects ; i++){
        int dummy;
        if(object[i]->infeasibility(this,dummy) > intTol)
	  {
	    integerSol=false;
	    break;
	  }
      }
    }
    else{//Only works for integer constraints
      int numcols = getNumCols();
      for(int i = 0 ; i < numcols ; i++){
        if(isInteger(i) || isBinary(i)){
          if(fabs(sol[i] - floor(sol[i]+0.5)) > intTol){
            integerSol = false;
            break;
          }
        }
      }
    }
    if(integerSol&&isProvenOptimal()){
      double help= problem_->evaluateUpperBoundingFunction(sol);


      OsiAuxInfo * auxInfo = getAuxiliaryInfo();
      Bonmin::AuxInfo * bonInfo = dynamic_cast<Bonmin::AuxInfo *>(auxInfo);
      if(bonInfo!=0)
	{
	
	  if(help<bonInfo->bestObj2())
	    {
	      bonInfo->setBestObj2(help);
	      bonInfo->setBestSolution2(getNumCols(), const_cast<double *>(getColSolution()));

	      messageHandler()->message(ALTERNATE_OBJECTIVE, messages_)
		<<help<<CoinMessageEol;
	    }
	}
      else {
        printf("\nWARNING: the algorithm selected does not consider the second objective function\n");
      }
    }
  }
  messageHandler()->message(IPOPT_SUMMARY, messages_)
    <<whereFrom<<optimizationStatus_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;

  if((nCallOptimizeTNLP_ % 20) == 1)
    messageHandler()->message(LOG_HEAD, messages_)<<CoinMessageEol;


  if ( (numIterationSuspect_ >= 0 && (getIterationCount()>numIterationSuspect_ || isAbandoned())) ||
       ( otherDisagree )){
    messageHandler()->message(SUSPECT_PROBLEM,
                              messages_)<<nCallOptimizeTNLP_<<CoinMessageEol;
    std::string subProbName;
    getStrParam(OsiProbName, subProbName);
    std::ostringstream os;
    os<<"_"<<nCallOptimizeTNLP_;
    subProbName+=os.str();
    problem_->outputDiffs(subProbName, NULL/*getVarNames()*/);
  }

}

////////////////////////////////////////////////////////////////////
// Solve Methods                                                  //
////////////////////////////////////////////////////////////////////
/// Solve initial continuous relaxation
void CouenneMINLPInterface::initialSolve()
{
  assert(IsValid(app_));
  assert(IsValid(problem_));

  // Discard warmstart_ if we had one
  delete warmstart_;
  warmstart_ = NULL;

  if(!hasPrintedOptions) {
    int printOptions;
    app_->options()->GetEnumValue("print_user_options",printOptions,"bonmin.");
    if(printOptions)
      app_->options()->SetStringValue("print_user_options","yes",true,true);
  }
  if(exposeWarmStart_)
    app_->disableWarmStart();
  solveAndCheckErrors(0,1,"initialSolve");

  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    app_->options()->SetStringValue("print_user_options","no");
    app_->options()->SetIntegerValue("print_level",0);
  }

  messageHandler()->message(LOG_FIRST_LINE, messages_)<<nCallOptimizeTNLP_
						      <<statusAsString()
                                                      <<getObjValue()
                                                      <<app_->IterationCount()
                                                      <<app_->CPUTime()
                                                      <<CoinMessageEol;

  int numRetry = firstSolve_ ? numRetryInitial_ : numRetryResolve_;
  if(isAbandoned()) {
    resolveForRobustness(numRetryUnsolved_);
  }
  else if(numRetry)
    {
      resolveForCost(numRetry, numRetryInitial_ > 0);
      /** Only do it once at the root.*/
      numRetryInitial_ = 0;
    }
  firstSolve_ = false;

  // if warmstart_ is not empty then had to use resolveFor... and that created
  // the warmstart at the end, and we have nothing to do here. Otherwise...
  if (! warmstart_ && ! isAbandoned()) {
    if (exposeWarmStart_) {
      warmstart_ = app_->getWarmStart(problem_);
    }
  }
}

/** Resolve the continuous relaxation after problem modification.
 * \note for Ipopt, same as resolve */
void
CouenneMINLPInterface::resolve()
{
  assert(IsValid(app_));
  assert(IsValid(problem_));

  int has_warmstart = warmstart_ == NULL ? 0 : 1;
  if(warmstart_ == NULL) has_warmstart = 0;
  else if(!app_->warmStartIsValid(warmstart_)) has_warmstart = 1;
  else has_warmstart = 2;
  if (has_warmstart < 2) {
    // empty (or unrecognized) warmstart
    initialSolve();
    return;
  }
  app_->setWarmStart(warmstart_, problem_);
  delete warmstart_;
  warmstart_ = NULL;

  if (INT_BIAS > 0.) {
    app_->options()->SetStringValue("warm_start_same_structure", "yes");
  }
  else {
    app_->options()->SetStringValue("warm_start_same_structure", "no");
  }

  if(problem_->duals_init() != NULL)
    app_->enableWarmStart();
  else app_->disableWarmStart();
  solveAndCheckErrors(1,1,"resolve");

  messageHandler()->message(LOG_FIRST_LINE, messages_)<<nCallOptimizeTNLP_
						      <<statusAsString()
                                                      <<getObjValue()
                                                      <<app_->IterationCount()
                                                      <<app_->CPUTime()<<CoinMessageEol;

  if(isAbandoned()) {
    resolveForRobustness(numRetryUnsolved_);
  }
  else if(numRetryResolve_ ||
	  (numRetryInfeasibles_ && isProvenPrimalInfeasible() ))
    resolveForCost(std::max(numRetryResolve_, numRetryInfeasibles_), 0);

  // if warmstart_ is not empty then had to use resolveFor... and that created
  // the warmstart at the end, and we have nothing to do here. Otherwise...
  if (! warmstart_ && ! isAbandoned()) {
    if (exposeWarmStart_) {
      warmstart_ = app_->getWarmStart(problem_);
    }
  }
}


////////////////////////////////////////////////////////////////
// Methods returning info on how the solution process terminated  //
////////////////////////////////////////////////////////////////
/// Are there a numerical difficulties?
bool CouenneMINLPInterface::isAbandoned() const
{
  return (
	  (optimizationStatus_==TNLPSolver::iterationLimit)||
	  (optimizationStatus_==TNLPSolver::computationError)||
	  (optimizationStatus_==TNLPSolver::illDefinedProblem)||
	  (optimizationStatus_==TNLPSolver::illegalOption)||
	  (optimizationStatus_==TNLPSolver::externalException)|
	  (optimizationStatus_==TNLPSolver::exception)
	  );
}

/// Is optimality proven?
bool CouenneMINLPInterface::isProvenOptimal() const
{
  return (optimizationStatus_==TNLPSolver::solvedOptimal) ||
    (optimizationStatus_==TNLPSolver::solvedOptimalTol);
}
/// Is primal infeasiblity proven?
bool CouenneMINLPInterface::isProvenPrimalInfeasible() const
{
  return (optimizationStatus_ == TNLPSolver::provenInfeasible);
}
/// Is dual infeasiblity proven?
bool CouenneMINLPInterface::isProvenDualInfeasible() const
{
  return (optimizationStatus_ == TNLPSolver::unbounded);
}
/// Is the given primal objective limit reached?
bool CouenneMINLPInterface::isPrimalObjectiveLimitReached() const
{
  (*handler_)<<"Warning : isPrimalObjectiveLimitReached not implemented yet"<<CoinMessageEol;
  return 0;
}
/// Is the given dual objective limit reached?
bool CouenneMINLPInterface::isDualObjectiveLimitReached() const
{
  //  (*messageHandler_)<<"Warning : isDualObjectiveLimitReached not implemented yet"<<CoinMessageEol;
  return (optimizationStatus_==TNLPSolver::unbounded);

}
/// Iteration limit reached?
bool CouenneMINLPInterface::isIterationLimitReached() const
{
  return (optimizationStatus_==TNLPSolver::iterationLimit);
}

void
CouenneMINLPInterface::extractInterfaceParams()
{
  if (IsValid(app_)) {
    int logLevel;
    app_->options()->GetIntegerValue("nlp_log_level", logLevel,"bonmin.");
    messageHandler()->setLogLevel(logLevel);

#ifdef COUENNE_HAS_FILTERSQP
    FilterSolver * filter = dynamic_cast<FilterSolver *>(GetRawPtr(app_));

    bool is_given =
#endif
      app_->options()->GetNumericValue("max_random_point_radius",maxRandomRadius_,"bonmin.");

#ifdef COUENNE_HAS_FILTERSQP
    if(filter && !is_given){
      // overwriting default for filter
      maxRandomRadius_ = 10.;
    }
#endif

    int oaCgLogLevel = 0;
    app_->options()->GetIntegerValue("oa_cuts_log_level", oaCgLogLevel,"bonmin.");
    oaHandler_->setLogLevel(oaCgLogLevel);

    int exposeWs = false;
    app_->options()->GetEnumValue("warm_start", exposeWs, "bonmin.");
    setExposeWarmStart(exposeWs > 0);

    app_->options()->GetIntegerValue("num_retry_unsolved_random_point", numRetryUnsolved_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_root", numRetryInitial_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_node", numRetryResolve_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_infeasibles", numRetryInfeasibles_,"bonmin.");
    app_->options()->GetIntegerValue("num_iterations_suspect", numIterationSuspect_,"bonmin.");
    app_->options()->GetEnumValue("nlp_failure_behavior",pretendFailIsInfeasible_,"bonmin.");
    app_->options()->GetNumericValue
      ("warm_start_bound_frac" ,pushValue_,"bonmin.");
    app_->options()->GetNumericValue("tiny_element",tiny_,"bonmin.");
    app_->options()->GetNumericValue("very_tiny_element",veryTiny_,"bonmin.");
    app_->options()->GetNumericValue("random_point_perturbation_interval",max_perturbation_,"bonmin.");
    app_->options()->GetEnumValue("random_point_type",randomGenerationType_,"bonmin.");
    int cut_strengthening_type;
    app_->options()->GetEnumValue("cut_strengthening_type", cut_strengthening_type,"bonmin.");

    if (cut_strengthening_type != CS_None) {
      // TNLP solver to be used in the cut strengthener
      cutStrengthener_ = new CutStrengthener(app_->clone(), app_->options());
    }
  }
}

void
CouenneMINLPInterface::SetStrongBrachingSolver(SmartPtr<StrongBranchingSolver> strong_branching_solver)
{
  strong_branching_solver_ = strong_branching_solver;
}

//#define STRONG_COMPARE
#ifdef STRONG_COMPARE
static double objorig;
#endif

void
CouenneMINLPInterface::markHotStart()
{
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::markHotStart();
    objorig = getObjValue();
#endif
    optimizationStatusBeforeHotStart_ = optimizationStatus_;
    strong_branching_solver_->markHotStart(this);
  }
  else {
    // Default Implementation
    OsiSolverInterface::markHotStart();
  }
}

void
CouenneMINLPInterface::solveFromHotStart()
{
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::solveFromHotStart();
    double obj_nlp = getObjValue() - objorig;
#endif
    optimizationStatus_ = strong_branching_solver_->solveFromHotStart(this);
    hasBeenOptimized_ = true;
#ifdef STRONG_COMPARE
    double obj_other = getObjValue() - objorig;
    printf("AWDEBUG: Strong Branching results: NLP = %15.8e Other = %15.8e\n",
	   obj_nlp, obj_other);
#endif
  }
  else {
    // Default Implementation
    OsiSolverInterface::solveFromHotStart();
  }
}

void
CouenneMINLPInterface::unmarkHotStart()
{
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::unmarkHotStart();
#endif
    strong_branching_solver_->unmarkHotStart(this);
    optimizationStatus_ = optimizationStatusBeforeHotStart_;
  }
  else {
    // Default Implementation
    OsiSolverInterface::unmarkHotStart();
  }
}

const double * CouenneMINLPInterface::getObjCoefficients() const
{
  const int n = getNumCols();
  delete [] obj_;
  obj_ = NULL;
  obj_ = new double[n];

  bool new_x = true;
  const double* x_sol = problem_->x_sol();
  bool retval = problem_->eval_grad_f(n, x_sol, new_x, obj_);

  if (!retval) {
    // Let's see if that happens - it will cause a crash
    fprintf(stderr, "ERROR WHILE EVALUATING GRAD_F in CouenneMINLPInterface::getObjCoefficients()\n");
    delete [] obj_;
    obj_ = NULL;
  }

  return obj_;
}


}/** end namespace Bonmin*/

#endif
