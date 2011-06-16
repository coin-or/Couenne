/* $Id$
 *
 * Name:    CouenneFeasPump.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Constructors and service methods of the Feasibility Pump class
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <string>

#include "CouenneConfig.h"
#include "CouenneFeasPump.hpp"
#include "CouenneFPpool.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprSum.hpp"

using namespace Couenne;

// common code for initializing ipopt application
void CouenneFeasPump::initIpoptApp () {

  // Although app_ is only used in CouenneFPSolveNLP, we need to have
  // an object lasting the program's lifetime as otherwise it appears
  // to delete the nlp pointer at deletion.

  if (!app_)
    app_ = IpoptApplicationFactory ();

  ApplicationReturnStatus status = app_ -> Initialize ();

  app_ -> Options () -> SetIntegerValue ("max_iter", 1000);
  app_ -> Options () -> SetIntegerValue // 0 for none, 4 for summary, 5 for iteration output
    ("print_level", (problem_ -> Jnlst () -> ProduceOutput (J_ERROR,         J_NLPHEURISTIC) ? 4 : 
		     problem_ -> Jnlst () -> ProduceOutput (J_STRONGWARNING, J_NLPHEURISTIC) ? 5 : 0));

  if (status != Solve_Succeeded)
    printf ("FP: Error in initialization\n");
}


// Constructor ////////////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (CouenneProblem *couenne,
				  CouenneCutGenerator *cg,
				  Ipopt::SmartPtr<Ipopt::OptionsList> options):

  CbcHeuristic         (),
  problem_             (couenne),
  couenneCG_           (cg),
  nlp_                 (NULL),
  app_                 (NULL),
  milp_                (NULL),
  postlp_              (NULL),
  pool_                (NULL),
  numberSolvePerLevel_ (-1),
  betaNLP_             (0.),
  betaMILP_            (0.),
  compDistInt_         (FP_DIST_INT),
  milpCuttingPlane_    (FP_CUT_NONE),
  nSepRounds_          (0),
  maxIter_             (COIN_INT_MAX),
  useSCIP_             (false),
  milpMethod_          (0) {

  if (IsValid (options)) {

    std::string s;

    int compareTerm;

    options -> GetIntegerValue ("feas_pump_iter",       maxIter_,             "couenne.");
    options -> GetIntegerValue ("feas_pump_level",      numberSolvePerLevel_, "couenne.");
    options -> GetIntegerValue ("feas_pump_milpmethod", milpMethod_,          "couenne.");

    options -> GetNumericValue ("feas_pump_beta_nlp",   betaNLP_,             "couenne.");
    options -> GetNumericValue ("feas_pump_beta_milp",  betaMILP_,            "couenne.");
				    
    options -> GetStringValue  ("feas_pump_convcuts", s, "couenne."); 

    milpCuttingPlane_ = 
      (s == "none")       ? FP_CUT_NONE       :
      (s == "integrated") ? FP_CUT_INTEGRATED : 
      (s == "postcut")    ? FP_CUT_POST       : FP_CUT_EXTERNAL;

    options -> GetIntegerValue ("feas_pump_nseprounds", nSepRounds_, "couenne."); 

    options -> GetStringValue  ("feas_pump_vardist",  s, "couenne."); 

    compDistInt_ = 
      (s == "integer") ? FP_DIST_INT : 
      (s == "all")     ? FP_DIST_ALL : FP_DIST_POST;

    options -> GetIntegerValue ("feas_pump_milpmethod", milpMethod_, "couenne."); 
    options -> GetIntegerValue ("feas_pump_poolcomp",   compareTerm, "couenne."); 

    pool_ = new CouenneFPpool ((enum what_to_compare) compareTerm);

    options -> GetStringValue  ("feas_pump_usescip",    s,           "couenne."); 

#ifdef COIN_HAS_SCIP
    useSCIP_ = (s == "yes");
    if (milpMethod_ < 0)
      milpMethod_ = 0;
#else
    if (s == "yes") 
      problem_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "Warning: you have set feas_pump_usescip to true, but SCIP is not installed.\n");
#endif

  } else
    pool_ = new CouenneFPpool (SUM_NINF);

  setHeuristicName ("Couenne Feasibility Pump");

  initIpoptApp ();
}

  
// Copy constructor ///////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (const CouenneFeasPump &other):

  CbcHeuristic         (other),
  problem_             (other. problem_),
  couenneCG_           (other. couenneCG_),
  nlp_                 (other. nlp_),
  app_                 (NULL),
  milp_                (other.milp_   ? other. milp_   -> clone () : NULL),
  postlp_              (other.postlp_ ? other. postlp_ -> clone () : NULL),
  pool_                (NULL),
  numberSolvePerLevel_ (other. numberSolvePerLevel_),
  betaNLP_             (other. betaNLP_),
  betaMILP_            (other. betaMILP_),
  compDistInt_         (other. compDistInt_),
  milpCuttingPlane_    (other. milpCuttingPlane_),
  nSepRounds_          (other. nSepRounds_),
  maxIter_             (other. maxIter_),
  useSCIP_             (other. useSCIP_),
  milpMethod_          (other. milpMethod_) {

  if (other. pool_)
    pool_ = new CouenneFPpool (*(other. pool_));

  for (std::set <CouenneFPsolution, compareSol>::iterator i = other.tabuPool_.begin (); 
       i != other.tabuPool_.end ();
       ++i)
    tabuPool_. insert (CouenneFPsolution (*i));

  initIpoptApp ();
}


// Clone //////////////////////////////////////////////////////// 
CbcHeuristic *CouenneFeasPump::clone () const
{return new CouenneFeasPump (*this);}


// Assignment operator ////////////////////////////////////////// 
CouenneFeasPump &CouenneFeasPump::operator= (const CouenneFeasPump & rhs) {

  if (this != &rhs) {

    CbcHeuristic::operator= (rhs);

    problem_             = rhs. problem_;
    couenneCG_           = rhs. couenneCG_;
    nlp_                 = rhs. nlp_;
    app_                 = NULL;
    milp_                = rhs. milp_   ? rhs. milp_   -> clone () : NULL;
    postlp_              = rhs. postlp_ ? rhs. postlp_ -> clone () : NULL;
    pool_                = NULL;
    numberSolvePerLevel_ = rhs. numberSolvePerLevel_;
    betaNLP_             = rhs. betaNLP_;
    betaMILP_            = rhs. betaMILP_;
    compDistInt_         = rhs. compDistInt_;
    milpCuttingPlane_    = rhs. milpCuttingPlane_;
    nSepRounds_          = rhs. nSepRounds_;
    maxIter_             = rhs. maxIter_;
    useSCIP_             = rhs. useSCIP_;
    milpMethod_          = rhs. milpMethod_;

    if (rhs. pool_)
      pool_ = new CouenneFPpool (*(rhs. pool_));

    for (std::set <CouenneFPsolution, compareSol>::iterator i = rhs.tabuPool_.begin (); 
	 i != rhs.tabuPool_.end ();
	 ++i)
      tabuPool_. insert (CouenneFPsolution (*i));

    initIpoptApp ();
  }

  return *this;
}


// Destructor /////////////////////////////////////////////////// 
CouenneFeasPump::~CouenneFeasPump () {

  if (pool_)   delete pool_;
  if (app_)    delete app_;
  if (milp_)   delete milp_;
  if (postlp_) delete postlp_;

  //if (nlp_) delete nlp_; // already deleted by "delete app_;"
}


/// Set new expression as the NLP objective function using
/// argument as point to minimize distance from. Return new
/// objective function.
expression *CouenneFeasPump::updateNLPObj (const double *iSol) {

  expression **list = new expression * [problem_ -> nVars ()];

  int nTerms = 0;

  if (betaNLP_ == 0.) {

    // here the objective function is ||x-x^0||_2^2

    // create the argument list (x_i - x_i^0)^2 for all i's
    for (int i=0; i<problem_ -> nVars (); i++) {

      if (compDistInt_ == FP_DIST_INT && !(problem_ -> Var (i) -> isInteger ()))
	continue;

      CouNumber iS = iSol [i];

      expression *base;

      if      (iS == 0.) base =              new exprClone (problem_ -> Var (i));
      else if (iS <  0.) base = new exprSum (new exprClone (problem_ -> Var (i)), new exprConst (-iS));
      else               base = new exprSub (new exprClone (problem_ -> Var (i)), new exprConst  (iS));

      list [nTerms++] = new exprPow (base, new exprConst (2.));
    }

  } else {

    // here the objective function is 
    //
    // ||P(x-x^0)||_2^2 = (x-x^0)' P'P (x-x^0)
    //
    // with P'P positive semidefinite
  }

  return new exprSum (list, nTerms);
}


/// Reads a (possibly fractional) solution and fixes the integer
/// components in the nonlinear problem for later re-solve
void CouenneFeasPump::fixIntVariables (double *sol) {

  assert (sol);

  for (int i = problem_ -> nVars (); i--;)

    if (problem_ -> Var (i) -> isInteger ()) {

      double 
	value = sol [i],
	rUp   = ceil  (value - COUENNE_EPS),
	rDn   = floor (value + COUENNE_EPS);

      // If numerics or sol[i] fractional, set to closest

      value = 
	(rUp < rDn + 0.5)           ? rUp : 
	(rUp - value < value - rDn) ? rUp : rDn;

      problem_ -> Lb (i) = 
      problem_ -> Ub (i) = value;
    }
}


/// initialize options
void CouenneFeasPump::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2
    ("feas_pump_heuristic",
     "Apply the nonconvex Feasibility Pump",
     "no",
     "no","",
     "yes","",
     "An implementation of the Feasibility Pump for nonconvex MINLPs");

  roptions -> AddLowerBoundedIntegerOption
    ("feas_pump_level",
     "Specify the logarithm of the number of feasibility pumps to perform" 
     " on average for each level of given depth of the tree.",
     -1,
     2, "Solve as many nlp's at the nodes for each level of the tree. "
     "Nodes are randomly selected. If for a "
     "given level there are less nodes than this number nlp are solved for every nodes. "
     "For example if parameter is 8, nlp's are solved for all node until level 8, " 
     "then for half the node at level 9, 1/4 at level 10.... "
     "Value -1 specify to perform at all nodes.");

  roptions -> AddLowerBoundedIntegerOption
    ("feas_pump_iter",
     "Number of iterations in the main Feasibility Pump loop",
     -1,
     10, "-1 means no limit");

  roptions -> AddBoundedNumberOption
    ("feas_pump_beta_nlp",
     "Weight of the Lagrangian Hessian in computing the objective function of the NLP problem",
     0., false,
     1., false,
     0., "0 for distance only, 1 for lagrangian hessian only");

  roptions -> AddBoundedNumberOption
    ("feas_pump_beta_milp",
     "Weight of the Lagrangian Hessian in computing the objective function of the MILP problem",
     0., false,
     1., false,
     0., "0 for distance only, 1 for lagrangian hessian only");

  // roptions -> AddStringOption2
  //   ("feas_pump_lincut",
  //    "Shortcut to linearization cutting plane applied to the MILP solution instead of solving NLPs",
  //    "no",
  //    "no","",
  //    "yes","",
  //    "");

  roptions -> AddStringOption3
    ("feas_pump_vardist",
     "Distance computed on integer-only or on both types of variables, in different flavors.",
     "integer",
     "integer",         "Only compute the distance based on integer coordinates (use post-processing if numerical errors occur)",
     "all",             "Compute the distance using continuous and integer variables",
     "int-postprocess", "Use a post-processing fixed-IP LP to determine a closest-point solution");

  roptions -> AddStringOption4
    ("feas_pump_convcuts",
     "Separate MILP-feasible, MINLP-infeasible solution during or after MILP solver.",
     "none",
     "integrated", "Done within the MILP solver in a branch-and-cut fashion",
     "external",   "Done after the MILP solver, in a Benders-like fashion",
     "postcut",    "Do one round of cuts and proceed with NLP",
     "none",       "Just proceed to the NLP");

  roptions -> AddBoundedIntegerOption
    ("feas_pump_nseprounds",
     "Number of rounds that separate convexification cuts. Must be at least 1",
     1, 1e5, 4,
     "");

  roptions -> AddStringOption3
    ("feas_pump_poolmgt",
     "Retrieval of MILP solutions when the one returned is unsatisfactory",
     "pool",
     "pool",       "Use a solution pool and replace unsatisfactory solution with Euclidean-closest in pool",
     "perturb",    "Randomly perturb unsatisfactory solution",
     "none",       "Bail out of feasibility pump");

  roptions -> AddStringOption2
    ("feas_pump_usescip",
     "Should SCIP be used to solve the MILPs?",
     "no",
     "no",  "Use Cbc's branch-and-cut to solve the MILP",
     "yes", "Use SCIP's branch-and-cut or heuristics (see feas_pump_milpmethod option) to solve the MILP",
     "");

  roptions -> AddBoundedIntegerOption
    ("feas_pump_milpmethod",
     "How should the integral solution be constructed?",
     -1, 4, -1, 
       "0: automatic, 1: aggressive heuristics, large node limit, 2: default, node limit, 3: RENS, 4: Objective Feasibility Pump,  -1: solve MILP completely");

  roptions -> AddBoundedIntegerOption
    ("feas_pump_poolcomp",
     "Priority field to compare solutions in FP pool",
     0, 2, 0, 
       "0: total number of infeasible objects (integer and nonlinear), 1: maximum infeasibility (integer or nonlinear), 2: objective value.");
}
