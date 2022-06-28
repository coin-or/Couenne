/*
 *
 * Name:    CouenneFeasPumpConstructors.cpp
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
#include "CouenneTNLP.hpp"
#include "CouenneSparseMatrix.hpp"

using namespace Ipopt;
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
    ("print_level", (problem_ -> Jnlst () -> ProduceOutput (J_ITERSUMMARY,  J_NLPHEURISTIC) ? 4 :
		     problem_ -> Jnlst () -> ProduceOutput (J_MOREDETAILED, J_NLPHEURISTIC) ? 5 : 0));

  app_ -> Options () -> SetStringValue ("fixed_variable_treatment", "make_parameter");

  // Suppress iteration output from nonlinear solver
  app_ -> Options () -> SetStringValue ("sb", "yes", false, true);

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

  numberSolvePerLevel_ (5), // if options are not valid, don't overuse FP

  multDistNLP_         (1.), // settings for classical FP
  multHessNLP_         (0.),
  multObjFNLP_         (0.),

  multDistMILP_        (1.),
  multHessMILP_        (0.),
  multObjFMILP_        (0.),

  compDistInt_         (FP_DIST_INT),
  milpCuttingPlane_    (FP_CUT_NONE),
  nSepRounds_          (0),
  maxIter_             (COIN_INT_MAX),
  useSCIP_             (false),
  milpMethod_          (0),
  tabuMgt_             (FP_TABU_NONE),
  nCalls_              (0),
  fadeMult_            (1) {

  int compareTerm = INTEGER_VARS;

  if (IsValid (options)) {

    std::string s;

    options -> GetIntegerValue ("feas_pump_iter",       maxIter_,             "couenne.");
    options -> GetIntegerValue ("feas_pump_level",      numberSolvePerLevel_, "couenne.");
    options -> GetIntegerValue ("feas_pump_milpmethod", milpMethod_,          "couenne.");

    options -> GetNumericValue ("feas_pump_mult_dist_nlp",  multDistNLP_,     "couenne.");
    options -> GetNumericValue ("feas_pump_mult_hess_nlp",  multHessNLP_,     "couenne.");
    options -> GetNumericValue ("feas_pump_mult_objf_nlp",  multObjFNLP_,     "couenne.");

    options -> GetNumericValue ("feas_pump_mult_dist_milp", multDistMILP_,    "couenne.");
    options -> GetNumericValue ("feas_pump_mult_hess_milp", multHessMILP_,    "couenne.");
    options -> GetNumericValue ("feas_pump_mult_objf_milp", multObjFMILP_,    "couenne.");

    options -> GetNumericValue ("feas_pump_fademult",  fadeMult_,     "couenne.");

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

    options -> GetStringValue  ("feas_pump_tabumgt", s, "couenne.");

    tabuMgt_ =
      (s == "pool")    ? FP_TABU_POOL    :
      (s == "perturb") ? FP_TABU_PERTURB :
      (s == "cut")     ? FP_TABU_CUT     : FP_TABU_NONE;

    options -> GetStringValue  ("feas_pump_usescip", s, "couenne.");

#ifdef COUENNE_HAS_SCIP
    useSCIP_ = (s == "yes");
    if (milpMethod_ < 0)
      milpMethod_ = 0;
#else
    if (s == "yes")
      problem_ -> Jnlst () -> Printf (J_ERROR, J_COUENNE, "Warning: you have set feas_pump_usescip to true, but SCIP is not installed.\n");
#endif

  }

  pool_ = new CouenneFPpool (problem_, (enum what_to_compare) compareTerm);

  //pool_ = new CouenneFPpool (SUM_NINF);
  //pool_ = new CouenneFPpool (problem_, INTEGER_VARS);

  setHeuristicName ("Couenne Feasibility Pump");

  initIpoptApp ();
}


// Copy constructor /////////////////////////////////////////////
CouenneFeasPump::CouenneFeasPump (const CouenneFeasPump &other)
{operator= (other);}


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

    multDistNLP_         = rhs. multDistNLP_;
    multHessNLP_         = rhs. multHessNLP_;
    multObjFNLP_         = rhs. multObjFNLP_;
			       	
    multDistMILP_        = rhs. multDistMILP_;
    multHessMILP_        = rhs. multHessMILP_;
    multObjFMILP_        = rhs. multObjFMILP_;

    compDistInt_         = rhs. compDistInt_;
    milpCuttingPlane_    = rhs. milpCuttingPlane_;
    nSepRounds_          = rhs. nSepRounds_;
    maxIter_             = rhs. maxIter_;
    useSCIP_             = rhs. useSCIP_;
    milpMethod_          = rhs. milpMethod_;
    tabuMgt_             = rhs. tabuMgt_;
    nCalls_              = rhs. nCalls_;
    fadeMult_            = rhs. fadeMult_;

    if (rhs. pool_)
      pool_ = new CouenneFPpool (*(rhs. pool_));

    for (std::set <CouenneFPsolution, compareSol>::const_iterator i = rhs.tabuPool_.begin ();
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

  expression **list = NULL;

  int nTerms = 0;

  const double *iS = iSol;

  // TODO:
  //
  // 1) resize list
  // 2) separate H norm from distance
  // 3

  if ((multHessNLP_ == 0.) ||
      (nlp_ -> optHessian () == NULL)) {

    list = new expression * [1 + problem_ -> nVars ()];

    // here the objective function is ||x-x^0||_2^2

    // create the argument list (x_i - x_i^0)^2 for all i's
    for (int i=0; i<problem_ -> nVars (); ++i, ++iS) {

      if (problem_ -> Var (i) -> Multiplicity () <= 0)
	continue;

      if (compDistInt_ == FP_DIST_INT &&
	  !(problem_ -> Var (i) -> isInteger ()))
	continue;

      expression *base;

      if      (*iS == 0.) base =              new exprClone (problem_ -> Var (i));
      else if (*iS <  0.) base = new exprSum (new exprClone (problem_ -> Var (i)), new exprConst (-*iS));
      else                base = new exprSub (new exprClone (problem_ -> Var (i)), new exprConst  (*iS));

      list [nTerms++] = new exprPow (base, new exprConst (2.));
    }

  } else {

    // possibly a much larger set of operands

    list = new expression * [problem_ -> nVars () *
			     problem_ -> nVars ()];

    // here the objective function is
    //
    // ||P(x-x^0)||_2^2 = (x-x^0)' P'P (x-x^0)
    //
    // with P'P positive semidefinite stored in CouenneTNLP::optHessian_

    // P is a convex combination, with weights multDistMILP_ and
    // multHessMILP_, of the distance and the Hessian respectively

    int    *row = nlp_ -> optHessian () -> row ();
    int    *col = nlp_ -> optHessian () -> col ();
    double *val = nlp_ -> optHessian () -> val ();

    int     num = nlp_ -> optHessian () -> num ();

    double
      trace_H = 0,
      nActualTerms = 0;

    // create the argument list (x_i - x_i^0)^2 for all i's
    for (int i=0; i<problem_ -> nVars (); ++i)
      if (!((problem_ -> Var (i) -> Multiplicity () <= 0) ||
	    (compDistInt_ == FP_DIST_INT &&
	     !(problem_ -> Var (i) -> isInteger ()))))
	nActualTerms += 1;

    nActualTerms = (nActualTerms == 0) ? 1 : (1 / sqrt (nActualTerms));

    // Add Hessian part -- only lower triangular part
    for (int i=0; i<num; ++i, ++val)
      if (*row++ == *col++)
	trace_H += *val * *val;

    trace_H = (trace_H < COUENNE_EPS) ? 1 : (1 / sqrt (trace_H));

    row = nlp_ -> optHessian () -> row ();
    col = nlp_ -> optHessian () -> col ();
    val = nlp_ -> optHessian () -> val ();

    // Add Hessian part -- only lower triangular part
    for (int i=0; i<num; ++i, ++row, ++col, ++val) {

      if ((problem_ -> Var (*row) -> Multiplicity () <= 0) ||
	  (problem_ -> Var (*col) -> Multiplicity () <= 0))
	continue;

      // check if necessary given options

      if (compDistInt_ == FP_DIST_INT &&
	  !(problem_ -> Var (*row) -> isInteger () &&
	    problem_ -> Var (*col) -> isInteger ()))
	continue;

      // second, only do subdiagonal elements

      if (*col < *row) { // that is, lower triangular

	if (2. * *val * trace_H == 1.) // check if this would have trivial coefficient when doubled (off-diagonal element)

	  list [nTerms++] = new exprMul (new exprSub (new exprClone (problem_ -> Var (*row)), new exprConst (iSol [*row])),
					 new exprSub (new exprClone (problem_ -> Var (*col)), new exprConst (iSol [*col])));

	else if (fabs (*val * trace_H) > COUENNE_EPS) { // we don't need extreme precision...

	  expression **mlist = new expression * [3];

	  mlist [0] = new exprConst (2. * *val * trace_H);  // twice elements off diagonal
	  mlist [1] = new exprSub (new exprClone (problem_ -> Var (*row)), new exprConst (iSol [*row]));
	  mlist [2] = new exprSub (new exprClone (problem_ -> Var (*col)), new exprConst (iSol [*col]));

	  list [nTerms++] = new exprMul (mlist, 3);
	}

      } else if (*col == *row) { // or diagonal elements

	//if (multDistNLP_ != 0.)
	//diag [*col] = true;

	if (trace_H * *val + multDistNLP () * nActualTerms == 1.) // the plus is for the distance term

	  list [nTerms++] = new exprPow (new exprSub (new exprClone (problem_ -> Var (*row)),
						      new exprConst (iSol [*row])),
					 new exprConst (2.));

	else if (fabs (trace_H * *val + nActualTerms * multDistNLP ()) > COUENNE_EPS)

	  list [nTerms++] = new exprMul (new exprConst (trace_H * *val + nActualTerms * multDistNLP ()),
					 new exprPow (new exprSub (new exprClone (problem_ -> Var (*row)),
								   new exprConst (iSol [*row])),
						      new exprConst (2.)));
      }
    }

    // third, add missing diagonal elements. NO! Already added above (fewer terms)
    /*
    if (multDistNLP () > 0.) {

      // create the argument list (x_i - x_i^0)^2 for all i's
      for (int i=0; i<problem_ -> nVars (); ++i, ++iS) {

	if (problem_ -> Var (i) -> Multiplicity () <= 0)
	  continue;
	
	if ((compDistInt_ == FP_DIST_INT &&
	     !(problem_ -> Var (i) -> isInteger ())) ||
	    diag [i])
	  continue;
	
	expression *base;

	if      (*iS == 0.) base =              new exprClone (problem_ -> Var (i));
	else if (*iS <  0.) base = new exprSum (new exprClone (problem_ -> Var (i)), new exprConst (-*iS));
	else                base = new exprSub (new exprClone (problem_ -> Var (i)), new exprConst  (*iS));

	base = new exprMul (base, new exprConst (nActualTerms));

	list [nTerms++] = new exprPow (base, new exprConst (2.));
      }

      delete [] diag;
    } */
  }

  // as per latest development: objective is multiplied by one to
  // normalize it with distance and Hessian-based objective

  if (multObjFNLP () != 0.)
    list [nTerms++] = new exprMul (new exprConst (multObjFNLP ()),
				   new exprClone (problem_ -> Obj (0) -> Body ()));

  // resize list

  expression **tmp = list;
  list = CoinCopyOfArray (tmp, nTerms);
  delete [] tmp;

  expression *retexpr = new exprSum (list, nTerms);

  // printf ("new objective: "); retexpr -> print (); printf ("\n");

  return retexpr;
}


/// Reads a (possibly fractional) solution and fixes the integer
/// components in the nonlinear problem for later re-solve
bool CouenneFeasPump::fixIntVariables (const double *sol) {

  assert (sol);

  t_chg_bounds *chg_bds = new t_chg_bounds [problem_ -> nVars ()];

  for (int i = problem_ -> nVars (); i--;)

    if ((problem_ -> Var (i) -> isInteger ()) &&
	(problem_ -> Var (i) -> Multiplicity () > 0)) {

      double
	value = sol [i],
	rUp   = ceil  (value - COUENNE_EPS),
	rDn   = floor (value + COUENNE_EPS);

      // If numerics or sol[i] fractional, set to closest

      value =
	(rUp < rDn + 0.5)           ? rUp :
	(rUp - value < value - rDn) ? rUp : rDn;

#define INT_NLP_BRACKET 1e-6

      problem_ -> Lb (i) = value - INT_NLP_BRACKET;
      problem_ -> Ub (i) = value + INT_NLP_BRACKET;

      chg_bds [i].setLower (t_chg_bounds::CHANGED);
      chg_bds [i].setUpper (t_chg_bounds::CHANGED);
    }

  // Now, to restrict the bounding box even more (and hopefully make
  // it easier) apply BT

  bool retval = problem_ -> btCore (chg_bds); // maybe fixing makes the nlp infeasible

  delete [] chg_bds;

  return retval;
}


/// initialize options
void CouenneFeasPump::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption4
    ("feas_pump_heuristic",
     "Apply the nonconvex Feasibility Pump",
     "no",
     "no",   "never called",
     "yes",  "called any time Cbc calls heuristics",
     "once", "call it at most once",
     "only", "Call it exactly once and then exit",
     "An implementation of the Feasibility Pump for nonconvex MINLPs");

  roptions -> AddBoundedNumberOption
    ("feas_pump_fademult",
     "decrease/increase rate of multipliers",
     0, false,
     1, false,
     1, "1 keeps initial multipliers from one call to the next; any <1 multiplies ALL of them");

  roptions -> AddLowerBoundedIntegerOption
    ("feas_pump_level",
     "Specify the logarithm of the number of feasibility pumps to perform"
     " on average for each level of given depth of the tree.",
     -1,
     3, "Solve as many nlp's at the nodes for each level of the tree. "
     "Nodes are randomly selected. If for a "
     "given level there are less nodes than this number nlp are solved for every nodes. "
     "For example, if parameter is 8 NLPs are solved for all node until level 8, "
     "then for half the node at level 9, 1/4 at level 10.... "
     "Set to -1 to perform at all nodes.");

  roptions -> AddLowerBoundedIntegerOption
    ("feas_pump_iter",
     "Number of iterations in the main Feasibility Pump loop (default: 10)",
     -1,
     10, "-1 means no limit");

  // six options

  char option [40];
  char help   [250];

  std::string terms [] = {"dist", "hess", "objf"};
  std::string types [] = {"nlp",  "milp"};

  for   (int j=0; j<3; j++)
    for (int i=0; i<2; i++) {

      sprintf (option, "feas_pump_mult_%s_%s",                                              terms [j].c_str (), types [i].c_str ());
      sprintf (help,       "Weight of the %s in the distance function of the %s problem",
	       !(strcmp ("dist", terms [j].c_str ())) ? "distance" :
	       !(strcmp ("hess", terms [j].c_str ())) ? "Hessian"  : "original objective function",             types [i].c_str ());

      roptions -> AddBoundedNumberOption
	(option, help,
	 -1., true,
	 +1., false,
	 0., "0: neglected; 1: full weight; a in ]0,1[: weight is a^k where k is the FP iteration; a in ]-1,0[: weight is 1-|a|^k");
    }

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
     "Number of rounds of convexification cuts. Must be at least 1",
     1, 1e5, 4,
     "");

  roptions -> AddStringOption4
    ("feas_pump_tabumgt",
     "Retrieval of MILP solutions when the one returned is unsatisfactory",
     "pool",
     "pool",       "Use a solution pool and replace unsatisfactory solution with Euclidean-closest in pool",
     "perturb",    "Randomly perturb unsatisfactory solution",
     "cut",        "Separate convexification cuts",
     "none",       "Bail out of feasibility pump");

  roptions -> AddStringOption2
    ("feas_pump_usescip",
     "Should SCIP be used to solve the MILPs?",
#ifdef COUENNE_HAS_SCIP
     "yes", // we want it by default if SCIP is available
#else
     "no",  // otherwise switch it off and warn the user if turned on
#endif
     "no",  "Use Cbc's branch-and-cut to solve the MILP",
     "yes", "Use SCIP's branch-and-cut or heuristics (see feas_pump_milpmethod option) to solve the MILP",
     "");

  roptions -> AddBoundedIntegerOption
    ("feas_pump_milpmethod",
     "How should the integral solution be constructed?",
     -1, 6, 0,
       "0: automatic, 1: aggressive heuristics, large node limit, 2: default, node limit, 3: RENS, 4: Objective Feasibility Pump, 5: MINLP rounding heuristic, 6: rounding, -1: solve MILP completely");

  roptions -> AddBoundedIntegerOption
    ("feas_pump_poolcomp",
     "Priority field to compare solutions in FP pool",
     0, 4, 4,
       "\
0: total number of infeasible objects (integer and nonlinear); \
1: maximum infeasibility (integer or nonlinear); \
2: objective value; \
3: compare value of all variables; \
4: compare value of all integers (RECOMMENDED).");
}
