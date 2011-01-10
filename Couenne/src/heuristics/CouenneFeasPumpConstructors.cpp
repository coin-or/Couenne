/* $Id$
 *
 * Name:    CouenneFeasPump.cpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the Feasibility Pump heuristic class
 *
 * This file is licensed under the Eclipse Public License (EPL) (EPL)
 */

#include <string>

#include "IpOptionsList.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprSum.hpp"
#include "CoinHelperFunctions.hpp"

using namespace Couenne;

// Constructor ////////////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump ():

  CbcHeuristic         (),
  problem_             (NULL),
  couenneCG_           (NULL),
  nlp_                 (NULL),
  milp_                (NULL),
  numberSolvePerLevel_ (-1),
  betaNLP_             (0.),
  betaMILP_            (0.),
  compDistInt_         (false),
  milpCuttingPlane_    (false),
  maxIter_             (COIN_INT_MAX) {

  setHeuristicName ("Couenne Feasibility Pump");
}


// Constructor ////////////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (CouenneProblem *couenne,
				  CouenneCutGenerator *cg,
				  Ipopt::SmartPtr<Ipopt::OptionsList> options):

  CbcHeuristic         (), //(model),
  problem_             (couenne -> clone ()),
  couenneCG_           (cg),
  nlp_                 (NULL),
  milp_                (NULL),
  numberSolvePerLevel_ (-1),
  betaNLP_             (0.),
  betaMILP_            (0.),
  compDistInt_         (false),
  milpCuttingPlane_    (false),
  maxIter_             (COIN_INT_MAX) {

  setHeuristicName ("Couenne Feasibility Pump");

  std::string s;

  options -> GetIntegerValue ("feas_pump_iter",      maxIter_,             "couenne.");
  options -> GetIntegerValue ("feas_pump_level",     numberSolvePerLevel_, "couenne.");
  options -> GetNumericValue ("feas_pump_beta_nlp",  betaNLP_,             "couenne.");
  options -> GetNumericValue ("feas_pump_beta_milp", betaMILP_,            "couenne.");
				    
  options -> GetStringValue  ("feas_pump_lincut",   s, "couenne."); milpCuttingPlane_ = (s == "yes");
  options -> GetStringValue  ("feas_pump_dist_int", s, "couenne."); compDistInt_      = (s == "yes");
}

  
// Copy constructor ///////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (const CouenneFeasPump &other):

  CbcHeuristic         (other),
  problem_             (other. problem_ -> clone ()),
  couenneCG_           (other. couenneCG_),
  nlp_                 (other. nlp_),
  milp_                (other. milp_),
  pool_                (other. pool_),
  numberSolvePerLevel_ (other. numberSolvePerLevel_),
  betaNLP_             (other. betaNLP_),
  betaMILP_            (other. betaMILP_),
  compDistInt_         (other. compDistInt_),
  milpCuttingPlane_    (other. milpCuttingPlane_),
  maxIter_             (other. maxIter_) {}


// Clone //////////////////////////////////////////////////////// 
CbcHeuristic *CouenneFeasPump::clone () const
{return new CouenneFeasPump (*this);}


// Assignment operator ////////////////////////////////////////// 
CouenneFeasPump &CouenneFeasPump::operator= (const CouenneFeasPump & rhs) {

  if (this != &rhs) {

    CbcHeuristic::operator= (rhs);

    problem_             = rhs. problem_ -> clone ();
    couenneCG_           = rhs. couenneCG_;
    nlp_                 = rhs. nlp_;
    milp_                = rhs. milp_;
    pool_                = rhs. pool_;
    numberSolvePerLevel_ = rhs. numberSolvePerLevel_;
    betaNLP_             = rhs. betaNLP_;
    betaMILP_            = rhs. betaMILP_;
    compDistInt_         = rhs. compDistInt_;
    milpCuttingPlane_    = rhs. milpCuttingPlane_;
    maxIter_             = rhs. maxIter_;
  }

  return *this;
}


// Destructor /////////////////////////////////////////////////// 
CouenneFeasPump::~CouenneFeasPump () {

  if (problem_)
    delete problem_;
}


/// Set new expression as the NLP objective function using
/// argument as point to minimize distance from. Return new
/// objective function.
expression *CouenneFeasPump::updateNLPObj (const double *iSol) {

  expression **list = new expression * [problem_ -> nVars ()];

  if (betaNLP_ == 0.) {

    // here the objective function is ||x-x^0||_2^2

    // create the argument list (x_i - x_i^0)^2 for all i's
    for (int i=0; i<problem_ -> nVars (); i++)
      list [i] = new exprPow (new exprSub (new exprClone (problem_ -> Var (i)), 
					   new exprConst (iSol [i])), 
			      new exprConst (2.));
  } else {

    // here the objective function is ||P(x-x^0)||_2^2 with P positive semidefinite

  }

  return new exprSum (list, problem_ -> nVars ());
}


/// Reads a (possibly fractional) solution and fixes the integer
/// components in the nonlinear problem for later re-solve
void CouenneFeasPump::fixIntVariables (double *sol) {

  problem_ -> domain () -> push (problem_ -> nVars (),
				 problem_ -> domain () -> x (),
				 problem_ -> domain () -> lb (),
				 problem_ -> domain () -> ub ());

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

  roptions -> AddStringOption2
    ("feas_pump_lincut",
     "Shortcut to linearization cutting plane applied to the MILP solution instead of solving NLPs",
     "no",
     "no","",
     "yes","",
     "");

  roptions -> AddStringOption2
    ("feas_pump_dist_int",
     "only compute the distance from integer coordinates (\"yes\") instead of all variables (\"no\")",
     "yes",
     "no","",
     "yes","",
     "");
}
