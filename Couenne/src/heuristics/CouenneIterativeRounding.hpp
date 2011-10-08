// (C) Copyright International Business Machines Corporation 2007 
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Giacomo Nannicini, Tepper School of Business, Carnegie Mellon University
//
// Date : 07/25/2010

#ifndef COUENNEITERATIVEROUNDING_HPP
#define COUENNEITERATIVEROUNDING_HPP

#include "CouenneConfig.h"
#include "CbcHeuristic.hpp"
#include "BonOsiTMINLPInterface.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#else
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicRandRound.hpp"
#endif

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

namespace Couenne{

  /** An iterative rounding heuristic, tailored for nonconvex
      MINLPs. It solves a sequence of MILPs and NLPs for a given
      number of iterations, or until a better solution is found.
  */

  class CouenneIterativeRounding : public CbcHeuristic{

  public:
    /** Default constructor.*/
    CouenneIterativeRounding();
    /** Constructor with model and Couenne problems.*/
    CouenneIterativeRounding(Bonmin::OsiTMINLPInterface* nlp, 
			     OsiSolverInterface* cinlp,
			     CouenneProblem* couenne,
			     Ipopt::SmartPtr<Ipopt::OptionsList> options);
    /** Copy constructor.*/
    CouenneIterativeRounding(const CouenneIterativeRounding &other);
    
    /** Destructor*/
    virtual ~CouenneIterativeRounding();
    
    /** Clone.*/
    virtual CbcHeuristic * clone() const;
    
    /** Assignment operator */
    CouenneIterativeRounding & operator=(const CouenneIterativeRounding &rhs);
    
    /** Set the minlp solver.*/
    void setNlp (Bonmin::OsiTMINLPInterface* nlp, OsiSolverInterface* cinlp);

    /** Set the couenne problem to use.*/
    void setCouenneProblem(CouenneProblem* couenne){
      couenne_ = couenne;
    }

    /** Does nothing. */
    void resetModel(CbcModel * model){}
    /** Run heuristic, return 1 if a better solution than the one
        passed is found and 0 otherwise. 
	\argument objectiveValue
        Best known solution in input and value of solution found in
        output 
	\argument newSolution 
	Solution found by heuristic.
      */
    int solution(double & objectiveValue, double * newSolution);

    /** Set maximum number of iterations for each rounding phase */
    void setMaxRoundingIter(int value){
      maxRoundingIter_ = value;
    }
    
    /** Set maximum number of points that we try to round in F-IR */
    void setMaxFirPoints(int value){
      maxFirPoints_ = value;
    }

    /** Set maximum CPU time for the heuristic at each node */
    void setMaxTime(double value){
      maxTime_ = value;
    }

    /** Set maximum CPU time for the heuristic at the root node only*/
    void setMaxTimeFirstCall(double value){
      maxTimeFirstCall_ = value;
    }

    /** Set the value for omega, the multiplicative factor for the minimum
	log-barrier parameter mu used by F-IR whenever we need to generate
	a new NLP feasible point (in the interior of the feasible region)
     */
    void setOmega(double value){
      omega_ = value;
    }

    /** Set the base value for the rhs of the local branching constraint
	in the I-IR heuristic. The actual rhs is then computed depending
	on current variable bounds */
    void setBaseLbRhs(int value){
      baseLbRhs_ = value;
    }

    /** Set aggressiveness of heuristic. Three levels, that sets
	various parameters accordingly.

	The levels are:
	0: maxRoundingIter = 5, maxTimeFirstCall = 300, maxFirPoints = 5,
	   maxTime = 60
	1: maxRoundingIter = 10, maxTimeFirstCall = 300, maxFirPoints = 5,
	   maxTime = 120
	2: maxRoundingIter = 20, maxTimeFirstCall = 1000, maxFirPoints = 5,
	   maxTime = 300
    */
    void setAggressiveness(int value);

    /// initialize options to be read later
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions>);

private:
    /** Pointer to a dynamic nlp solver interface.*/
    Bonmin::OsiTMINLPInterface * nlp_;
    /** Pointer to the original NLP solver interface*/
    OsiSolverInterface * cinlp_;
    /** Pointer to a milp solver interface.*/
#ifdef COIN_HAS_CPX
    OsiCpxSolverInterface * milp_;
#else
    OsiClpSolverInterface * milp_;
#endif
    /** Maximum number of iterations in the main loop*/
    int maxRoundingIter_;
    /** Maximum number of iterations in the outer loop for feasibility*/
    int maxFirPoints_;
    /** Max CPU time to run the heuristic */
    double maxTime_;
    /** Max CPU time to run the heuristic when no other solution is known */
    double maxTimeFirstCall_;
    /** Number of rows in the original convexification */
    int numInitialRows_;
    /** Number of solutions last time the heuristic was called*/
    int numSol_;
    /** Number of integer variables in the original model */
    int numIntegers_;
    /** Pointer to original column lower and upper bounds
        of the reformulated problem, i.e. the linearization */
    double* colLower_;
    double* colUpper_;
    /** Pointer to column lower and upper bounds
        of the original problem, i.e. the MINLP */
    double* colLowerNlp_;
    double* colUpperNlp_;
#ifndef COIN_HAS_CPX
    /** Heuristics for the MILP */
    CbcHeuristic** heuristics_;
    int numHeuristics_;
#endif
    /** Starting time for the heuristics */
    double startTime_;
    /** Maximum allowed time for current run */
    double endTime_;
    /** Multiplication factor for log barrier parameter in F-IR; see \omega
	in the paper */
    double omega_;
    /** Base value for the rhs of the local branching constraint it I-IR */
    int baseLbRhs_;

    /** Pointer to a couenne representation of the problem. */
    CouenneProblem * couenne_;

    /** Check if two double precision numbers are equal, up to some tolerance*/
    inline bool areEqual(double a, double b){
      return (fabs(a-b) <= COUENNE_EPS);
    }

    /* Try to find a first feasible solution */
    int feasibilityIR(double& objectiveValue, double* newSolution);

    /* Try to improve a given solution */
    int improvementIR(double& objectiveValue, double* newSolution,
		      const double* startingSolution);

    /** Set the milp solver at the root */
    void setMilp();

    /** Compute number of integer variables at one of their bounds
	at a given point x */
    int computeIntAtBound(const double* x);

    /** Compute number of integer variables at one of their bounds
	at a given point x and average distance between the bounds of
        variables at one of their bounds */
    int computeIntAtBound(const double* x, double& avgBoundSize);

    /** Write down a local branching constraint, as an OsiRowCut */
    void writeLB(OsiRowCut& cut, const double* x, char sense, double rhs);

    /** Branch on a random variable to cut given solution; returns var index */
    int branchToCut(const double* x, OsiSolverInterface* solver,
		    std::vector<int>& previousBranches);

    /** Solve the MILP contained in milp to feasibility, or report failure */
    bool solveMilp(OsiSolverInterface* milp, double maxTime);
  };
}/* Ends namespace Couenne. */

#endif

