/* $Id$
 *
 * Name:    CouenneFeasPump.hpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Define the Feasibility Pump heuristic class
 * Created: August 5, 2009
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneFeasPump_HPP
#define CouenneFeasPump_HPP

#include <queue>

#include "CouenneTypes.hpp"
#include "CbcHeuristic.hpp"

#include "CouenneConfig.h"

#include "CouenneFPpool.hpp"

#include "IpOptionsList.hpp"

#ifdef COIN_HAS_SCIP
/* general SCIP includes */
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#endif

namespace Osi {
  class OsiSolverInterface;
}

namespace Ipopt {
  class IpoptApplication;
}

namespace Bonmin {
  class RegisteredOptions;
}

namespace Couenne {

  class expression;
  class CouenneProblem;
  class CouenneCutGenerator;
  class CouenneTNLP;

  /// An implementation of the Feasibility pump that uses
  /// linearization and Ipopt to find the two sequences of points.
  
  class CouenneFeasPump: public CbcHeuristic {

  public:

    /// Constructor with (optional) MINLP pointer
    CouenneFeasPump (CouenneProblem *couenne                     = NULL,
		     CouenneCutGenerator *cg                     = NULL,
		     Ipopt::SmartPtr<Ipopt::OptionsList> options = NULL);

    /// Copy constructor
    CouenneFeasPump (const CouenneFeasPump &other);
    
    /// Destructor
    virtual ~CouenneFeasPump();

    /// Clone
    virtual CbcHeuristic *clone () const;
    
    /// Assignment operator 
    CouenneFeasPump &operator= (const CouenneFeasPump &rhs);

    /// Does nothing, but necessary as CbcHeuristic declares it pure virtual
    virtual void resetModel (CbcModel *model) {}

    /// Run heuristic, return 1 if a better solution than the one
    /// passed is found and 0 otherwise.
    ///
    /// \argument objectiveValue Best known solution in input and
    /// value of solution found in output
    ///
    /// \argument newSolution Solution found by heuristic.
    virtual int solution (double &objectiveValue, double *newSolution);

    /// set number of nlp's solved for each given level of the tree
    void setNumberSolvePerLevel (int value)
    {numberSolvePerLevel_ = value;}

#ifdef COIN_HAS_SCIP
    /// checks if val is above a threshold for finiteness
    void checkInfinity (SCIP *scip, SCIP_Real val, double infinity);
#endif

    /// find integer (possibly NLP-infeasible) point isol closest
    /// (according to the l-1 norm of the hessian) to the current
    /// NLP-feasible (but fractional) solution nsol
    CouNumber solveMILP (CouNumber *nSol, CouNumber *&iSol); 

    /// obtain solution to NLP
    CouNumber solveNLP  (CouNumber *nSol, CouNumber *&iSol);

    /// set new expression as the NLP objective function using
    /// argument as point to minimize distance from. Return new
    /// objective function
    expression *updateNLPObj (const double *);

    /// admits a (possibly fractional) solution and fixes the integer
    /// components in the nonlinear problem for later re-solve
    void fixIntVariables (double *sol);

    /// initialize options to be read later
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions>);

    /// find feasible solution (called by solveMILP ())
    double findSolution (double *&sol);

    /// initialize all solvers at the first call, where the initial
    /// MILP is built
    void init_MILP ();

    /// Common code for initializing non-smartptr ipopt application
    void initIpoptApp ();

  private:

    //
    // Essential tools for the FP: a problem pointer and one for the
    // linearization cut generator
    //

    /// Couenne representation of the problem. 
    CouenneProblem *problem_;

    /// CouenneCutGenerator for linearization cuts
    CouenneCutGenerator *couenneCG_;

    //
    // PERSISTENT OBJECTS -- not necessary to identify FP, but it's
    // useful to keep them between calls
    //

    /// Continuous relaxation of the problem, with an interface for
    /// Ipopt only
    CouenneTNLP *nlp_;

    /// Ipopt Application pointer for solving NLPs
    Ipopt::IpoptApplication *app_;

    /// MILP relaxation of the MINLP (used to find integer
    /// non-NLP-feasible solution)
    OsiSolverInterface *milp_;

    /// Pool of solutions
    CouenneFPpool *pool_;

    /// Solutions to avoid
    std::set <CouenneFPsolution, compareSol> tabuPool_;

    /// These methods can be used to solve the MILP, from the most
    /// expensive, exact method to a cheap, inexact one:
    ///
    /// 1. Solve a MILP relaxation with Manhattan distance as objective
    /// 2. Apply RENS on 1.
    /// 3. Use Objective FP 2.0 for MILPs
    /// 4. round-and-propagate
    /// 5. choose from pool, see 4
    /// 6. random pertubation

    enum Methods {SOLVE_MILP,       APPLY_RENS,  USE_FP_2, ROUND_N_PROPAGATE, 
		  CHOOSE_FROM_POOL, RND_PERTURB, NUM_METHODS};

    /// array of NUM_METHODS elements containing a number indicating
    /// the recent successes of the corresponding algorithm (see enum
    /// right above). The larger method_success_ [i], the more often
    /// the i-th method should be used.

    int *method_success_;

    //
    // PARAMETERS
    //

    /// Number of NLPs solved for each given level of the tree
    int numberSolvePerLevel_;

    /// weight of the Hessian in computing the objective functions of NLP
    double betaNLP_;

    /// weight of the Hessian in computing the objective functions of MILP
    double betaMILP_; // decrease it over time

    /// compute distance from integer variables only, not all variables;
    bool compDistInt_;

    /// Skip NLP solver if found integer but MINLP-infeasible solution 
    bool milpCuttingPlane_;

    /// maximum iterations per call
    int maxIter_;

    /// use SCIP instead of Cbc for solving MILPs
    bool useSCIP_;

    /// Which SCIP MILP method to use
    int milpMethod_;
  };
}

void printCmpSol (int n, double *iSol, double *nSol, int direction);

#endif
