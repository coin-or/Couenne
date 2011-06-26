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

struct Scip;
class OsiSolverInterface;

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
  class CouenneSparseMatrix;

  /// An implementation of the Feasibility pump that uses
  /// linearization and Ipopt to find the two sequences of points.
  
  class CouenneFeasPump: public CbcHeuristic {

  public:

    enum fpCompDistIntType {FP_DIST_INT,  FP_DIST_ALL,       FP_DIST_POST};
    enum fpCutPlane        {FP_CUT_NONE,  FP_CUT_INTEGRATED, FP_CUT_EXTERNAL, FP_CUT_POST}; 
    enum fpTabuMgtPolicy   {FP_TABU_NONE, FP_TABU_POOL,      FP_TABU_PERTURB, FP_TABU_CUT};

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
    void checkInfinity (struct Scip *scip, double val, double infinity);
#endif

    /// find integer (possibly NLP-infeasible) point isol closest
    /// (according to the l-1 norm of the hessian) to the current
    /// NLP-feasible (but fractional) solution nsol
    CouNumber solveMILP (CouNumber *nSol, CouNumber *&iSol, int niter, int* nsuciter); 

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
    double findSolution (double *&sol, int niter, int* nsuciter);

    /// initialize all solvers at the first call, where the initial
    /// MILP is built
    void init_MILP ();

    /// Common code for initializing non-smartptr ipopt application
    void initIpoptApp ();

    /// return pointer to problem
    CouenneProblem *Problem () const
    {return problem_;}

    /// return type of MILP solved
    enum fpCompDistIntType compDistInt () const
    {return compDistInt_;}

    /// Return Weights in computing distance, in both MILP and NLP (must sum
    /// up to 1 for MILP and for NLP):

    double multDistNLP  () const {return multDistNLP_;}  ///< weight of distance  in NLP
    double multHessNLP  () const {return multHessNLP_;}  ///< weight of Hessian   in NLP
    double multObjFNLP  () const {return multObjFNLP_;}  ///< weight of objective in NLP

    double multDistMILP () const {return multDistMILP_;} ///< weight of distance  in MILP
    double multHessMILP () const {return multHessMILP_;} ///< weight of Hessian   in MILP
    double multObjFMILP () const {return multObjFMILP_;} ///< weight of objective in MILP

    /// return NLP
    CouenneTNLP *nlp () const
    {return nlp_;}

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
    // PERSISTENT OBJECTS
    //
    // (not necessary to identify FP, but useful to keep between
    // calls)
    //

    /// Continuous relaxation of the problem, with an interface for
    /// Ipopt only
    CouenneTNLP *nlp_;

    /// Ipopt Application pointer for solving NLPs
    Ipopt::IpoptApplication *app_;

    /// MILP relaxation of the MINLP (used to find integer
    /// non-NLP-feasible solution)
    OsiSolverInterface *milp_;

    /// LP relaxation of the MINLP used when fixing integer variables
    /// (used for compDistInt_ in FP_DIST_POST and possibly
    /// FP_DIST_INT)
    OsiSolverInterface *postlp_;

    /// Pool of solutions
    CouenneFPpool *pool_;

    /// Solutions to avoid
    std::set <CouenneFPsolution, compareSol> tabuPool_;

    //
    // PARAMETERS
    //

    /// Number of NLPs solved for each given level of the tree
    int numberSolvePerLevel_;

    /// Weights in computing distance, in both MILP and NLP (must sum
    /// up to 1 for MILP and for NLP):

    double multDistNLP_;  ///< weight of distance  in NLP
    double multHessNLP_;  ///< weight of Hessian   in NLP
    double multObjFNLP_;  ///< weight of objective in NLP

    double multDistMILP_; ///< weight of distance  in MILP
    double multHessMILP_; ///< weight of Hessian   in MILP
    double multObjFMILP_; ///< weight of objective in MILP

    /// compute distance from integer variables only, not all variables;
    enum fpCompDistIntType compDistInt_;

    /// Separate convexification cuts during or after MILP
    enum fpCutPlane milpCuttingPlane_;

    /// Number of separation rounds for MILP convexification cuts
    int nSepRounds_;

    /// maximum iterations per call
    int maxIter_;

    /// use SCIP instead of Cbc for solving MILPs
    bool useSCIP_;

    /// Which SCIP MILP method to use
    int milpMethod_;

    /// Tabu management policy: none, use from pool, random perturbation of current solution
    enum fpTabuMgtPolicy tabuMgt_;
  };
}

#endif
