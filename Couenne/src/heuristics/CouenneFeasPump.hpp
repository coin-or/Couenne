/* $Id$
 *
 * Name:    CouenneFeasPump.hpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Define the Feasibility Pump heuristic class
 * Created: August 5, 2009
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef CouenneFeasPump_HPP
#define CouenneFeasPump_HPP

#include "CouenneTypes.hpp"
#include "CbcHeuristic.hpp"

#include <queue>

namespace Couenne {

  class expression;
  class CouenneMINLPInterface;
  class CouenneProblem;
  class CouenneCutGenerator;
  class CouenneInterface;

  //const double maxNlpInf_0 = 1e-5;

  /// An implementation of the Feasibility pump that uses
  /// linearization and Ipopt to find the two sequences of points.
  
  class CouenneFeasPump: public CbcHeuristic{

  public:

    // Default constructor
    CouenneFeasPump ();

    /// Constructor with model and Ipopt problems.
    CouenneFeasPump (CbcModel & mip, CouenneMINLPInterface &nlp, 
		     bool cloneNlp = false, CouenneProblem * couenne = NULL);

    /// Copy constructor.
    CouenneFeasPump (const CouenneFeasPump &other);
    
    /// Destructor
    virtual ~CouenneFeasPump();
    
    /// Clone.
    virtual CbcHeuristic *clone () const;
    
    /// Assignment operator 
    CouenneFeasPump &operator= (const CouenneFeasPump &rhs);
    
    /// Set the nlp solver.
    void setNlp (CouenneInterface &nlp, bool cloneNlp = true);
    
    /// set the couenne problem to use.
    void setCouenneProblem (CouenneProblem *);

    /// Does nothing. 
    virtual void resetModel (CbcModel * model) {}

    /// Run heuristic, return 1 if a better solution than the one passed is found and 0 otherwise.
    /// \argument objectiveValue Best known solution in input and value of solution found in output
    /// \argument newSolution Solution found by heuristic.
    /// \todo Find a quicker way to get to Couenne objects, store them or something      
    virtual int solution (double & objectiveValue, double * newSolution);

    /// set maxNlpInf. 
    //void setMaxNlpInf (double value)
    //{maxNlpInf_ = value;}

    /// set number of nlp's solved for each given level of the tree
    void setNumberSolvePerLevel (int value)
    {numberSolvePerLevel_ = value;}

    /// find integer (possibly NLP-infeasible) point isol closest
    /// (according to the l-1 norm of the hessian) to the current
    /// NLP-feasible (but fractional) solution nsol
    CouNumber getMILPSolution (CouNumber *iSol, CouNumber *nSol); 

    /// obtain solution to NLP
    CouNumber *getContSolution ();

    /// set new expression as the NLP objective function using
    /// argument as point to minimize distance from. Return new
    /// objective function
    expression *updateNLPObj (double *);

    /// admits a (possibly fractional) solution and fixes the integer
    /// components in the nonlinear problem for later re-solve
    void fixIntVariables (double *sol);

  private:

    /// Pointer to an nlp solver interface.
    CouenneMINLPInterface *nlp_;

    /// is nlp_ cloned or just a pointer?
    bool hasCloned_;

    /// maximum nlp infeasibility under which try to solve problem with Ipopt.
    //double maxNlpInf_;

    /// Number of nlp's solved for each given level of the tree
    int numberSolvePerLevel_;

    /// Pointer to a couenne representation of the problem. 
    CouenneProblem *problem_;

    /// Pointer to a CouenneCutGenerator for linearization cuts
    CouenneCutGenerator *couenneCG_;

    /// Pool of solutions
    std::priority_queue <std::pair <CouNumber *, CouNumber> > pool_;

    /// Skip NLP solver if found integer but MINLP-infeasible solution 
    bool milpCuttingPlane_;

    /// interface to MILP problem
    OsiSolverInterface *milp_;

    /// maximum iterations per call
    int maxIter_;

    /// original objective function (to be restored at the end)
    expression *originalObjective_;
  };
}

#endif
