/* $Id$
 *
 * Name:    CouenneNlpRoundOne.hpp
 * Author:  Pietro Belotti
 * Purpose: Definition of class for rounding heuristic that simply
 *          calls a NLP solver on a (possibly fractional) solution to
 *          obtain an integer feasible one
 */

#ifndef CouenneNlpRoundOne_hpp
#define CouenneNlpRoundOne_hpp

#include "BonOsiTMINLPInterface.hpp"
#include "CbcHeuristic.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "CouenneProblem.hpp"

namespace Couenne {

  // /** A heuristic to call an NlpSolver if all CouenneObjects are
  //     close to be satisfied (for other integer objects, rounding is
  //     performed, is SOS are not satisfied does not run).  */

  // const double maxNlpInf_0 = 1e-5;

  class CouenneNlpRoundOne: public CbcHeuristic {

  public:

    CouenneNlpRoundOne ();                                                                                                       ///< Empty constructor
    CouenneNlpRoundOne (CbcModel &mip, Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp = false, CouenneProblem * couenne = NULL); ///< Constructor with model and Ipopt problems.
    CouenneNlpRoundOne (const CouenneNlpRoundOne &other);                                                                        ///< Copy constructor
    virtual ~CouenneNlpRoundOne ();                                                                                              ///< Destructor
    virtual CbcHeuristic * clone () const;                                                                                       ///< Clone
    CouenneNlpRoundOne & operator= (const CouenneNlpRoundOne &rhs);                                                              ///< Assignment

    void setNlp (Bonmin::OsiTMINLPInterface &nlp, bool cloneNlp = true); ///< set the NLP solver
    void setCouenneProblem(CouenneProblem *);                            ///< set CouenneProblem pointer

    virtual void resetModel (CbcModel * model) {}

    /** Run heuristic

	Return 1 if a better solution than the one passed is found and 0 otherwise.
        \argument objectiveValue Best known solution in input and value of solution found in output
        \argument newSolution Solution found by heuristic.
    */
    virtual int solution (double &objectiveValue, double * newSolution);

    /// initialize options
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions>);

    /** set number of NLPs solved for each given level of the tree*/
    void setNumberSolvePerLevel (int value)
    {numberSolvePerLevel_ = value;}

  private:

    Bonmin::OsiTMINLPInterface * nlp_; ///< Pointer to an NLP solver interface
    bool hasCloned_;                   ///< Is nlp_ cloned or just a pointer?
    CouenneProblem * couenne_;         ///< Pointer to CouenneProblem
    int numberSolvePerLevel_;          ///< Number of nlp's solved for each given level of the tree
  };
}

#endif
