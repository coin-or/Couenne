/*
 *
 * Name:    CouenneChooseStrong.hpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching object for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNECHOOSESTRONG_HPP
#define COUENNECHOOSESTRONG_HPP

#include "BonChooseVariable.hpp"
#include "CouenneJournalist.hpp"
#include "CouenneConfig.h"

namespace Couenne {

class CouenneProblem;

template <class T> class CouenneSolverInterface;

class COUENNELIB_EXPORT CouenneChooseStrong : public Bonmin::BonChooseVariable {

public:

  /// Constructor from solver (so we can set up arrays etc)
  CouenneChooseStrong (Bonmin::BabSetupBase& b, CouenneProblem* problem, JnlstPtr jnlst);

  /// Copy constructor
  CouenneChooseStrong (const CouenneChooseStrong &);

  /// Assignment operator
  CouenneChooseStrong & operator= (const CouenneChooseStrong& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~CouenneChooseStrong ();

  /// Sets up strong list and clears all if initialize is true.
  /// Returns number of infeasibilities.
  virtual int setupList (OsiBranchingInformation *info, bool initialize);

  // actually setting up the list
  int gutsOfSetupList(OsiBranchingInformation *info, bool initialize);

  /**  This is a utility function which does strong branching on a
       list of objects and stores the results in OsiHotInfo.objects.
       On entry the object sequence is stored in the OsiHotInfo
       object and maybe more.
       It returns -
       -1 - one branch was infeasible both ways
       0 - all inspected - nothing can be fixed
       1 - all inspected - some can be fixed                         (returnCriterion==0)
       2 - may be returning early - one can be fixed (last one done) (returnCriterion==1)
       3 - returning because max time
  */
  virtual int doStrongBranching (OsiSolverInterface * OsiSolver,
				 OsiBranchingInformation *info,
				 int numberToDo, int returnCriterion);

  /// Returns true if solution looks feasible against given objects
  virtual bool feasibleSolution (const OsiBranchingInformation * info,
				 const double * solution,
				 int numberObjects,
				 const OsiObject ** objects);

  /// choose object to branch based on earlier setup
  virtual int chooseVariable (OsiSolverInterface * solver,
			      OsiBranchingInformation *info,
			      bool fixVariables);

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

private:

  /** Default Constructor, forbidden for some reason.*/
  CouenneChooseStrong ();

  // print object violations
  void printObjViol(OsiBranchingInformation *info);

  // Due to possible fixing during strong branching, must check
  // that variable reference vInd for the selected object (if not -1)
  // is not fixed or has a current value inside bounds.
  // Return value is:
  // 0: OsiSimpleInteger with upper[vInd] == lower[vInd]
  // 1: OsiSimpleInteger with upper[vInd] > lower[vInd] and
  //    info->solution[vInd] outside bounds
  // 2: CouenneBranchingObject with upper[vInd] == lower[vInd]
  // 3: otherwise (meaning good object)
  int goodCandidate(OsiSolverInterface *solver,
		    OsiBranchingInformation *info,
		    OsiObject **object, const int iObject, const double prec);

  /// Save best candidate
  bool saveBestCand(OsiObject **object, const int iObject,
		    const double value,
		    const double upEstimate,
		    const double downEstimate,
		    double &bestVal1,
		    double &bestVal2, int &bestIndex,
		    int &bestWay);

protected:

  /// does one side of the branching
  int simulateBranch (OsiObject *Object,
		      OsiBranchingInformation *info,
		      OsiBranchingObject *branch,
		      OsiSolverInterface *solver,
		      Bonmin::HotInfo * result,
		      int direction);

  /// Pointer to the associated MINLP problem
  CouenneProblem *problem_;

  /// should we update the pseudocost multiplier with the distance
  /// between the LP point and the solution of the resulting
  /// branches' LPs? If so, this only happens in strong branching
  bool pseudoUpdateLP_;

  /// Normally, a convex combination of the min/max lower bounds'
  /// estimates is taken to select a branching variable, as in the
  /// original definition of strong branching. If this option is set
  /// to true, their product is taken instead:
  ///
  /// (1e-6+min) * max
  ///
  /// where the 1e-6 is used to ensure that even those with one
  /// subproblem with no improvement are compared.
  bool estimateProduct_;

  /// pointer to journalist for detailed information
  JnlstPtr jnlst_;

  /// total time spent in strong branching
  double branchtime_;
};

}

#endif
