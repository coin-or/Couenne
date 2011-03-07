/* $Id$
 *
 * Name:    CouenneVarObject.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for branching on variables
 *
 * (C) Carnegie-Mellon University, 2008-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEVAROBJECT_HPP
#define COUENNEVAROBJECT_HPP

#include "CouenneObject.hpp"

namespace Couenne {

class CouenneProblem;


/// OsiObject for variables in a MINLP
class CouenneVarObject: public CouenneObject {

public:

  /// Constructor with information for branching point selection strategy
  CouenneVarObject (CouenneCutGenerator *c,
		    CouenneProblem *p, 
		    exprVar *ref, 
		    Bonmin::BabSetupBase *base, 
		    JnlstPtr jnlst,
		    int varSelection); //< either OSI_SIMPLE or OSI_STRONG

  /// Copy constructor
  CouenneVarObject (const CouenneVarObject &src):
    CouenneObject (src),
    varSelection_ (src.varSelection_) {}

  /// Destructor
  ~CouenneVarObject () {}

  /// Cloning method
  virtual CouenneObject *clone () const
  {return new CouenneVarObject (*this);}

  /// compute infeasibility of this variable x as the sum/min/max of
  /// all infeasibilities of auxiliaries w whose defining function
  /// depends on x |w - f(x)|
  ///
  /// TODO: suggest way
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// compute infeasibility of this variable, |w - f(x)|, where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double checkInfeasibility (const OsiBranchingInformation * info) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject *createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

  /// fix nonlinear coordinates of current integer-nonlinear feasible solution
  virtual double feasibleRegion (OsiSolverInterface *, const OsiBranchingInformation *) const;

  /// are we on the bad or good side of the expression?
  virtual bool isCuttable () const;

protected:

  /// branching scheme used. Experimental: still figuring out why
  /// plain LP branching doesn't work with strong/reliability
  /// branching
  int varSelection_;

  /// Method computing the branching point
  CouNumber computeBranchingPoint (const OsiBranchingInformation *info, 
				   int& bestWay, 
				   const CouenneObject *&criticalObject) const;
};

}

#endif
