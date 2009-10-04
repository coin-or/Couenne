/* $Id$
 *
 * Name:    CouenneComplBranchingObject.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Branching object for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECOMPLBRANCHINGOBJECT_HPP
#define COUENNECOMPLBRANCHINGOBJECT_HPP

#include "CouenneBranchingObject.hpp"

/** "Spatial" branching object for complementarity constraints. 
 *
 *  Branching on such an object x_1 x_2 = 0 is performed by setting
 *  either x_1=0 or x_2=0
 */

class CouenneComplBranchingObject: public CouenneBranchingObject {

public:

  /// Constructor
  CouenneComplBranchingObject (OsiSolverInterface *solver,
			       const OsiObject *originalObject,
			       JnlstPtr jnlst, 
			       CouenneCutGenerator *c,
			       CouenneProblem *p,
			       expression *var, 
			       expression *var2,
			       int way, 
			       CouNumber brpoint, 
			       bool doFBBT, 
			       bool doConvCuts,
			       int sign);

  /// Copy constructor
  CouenneComplBranchingObject (const CouenneComplBranchingObject &src):
    CouenneBranchingObject (src),
    variable2_ (src.variable2_),
    sign_      (src.sign_) {}

  /// cloning method
  virtual OsiBranchingObject *clone () const
  {return new CouenneComplBranchingObject (*this);}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state. 
	     Returns change in guessed objective on next branch
  */
  virtual double branch (OsiSolverInterface * solver = NULL);

protected:

  /// use CouenneBranchingObject::variable_ as the first variable to set to 0,
  /// and this one as the second
  expression *variable2_;

  /// -1 if object is for xi * xj <= 0
  /// +1 if object is for xi * xj <= 0
  ///  0 if object is for xi * xj  = 0 (classical)
  int sign_;
};

#endif
