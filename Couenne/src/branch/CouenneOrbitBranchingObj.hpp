/* $Id$
 *
 * Name:    CouenneOrbitBranchingObj.hpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 * Purpose: Branching object for auxiliary variables
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEORBITBRANCHINGOBJ_HPP
#define COUENNEORBITBRANCHINGOBJ_HPP

#include "CouenneExprAux.hpp"
#include "CouenneJournalist.hpp"
#include "OsiBranchingObject.hpp"
#include "CouenneBranchingObject.hpp"

namespace Couenne {

//class CouenneCutGenerator;
//class CouenneProblem;

//#define COUENNE_CROP 1
//#define COUENNE_LCROP (1e2*COUENNE_CROP)

//#define COUENNE_LARGE_INTERVAL 1e4
//#define COUENNE_NEAR_BOUND 1e-2


/** "Spatial" branching object. 
 *
 *  Branching can also be performed on continuous variables.
 */

class CouenneOrbitBranchingObj: public CouenneBranchingObject {

public:

  /// Constructor
  CouenneOrbitBranchingObj (OsiSolverInterface *solver,
			    const OsiObject *originalObject,
			    JnlstPtr jnlst, 
			    CouenneCutGenerator *c,
			    CouenneProblem *p,
			    expression *var, 
			    int way, 
			    CouNumber brpoint, 
			    bool doFBBT, 
			    bool doConvCuts);

  /// Copy constructor
  CouenneOrbitBranchingObj (const CouenneOrbitBranchingObj &src):
    CouenneBranchingObject (src) {}


  /// cloning method
  virtual OsiBranchingObject * clone () const
  {return new CouenneOrbitBranchingObj (*this);}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state. 
	     Returns change in guessed objective on next branch
  */
  virtual double branch (OsiSolverInterface * solver = NULL);

  /// does this branching object only change variable bounds?
  virtual bool boundBranch () const
  {return !doConvCuts_;} // iff it does not add convexification cuts

  /// set simulate_ field below
  void setSimulate (bool s)
  {simulate_ = s;}

protected:

};

}

#endif
