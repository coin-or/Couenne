/*
 *
 * Name:    CouenneBranchingObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEBRANCHINGOBJECT_HPP
#define COUENNEBRANCHINGOBJECT_HPP

#include "stdio.h"
#include "CouenneExprAux.hpp"
#include "CouenneJournalist.hpp"
#include "OsiBranchingObject.hpp"

namespace Couenne {

class CouenneCutGenerator;
class CouenneProblem;

#define COUENNE_CROP 1
#define COUENNE_LCROP (1e2*COUENNE_CROP)

#define COUENNE_LARGE_INTERVAL 1e4
#define COUENNE_NEAR_BOUND 1e-2


/** "Spatial" branching object.
 *
 *  Branching can also be performed on continuous variables.
 */

class COUENNELIB_EXPORT CouenneBranchingObject: public OsiTwoWayBranchingObject {

public:

  /// Constructor
  CouenneBranchingObject (OsiSolverInterface *solver,
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
  CouenneBranchingObject (const CouenneBranchingObject &src):

    OsiTwoWayBranchingObject (src),
    cutGen_       (src.cutGen_),
    problem_      (src.problem_),
    variable_     (src.variable_),
    jnlst_        (src.jnlst_),
    doFBBT_       (src.doFBBT_),
    doConvCuts_   (src.doConvCuts_),
    downEstimate_ (src.downEstimate_),
    upEstimate_   (src.upEstimate_),
    simulate_     (src.simulate_) {}

  /// cloning method
  virtual OsiBranchingObject * clone () const
  {return new CouenneBranchingObject (*this);}

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

  /// return branching variable
  expression *variable ()
  {return variable_;}

  /// Perform branching step
  void branchCore (OsiSolverInterface *, int, int, bool, double, t_chg_bounds *&);

  // FIXME: horrible global variables. Brrr.
  static int nOrbBr;
  static int maxDepthOrbBranch;
  static int nSGcomputations;

protected:

  /// Pointer to CouenneCutGenerator (if any); if not NULL, allows to
  /// do extra cut generation during branching
  CouenneCutGenerator *cutGen_;

  /// Pointer to CouenneProblem (necessary to allow FBBT)
  CouenneProblem *problem_;

  /// The index of the variable this branching object refers to. If
  /// the corresponding CouenneObject was created on w=f(x,y), it is
  /// either x or y, chosen previously with a call to getFixVar()
  /// expression *reference_;
  expression *variable_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// shall we do Feasibility based Bound Tightening (FBBT) at branching?
  bool doFBBT_;

  /// shall we add convexification cuts at branching?
  bool doConvCuts_;

  /// down branch estimate (done at selectBranch with reduced costs)
  double downEstimate_;

  /// up branch estimate
  double upEstimate_;

  /// are we currently in strong branching?
  bool simulate_;
};

}

#endif
