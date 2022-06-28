/*
 *
 * Name:    CouenneSolverInterface.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: OsiSolverInterface with a pointer to a CouenneCutGenerator object
 *
 * (C) Carnegie-Mellon University, 2007-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNESOLVERINTERFACE_HPP
#define COUENNESOLVERINTERFACE_HPP

class OsiSolverInterface;

namespace Couenne {

class CouenneCutGenerator;

/// Solver interface class with a pointer to a Couenne cut
/// generator. Its main purposes are:
///
/// 1) to apply bound tightening before re-solving
/// 2) to replace OsiSolverInterface::isInteger () with problem_ -> [expression] -> isInteger ()
/// 3) to use NLP solution at branching

template <class T> class COUENNELIB_EXPORT CouenneSolverInterface: public T {

public:

  /// Constructor
  CouenneSolverInterface (CouenneCutGenerator *cg = NULL);

  /// Copy constructor
  CouenneSolverInterface (const CouenneSolverInterface &src);

  /// Destructor
  ~CouenneSolverInterface ();

  /// Clone
  virtual OsiSolverInterface * clone (bool copyData = true) const
  {return new CouenneSolverInterface (*this);}

  /// we need to overwrite this since we might have internal knowledge
  virtual bool isProvenPrimalInfeasible () const;

  /// we need to overwrite this since we might have internal knowledge
  virtual bool isProvenOptimal () const;

  /// Return cut generator pointer
  CouenneCutGenerator *CutGen ()
  {return cutgen_;}

  /// Set cut generator pointer after setup, to avoid changes in the
  /// pointer due to cut generator cloning (it happens twice in the
  /// algorithm)
  void setCutGenPtr (CouenneCutGenerator *cg) {
    cutgen_ = cg;
    //if (cutgen_ && !(cutgen_ -> enableLpImpliedBounds ()))
    //specialOptions_ = specialOptions_ | 262144;
  }

  /// Solve initial LP relaxation
  virtual void initialSolve ();

  /// Resolve an LP relaxation after problem modification
  virtual void resolve ();

  /// Resolve an LP without applying bound tightening beforehand
  virtual void resolve_nobt ()
  {T::resolve ();}

  /** @name Methods for strong branching.
   */
  //@{
  /// Create a hot start snapshot of the optimization process.
  virtual void markHotStart();

  /// Optimize starting from the hot start snapshot.
  virtual void solveFromHotStart();

  /// Delete the hot start snapshot.
  virtual void unmarkHotStart();
  //@}

  /// Tighten bounds on all variables (including continuous).
  virtual int tightenBounds (int lightweight);

  /// set doingResolve_
  //bool &doingResolve ()
  //{return doingResolve_;}

  /// is this problem unbounded?
  bool isProvenDualInfeasible () const;
  //{return knowDualInfeasible_;}

  /// Get the objective function value. Modified due to possible
  /// constant objectives passed to Couenne
  virtual double getObjValue() const;

  // /// returns LP optimum at root node
  // double rootLB () const
  // {return rootLB_;}

protected:

  /// Copy of the Clp version --- not light version
  virtual int tightenBoundsCLP (int lightweight);

  /// Copy of the Clp version --- light version
  virtual int tightenBoundsCLP_Light (int lightweight);

  /// The pointer to the Couenne cut generator. Gives us a lot of
  /// information, for instance the nlp solver pointer, and the chance
  /// to do bound tightening before resolve ().
  CouenneCutGenerator *cutgen_;

  /// Flag indicating that infeasibility was detected during solveFromHotStart
  bool knowInfeasible_;

  /// Flag indicating that optimality was detected during solveFromHotStart
  bool knowOptimal_;

  /// Flag indicating this problem's continuous relaxation is unbounded
  bool knowDualInfeasible_;

  // /// First (root) LP not solved yet -- use this flag to obtain root
  // /// lower bound and print an alternative gap
  // bool beforeFirstRootLP_;

  // /// First (root) LP bound -- used to get an alternative gap
  // double rootLB_;
};

}

// These source files are #included due to the template classes
// defined in there

#include "CouenneSolverInterface.cpp"
#include "CouenneLPtightenBounds.cpp"
#include "CouenneLPtightenBoundsCLP-light.cpp"
#include "CouenneLPtightenBoundsCLP.cpp"

#endif
