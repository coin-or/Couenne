/* $Id: CouenneSolverInterface.hpp 154 2009-06-16 18:52:53Z pbelotti $ */
/*
 * Name:    CouenneSolverInterface.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: OsiSolverInterface with a pointer to a CouenneCutGenerator object
 *
 * (C) Carnegie-Mellon University, 2007-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNESOLVERINTERFACE_HPP
#define COUENNESOLVERINTERFACE_HPP

#include "CoinFinite.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CouenneCutGenerator.hpp"

class ClpSimplex;

/// Solver interface class with a pointer to a Couenne cut
/// generator. Its main purposes are:
///
/// 1) to apply bound tightening before re-solving
/// 2) to replace OsiSolverInterface::isInteger () with problem_ -> [expression] -> isInteger ()
/// 3) to use NLP solution at branching
 
class CouenneSolverInterface: public OsiClpSolverInterface {

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
    if (cutgen_ && !(cutgen_ -> enableLpImpliedBounds ()))
      specialOptions_ = specialOptions_ | 262144; 
  }

  /// Solve initial LP relaxation 
  virtual void initialSolve (); 

  /// Resolve an LP relaxation after problem modification
  virtual void resolve ();

  /// Resolve an LP without applying bound tightening beforehand
  virtual void resolve_nobt ()
  {OsiClpSolverInterface::resolve ();}

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

  /// Copy of the Clp version --- not light version
  virtual int tightenBoundsCLP (int lightweight);

  /// Copy of the Clp version --- light version
  virtual int tightenBoundsCLP_Light (int lightweight);

  /// Returns pointer to CLP structure.
  ClpSimplex *continuousModel ()
  {return continuousModel_;}

  /// set doingResolve_
  bool &doingResolve () 
  {return doingResolve_;}

  /// is this problem unbounded?
  bool isProvenDualInfeasible () const;
  //{return knowDualInfeasible_;}

protected:

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

  /// flag to indicate this is an LP for the BB, not for (e.g.) strong
  /// branching or OBBT
  bool doingResolve_;
};

#endif
