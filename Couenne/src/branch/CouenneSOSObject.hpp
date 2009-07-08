/* $Id$
 *
 * Name:    CouenneSOSObject.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: SOS Object
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNESOSOBJECT_HPP
#define COUENNESOSOBJECT_HPP

#include "OsiBranchingObject.hpp"
#include "CouenneJournalist.hpp"

class CouenneProblem;
class CouenneSOSObject;
class exprVar;


// TODO: SOS of the form sum x_i \le k with k small. Instead of
// branching on a single variable do a SOS-like branching

class CouenneSOSBranchingObject: public OsiSOSBranchingObject {

protected:

  /// pointer to Couenne problem
  CouenneProblem *problem_;

  /// The (auxiliary) variable this branching object refers to. If the
  /// expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprVar *reference_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// shall we do Feasibility based Bound Tightening (FBBT) at branching?
  bool doFBBT_;

  /// shall we add convexification cuts at branching?
  bool doConvCuts_;

public:

  // Default Constructor 
  CouenneSOSBranchingObject () {}

  // Useful constructor
  CouenneSOSBranchingObject (CouenneProblem *p,
			     exprVar *ref,
			     OsiSolverInterface *solver,  
			     const OsiSOS *originalObject,
			     int way, 
			     double separator,
			     JnlstPtr jnlst,
			     bool doFBBT,
			     bool doConvCuts):

    OsiSOSBranchingObject (solver, originalObject, way, separator),
    problem_   (p),
    reference_ (ref),
    jnlst_ (jnlst),
    doFBBT_ (doFBBT),
    doConvCuts_ (doConvCuts) {}

  
  // Copy constructor 
  CouenneSOSBranchingObject (const CouenneSOSBranchingObject &src):
    OsiSOSBranchingObject (src),
    problem_    (src.problem_),
    reference_  (src.reference_),
    jnlst_      (src.jnlst_),
    doFBBT_     (src.doFBBT_),
    doConvCuts_ (src.doConvCuts_) {}

   
  /// Clone
  virtual OsiBranchingObject * clone() const
  {return new CouenneSOSBranchingObject (*this);}

  /// Does next branch and updates state
  virtual double branch (OsiSolverInterface * solver);
};


///
///
///

class CouenneSOSObject: public OsiSOS {

protected:

  /// pointer to Couenne problem
  CouenneProblem *problem_;

  /// The (auxiliary) variable this branching object refers to. If the
  /// expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprVar *reference_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// shall we do Feasibility based Bound Tightening (FBBT) at branching?
  bool doFBBT_;

  /// shall we add convexification cuts at branching?
  bool doConvCuts_;

public:

  CouenneSOSObject (OsiSolverInterface *solver, int nelem, int *indices, double *weights, int type,
		    CouenneProblem *problem,
		    exprVar *ref,
		    JnlstPtr jnlst,
		    bool doFBBT,
		    bool doConvCuts):

    OsiSOS      (solver, nelem, indices, weights, type),
    problem_    (problem),
    reference_  (ref),
    jnlst_      (jnlst),
    doFBBT_     (doFBBT),
    doConvCuts_ (doConvCuts) {}


  /// Copy constructor
  CouenneSOSObject (const CouenneSOSObject &src):
    OsiSOS      (src),
    problem_    (src.problem_),
    reference_  (src.reference_),
    jnlst_      (src.jnlst_),
    doFBBT_     (src.doFBBT_),
    doConvCuts_ (src.doConvCuts_) {}
    
  /// Cloning method
  virtual OsiObject * clone () const
  {return new CouenneSOSObject (*this);}

  /// create branching objects
  OsiBranchingObject *createBranch (OsiSolverInterface* si, 
				    const OsiBranchingInformation* info, int way) const;

  /// return reference auxiliary variable
  //exprVar *Reference () const
  //{return reference_;}
};

#endif
