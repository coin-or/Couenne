/* $Id$
 *
 * Name:    CouenneComplObject.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Branching object for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECOMPLOBJECT_HPP
#define COUENNECOMPLOBJECT_HPP

#include "CouenneObject.hpp"

/// OsiObject for complementarity constraints $x_1 x_2 <=/>=/= 0$. 
///
/// Associated with two variables x_1 and x_2, branches with either x_1=0 or x_2=0

class CouenneComplObject: public CouenneObject {

public:

  /// Constructor with information for branching point selection strategy
  CouenneComplObject (CouenneCutGenerator *c,
		      CouenneProblem *p, 
		      exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst,
		      int sign);

  /// Constructor with lesser information, used for infeasibility only
  CouenneComplObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst,
		      int sign);

  /// Destructor
  ~CouenneComplObject () {}

  /// Copy constructor
  CouenneComplObject (const CouenneComplObject &src);

  /// Cloning method
  virtual CouenneObject * clone () const
  {return new CouenneComplObject (*this);}

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// compute infeasibility of this variable, |w - f(x)|, where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double checkInfeasibility (const OsiBranchingInformation * info) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject *createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, 
					    int way) const;
private:

  /// -1 if object is for xi * xj <= 0
  /// +1 if object is for xi * xj <= 0
  ///  0 if object is for xi * xj  = 0 (classical)
  int sign_;

};

#endif
