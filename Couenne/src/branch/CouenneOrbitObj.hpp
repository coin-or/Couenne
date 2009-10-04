/* $Id$
 *
 * Name:    CouenneOrbitObj.hpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 * Purpose: Object for auxiliary variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEORBITOBJ_HPP
#define COUENNEORBITOBJ_HPP

#include "BonBabSetupBase.hpp"
#include "CoinFinite.hpp"

#include "exprVar.hpp"
#include "CouenneJournalist.hpp"
#include "OsiBranchingObject.hpp"

#include "CouenneOrbitObj.hpp"


/// OsiObject for Orbital Branching

class CouenneOrbitObj: public CouenneObject {

public:

  /// empty constructor (for unused objects)
  CouenneOrbitObj ();

  /// Constructor with information for branching point selection strategy
  CouenneOrbitObj (CouenneCutGenerator *cutgen,
		   CouenneProblem *p, 
		   exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Constructor with lesser information, used for infeasibility only
  CouenneOrbitObj (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Destructor
  ~CouenneOrbitObj () {}

  /// Copy constructor
  CouenneOrbitObj (const CouenneOrbitObj &src);

  /// Cloning method
  virtual CouenneObject * clone () const
  {return new CouenneOrbitObj (*this);}

  /// set object parameters by reading from command line
  void setParameters (Bonmin::BabSetupBase *base);

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// compute infeasibility of this variable, |w - f(x)|, where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double checkInfeasibility (const OsiBranchingInformation * info) const;

  /// fix (one of the) arguments of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject *createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

protected:

};

#endif
