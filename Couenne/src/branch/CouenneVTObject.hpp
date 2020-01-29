/* $Id$
 *
 * Name:    CouenneVTObject.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for branching on variables using violation transfer
 *
 * (C) Carnegie-Mellon University, 2008-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEVTOBJECT_HPP
#define COUENNEVTOBJECT_HPP

#include "CouenneVarObject.hpp"

namespace Couenne {

/// OsiObject for violation transfer on variables in a MINLP
class COUENNELIB_EXPORT CouenneVTObject: public CouenneVarObject {

public:

  /// Constructor with information for branching point selection strategy
  CouenneVTObject (CouenneCutGenerator *c,
		   CouenneProblem *p,
		   exprVar *ref, 
		   Bonmin::BabSetupBase *base, 
		   JnlstPtr jnlst,
		   int varSelection // either OSI_SIMPLE or OSI_STRONG
		   ):

    CouenneVarObject (c, p, ref, base, jnlst, varSelection) {}

  /// Copy constructor
  CouenneVTObject (const CouenneVTObject &src):
    CouenneVarObject (src) {}

  /// Destructor
  ~CouenneVTObject () {}

  /// Cloning method
  virtual CouenneObject *clone () const
  {return new CouenneVTObject (*this);}

  /// compute infeasibility of this variable x as the sum/min/max of
  /// all infeasibilities of auxiliaries w whose defining function
  /// depends on x |w - f(x)|
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;
};

}

#endif
