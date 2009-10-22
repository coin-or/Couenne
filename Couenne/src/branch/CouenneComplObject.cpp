/* $Id$
 *
 * Name:    CouenneComplObject.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Implementation of branching rules for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneComplObject.hpp"
#include "CouenneComplBranchingObject.hpp"


/// Constructor with information for branching point selection strategy
CouenneComplObject::CouenneComplObject (CouenneCutGenerator *c,
					CouenneProblem *p, 
					exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst,
					int sign):
  CouenneObject (c, p, ref, base, jnlst),
  sign_ (sign) {
  jnlst -> Printf (J_DETAILED, J_BRANCHING, 
		   "[created Complementarity constraint object with sign %d]\n", sign);
}


/// Constructor with lesser information, used for infeasibility only
CouenneComplObject::CouenneComplObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst,
					int sign):
  CouenneObject (ref, base, jnlst),
  sign_ (sign) {}


/// Copy constructor
CouenneComplObject::CouenneComplObject (const CouenneComplObject &src): 
  CouenneObject (src),
  sign_ (src.sign_) {}


/// compute infeasibility of this variable, |w - f(x)| (where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  expression **arglist = reference_ -> Image () -> ArgList ();

  int index0 = arglist [0] -> Index (),
      index1 = arglist [1] -> Index ();

  if (sign_) { // it is a xy <= 0 or a xy >= 0 object

    CouNumber 
      x0 = info -> solution_ [index0],
      x1 = info -> solution_ [index1],
      prod = x0*x1;

    if (sign_ < 0) {

      if (prod <= 0) return 0; // object feasible

      way = (x1<=x0); // zero if closer to second orthant (prefer to
                      // branch on "more feasible" variable)

    } else {

      if (prod >= 0) return 0; // object feasible

      way = (x1<=-x0); // zero if closer to second orthant (prefer to
                       // branch on "more feasible" variable)
    }

    // if x1 < x0, it is preferrable to branch with x1=0 instead of x0=0
    // as this is closer to the point
    //way = (x1 < x0) ? 1 : 0;
    return fabs (prod);

  } else { // it is a xy=0 object, use old infeasibility

    CouNumber 
      x0 = fabs (info -> solution_ [index0]),
      x1 = fabs (info -> solution_ [index1]);

    // if x1 < x0, it is preferrable to branch with x1=0 instead of x0=0
    // as this is closer to the point
    way = (x1 < x0) ? 1 : 0;

    return x0 * x1;
  }
}


/// compute infeasibility of this variable, |w - f(x)|, where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::checkInfeasibility (const OsiBranchingInformation * info) const {

  expression **arglist = reference_ -> Image () -> ArgList ();

  int index0 = arglist [0] -> Index (),
      index1 = arglist [1] -> Index ();

  CouNumber
    x0 = info -> solution_ [index0],
    x1 = info -> solution_ [index1],
    prod = x0*x1;

  if (!sign_)
    return fabs (prod);
  else return
	 ((sign_ < 0) && (prod >= 0)) ||
	 ((sign_ > 0) && (prod <= 0)) ? fabs (prod) : 0.;
}


/// create CouenneBranchingObject or CouenneThreeWayBranchObj based
/// on this object
OsiBranchingObject *CouenneComplObject::createBranch (OsiSolverInterface *solver, 
						      const OsiBranchingInformation *info, 
						      int way) const {

  expression **args = reference_ -> Image () -> ArgList ();

  /*  printf ("creating CCobj: %d %d.%d\n", reference_ -> Index (), 
	  args [0] -> Index (),
	  args [1] -> Index ());*/

  return new CouenneComplBranchingObject (solver, this, jnlst_,
					  cutGen_,
					  problem_,
					  args [0],
					  args [1], 
					  way, 0, doFBBT_, doConvCuts_, sign_);
}
