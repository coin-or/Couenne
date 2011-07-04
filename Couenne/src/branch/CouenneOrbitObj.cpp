/* $Id$
 *
 * Name:    CouenneOrbitObj.cpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 * Purpose: Base object for variables (to be used in branching)
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

#include "CouenneProblem.hpp"
#include "CouenneOrbitObj.hpp"
#include "CouenneBranchingObject.hpp"

const CouNumber default_clamp  = 0.2;
const CouNumber max_pseudocost = 1000.;

/// Empty constructor
CouenneOrbitObj::CouenneOrbitObj ():

  CouenneObject () {}


/// Constructor with information for branching point selection strategy
CouenneOrbitObj::CouenneOrbitObj (CouenneCutGenerator *cutgen,
				  CouenneProblem *p, 
				  exprVar *ref, 
				  Bonmin::BabSetupBase *base, 
				  JnlstPtr jnlst):
  CouenneObject (cutgen, p, ref, base, jnlst) {}


/// Constructor with lesser information, used for infeasibility only
CouenneOrbitObj::CouenneOrbitObj (exprVar *ref, 
				  Bonmin::BabSetupBase *base, 
				  JnlstPtr jnlst):

  CouenneObject (ref, base, jnlst) {}


/// Copy constructor
CouenneOrbitObj::CouenneOrbitObj (const CouenneOrbitObj &src):
  CouenneObject       (src) {}


/// apply the branching rule
OsiBranchingObject *CouenneOrbitObj::createBranch (OsiSolverInterface *si,
						   const OsiBranchingInformation *info,
						   int way) const {

  return NULL;
}


// set point at current LP solution
double CouenneOrbitObj::feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const {
  return 0;
}



/// non-linear infeasibility -- not called by independent's CouenneVarObject
double CouenneOrbitObj::infeasibility (const OsiBranchingInformation *info, int &way) const {

  return 0;
}


/// non-linear infeasibility -- no need for the domain's push
/// instruction as this is called from
/// CouenneVarObject::infeasibility()
double CouenneOrbitObj::checkInfeasibility (const OsiBranchingInformation *info) const {

  return 0;

}
