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

#include "CouenneExprGroup.hpp"

using namespace Couenne;

/// Empty constructor
CouenneOrbitObj::CouenneOrbitObj ():

  CouenneObject () {}


/// Constructor with information for branching point selection strategy
CouenneOrbitObj::CouenneOrbitObj (CouenneCutGenerator *cutgen,
				  CouenneProblem *p, 
				  exprVar *ref, 
				  Bonmin::BabSetupBase *base, 
				  JnlstPtr jnlst):
  CouenneObject (cutgen, p, ref, base, jnlst) {

  // nautyGraph_ = createNautyGraph (p);

  // create graph

  // create p -> nVars () nodes for (aux+orig) variables.


  // join them

  for (std::vector <exprVar *>:: iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end (); ++i) {

    if ((*i) -> Type () == AUX) {

      // this is an auxiliary variable

      // add node in nauty graph for its index, (*i) -> Index ()

      if        ((*i) -> Image () -> Type () == N_ARY) {

	for (int a=0; a < (*i) -> Image () -> nArgs(); a++) {

	  expression *arg = (*i) -> Image () -> ArgList () [a];

	  if (arg -> Type () == AUX) {

	    // a-th argument is an auxiliary variable

	  } else if (arg -> Type () == VAR) {

	    CouNumber lb, ub;

	    arg -> getBounds (lb, ub);

	    // a-th argument is an original variable

	  } else {

	    assert (arg -> Type () == CONST);

	    // this is a constant.

	  }
	}


	if ((*i) -> Image () -> code () == COU_EXPRGROUP) {

	  // dynamic_cast it to an exprGroup

	  exprGroup *e = dynamic_cast <exprGroup *> ((*i) -> Image ());

	  // add a node for e -> getC0 ();

	  // for each term add nodes for their non-one coefficients and their variable

	  for (exprGroup::lincoeff::iterator el = e ->lcoeff().begin (); el != e -> lcoeff ().end (); ++el) {

	    // coefficient = el -> second

	    // variable index is el -> first -> Index ()
	  }
	  
	}

      } else if ((*i) -> Image () -> Type () == UNARY) {

      }

    } else {

      // this is an original variable

    }
  }

  // create partitions for log only

  for (std::vector <exprVar *>:: iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end (); ++i) {

    if ((*i) -> Type () == AUX) {

      // this is an auxiliary variable

      if ((*i) -> Image () -> code () == COU_EXPRLOG) {
	printf ("log is %d\n", (*i) -> Index ());
      }

    } else {

      // this is an original variable

    }
  }
}


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
