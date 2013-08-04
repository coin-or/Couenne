/* $Id$
 *
 * Name:    fillDependence.cpp
 * Author:  Pietro Belotti
 * Purpose: fill in inverse dependence structure, for CouenneObject
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>
#include <set>

#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

/// fill in inverse dependence structure: for each variable x give set
/// of auxiliary variables (or better, their indices) whose expression
/// depends on x

void CouenneProblem::fillDependence (Bonmin::BabSetupBase *base, CouenneCutGenerator *cg) {

  // initialize vector of empty sets
  for (int i=nVars (); i--;)
    dependence_.push_back (std::set <int> ());

  for (std::vector <exprVar *>::iterator i = variables_.begin (); 
       i != variables_.end (); ++i) {

    if (((*i) -> Type () == AUX)                           // consider auxs only
	&& ((*i) -> Image () -> Linearity () > LINEAR)) {  // and nonlinear

      CouenneObject *infeasObj = (*i) -> properObject (cg, this, base, jnlst_);

      if (!(infeasObj -> Reference ())) // found something that will never be infeasible
	continue;

      // add object for this variable
      objects_.push_back (infeasObj);

      std::set <int> deplist;

      // fill the set of independent variables on which the expression
      // associated with *i depends; if empty (should not happen...), skip
      if ((*i) -> Image () -> DepList (deplist, STOP_AT_AUX) == 0)
	continue;

      // build dependence set for this variable
      for (std::set <int>::iterator j = deplist.begin (); j != deplist.end (); ++j) {

	std::set <int> &obj = dependence_ [*j];
	int ind = (*i) -> Index ();
	if (obj.find (ind) == obj.end ())
	  obj.insert (ind);
      }

    } else objects_.push_back (new CouenneObject ()); 
    // null object for original and linear auxiliaries
  }
}
