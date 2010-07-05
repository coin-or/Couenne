/* $Id$
 *
 * Name:    CouenneFixPoint.cpp
 * Author:  Pietro Belotti
 * Purpose: fixpoint bound tightener -- constructors
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneFixPoint.hpp"
#include <string>

using namespace Couenne;

/// constructor
CouenneFixPoint::CouenneFixPoint (CouenneProblem *p,
				  const Ipopt::SmartPtr<Ipopt::OptionsList> options):
  problem_ (p) {

  std::string s;
  options -> GetStringValue ("fixpoint_bt_model", s, "couenne."); 
  extendedModel_ = (s == "extended");
}


/// copy constructor
CouenneFixPoint::CouenneFixPoint  (const CouenneFixPoint &rhs):
  extendedModel_ (rhs.extendedModel_),
  problem_       (rhs.problem_) {}


/// destructor
CouenneFixPoint::~CouenneFixPoint () {}


/// Add list of options to be read from file
void CouenneFixPoint::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("fixpoint_bt",
     "The frequency (in terms of nodes) at which Fix Point Bound Tightening is performed.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. "
     "Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. "
     "A negative number -n means that generation should be attempted at the root node, and if successful it can be repeated at every n nodes, otherwise it is stopped altogether."
     );

  roptions -> AddStringOption2
    ("fixpoint_bt_model",
     "Choose whether to add an extended fixpoint LP model or a more compact one.",
     "extended",
     "extended", "",
     "compact", "",
     "");
}
