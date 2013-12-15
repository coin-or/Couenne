/* $Id$
 *
 * Name:    FixPointConstructors.cpp
 * Author:  Pietro Belotti
 * Purpose: fixpoint bound tightener -- constructors
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneFixPoint.hpp"
#include <string>

using namespace Couenne;

/// constructor
CouenneFixPoint::CouenneFixPoint (CouenneProblem *p,
				  const Ipopt::SmartPtr<Ipopt::OptionsList> options):
  problem_    (p),
  firstCall_  (true),
  CPUtime_    (0.),
  nTightened_ (0),
  perfIndicator_ (p, "Fixed Point LP") {

  std::string s;
  options -> GetStringValue ("fixpoint_bt_model", s, "couenne."); 
  extendedModel_ = (s == "extended");
}


/// copy constructor
CouenneFixPoint::CouenneFixPoint (const CouenneFixPoint &rhs):
  extendedModel_ (rhs.extendedModel_),
  problem_       (rhs.problem_),
  firstCall_     (rhs.firstCall_),
  CPUtime_       (rhs.CPUtime_),
  nTightened_    (rhs.nTightened_),
  perfIndicator_ (rhs.perfIndicator_) {}


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
     "compact",
     "extended", "Extended model with variables for lower/upper bounds of right-hand sides (see paper by Belotti, Cafieri, Lee, Liberti)",
     "compact", "Compact equivalent model obtained by projecting lower/upper bounds of rhs",
     "The \"extended\" option is for debugging purposes; the compact formulation is equivalent and this option should be used");
}
