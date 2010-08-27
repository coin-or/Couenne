/* $Id$
 *
 * Name:    TwoImpliedConstructors.cpp
 * Author:  Pietro Belotti
 * Purpose: Bound Tightening on pairs of linear constraints
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#include "BonRegisteredOptions.hpp"
#include "CouenneTwoImplied.hpp"

#include "CglCutGenerator.hpp"
#include "CouenneJournalist.hpp"

#include "IpOptionsList.hpp"

using namespace Couenne;

/// constructor
CouenneTwoImplied::CouenneTwoImplied (JnlstPtr jnlst,
				      const Ipopt::SmartPtr <Ipopt::OptionsList> options):
  jnlst_ (jnlst) {

  options -> GetIntegerValue ("two_implied_max_trials", nMaxTrials_, "couenne.");
}


/// copy constructor
CouenneTwoImplied::CouenneTwoImplied (const CouenneTwoImplied &src):

  CglCutGenerator (src),
  jnlst_          (src.jnlst_),
  nMaxTrials_     (src.nMaxTrials_) {}


/// destructor
CouenneTwoImplied::~CouenneTwoImplied () {}


/// Add list of options to be read from file
void CouenneTwoImplied::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("two_implied_bt",
     "The frequency (in terms of nodes) at which Couenne two-implied bounds are tightened.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. \
Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. \
A negative number -n means that generation should be attempted at the root node, and if \
successful it can be repeated at every n nodes, otherwise it is stopped altogether."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("two_implied_max_trials",
     "The number of iteration at each call to the cut generator.",
     1, 3,
     "");
}
