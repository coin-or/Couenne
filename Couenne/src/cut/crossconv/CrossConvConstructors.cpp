/* $Id$
 *
 * Name:    CrossConvConstructors.cpp
 * Author:  Pietro Belotti
 * Purpose: Convexification cuts on redundant relationships between auxiliaries
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "BonRegisteredOptions.hpp"
#include "CouenneCrossConv.hpp"

#include "CglCutGenerator.hpp"
#include "CouenneJournalist.hpp"
#include "IpOptionsList.hpp"

using namespace Couenne;

/// constructor
CouenneCrossConv::CouenneCrossConv (CouenneProblem *p,
				    JnlstPtr,
				    const Ipopt::SmartPtr <Ipopt::OptionsList>) {
  setup ();
}

/// copy constructor
CouenneCrossConv::CouenneCrossConv  (const CouenneCrossConv &src):

  CglCutGenerator (src) {
}

/// destructor
CouenneCrossConv::~CouenneCrossConv () {


}

/// Add list of options to be read from file
void CouenneCrossConv::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("crossconv_cuts",
     "The frequency (in terms of nodes) at which Couenne cross-aux convexification cuts are generated.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. "
     "Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. "
     "A negative number -n means that generation should be attempted at the root node, and if successful it can be repeated at every n nodes, otherwise it is stopped altogether."
    );
}
