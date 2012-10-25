/* $Id$
 *
 * Name:    CouenneSdpCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: wrapper for Couenne to insert sdpcuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneSdpCuts.hpp"
#include "CouenneProblem.hpp"

#include "IpOptionsList.hpp"

using namespace Couenne;

/// Constructor
CouenneSdpCuts::CouenneSdpCuts (CouenneProblem *p) {

  // 0) Search for X \succeq 0 constraints, add matrix to minors for
  //    each such constraint

  // 1) Construct matrix with entries x_i, x_j
  // 2) Block-partition it (optional), obtain matrices
}

/// Destructor
CouenneSdpCuts::~CouenneSdpCuts () {

  // Destroy matrix structures
}

/// Copy constructor
CouenneSdpCuts::CouenneSdpCuts (const CouenneSdpCuts &rhs) {

  // Copy matrices
}

/// Assignment
CouenneSdpCuts &CouenneSdpCuts::operator= (const CouenneSdpCuts &rhs) {

  return *this;
  //  CouenneSdpCouenneSdpCuts (rhs);
}

/// Assignment
CglCutGenerator *CouenneSdpCuts::clone () const
{return new CouenneSdpCuts (*this);}


/// The main CglCutGenerator
void CouenneSdpCuts::generateCuts (const OsiSolverInterface &si, 
				   OsiCuts &cs, 
				   const CglTreeInfo treeinfo) const {

  // for each matrix, call CutGen from Andrea's code
}

/// Add list of options to be read from file
void CouenneSdpCuts::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("sdp_cuts",
     "The frequency (in terms of nodes) at which Couenne sdp cuts are generated.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. "
     "Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. "
     "A negative number -n means that generation should be attempted at the root node, and if successful it can be repeated at every n nodes, otherwise it is stopped altogether."
    );
}
