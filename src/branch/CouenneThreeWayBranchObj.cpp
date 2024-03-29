/*
 *
 * Name:    CouenneThreeWayBranchObj.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneThreeWayBranchObj.hpp"

using namespace Ipopt;
using namespace Couenne;

/// Constructor. Get a variable as an argument and set value_ through
/// a call to operator () of that exprAux.
CouenneThreeWayBranchObj::CouenneThreeWayBranchObj (JnlstPtr jnlst,
						    expression *brVar,
						    CouNumber lcrop,
						    CouNumber rcrop,
						    int way
						    //bool isint
						    ):
  brVar_   (brVar),
  lcrop_   (lcrop),
  rcrop_   (rcrop),
  jnlst_   (jnlst) {

  numberBranches_ = 3;

  // if none of these three, set to 0 and do code below
  firstBranch_ = (way == THREE_LEFT)   ? 0 :
                 (way == THREE_CENTER) ? 1 :
                 (way == THREE_RIGHT)  ? 2 : 0;

  if (way == THREE_RAND) { // pick first branch randomly
    CouNumber rnd = 3. * CoinDrand48 ();
    firstBranch_ = (rnd < 1.) ? 0 : (rnd < 2.) ? 1: 2;
  }
}


/** \brief Execute the actions required to branch, as specified by the
	   current state of the branching object, and advance the
	   object's state.
	   Returns change in guessed objective on next branch
*/

double CouenneThreeWayBranchObj::branch (OsiSolverInterface * solver) {

  //       -1 if "<= a"  node
  // way =  0 if "[a,b]" node
  //        1 if ">= b"  node

  int way = 0;

  switch (branchIndex_) {
    // if first offspring, let firstBranch_ decide who's first
  case 0: way = firstBranch_;                break;
  case 1: way = (firstBranch_ == 0) ? 1 : 0; break;
  case 2: way = (firstBranch_ == 2) ? 1 : 2; break;
  default: jnlst_->Printf(J_WARNING, J_BRANCHING,
			  "Warning, branchIndex_ has a strange value (%d)\n", branchIndex_);
  }

  // set lower or upper bound (round if this variable is integer)

  int  index   = brVar_ -> Index ();
  bool integer = brVar_ -> isInteger ();

  CouNumber
    l = solver -> getColLower () [index],
    u = solver -> getColUpper () [index];

  if (lcrop_ < l) lcrop_ = l;
  if (rcrop_ > u) rcrop_ = u;

  switch (--way) { // from {0,1,2} to {-1,0,1}

  case -1: solver -> setColUpper (index, integer ? floor (lcrop_) : lcrop_); break; // left
  case  0: solver -> setColLower (index, integer ? ceil  (lcrop_) : lcrop_);
           solver -> setColUpper (index, integer ? floor (rcrop_) : rcrop_); break; // central
  case  1: solver -> setColLower (index, integer ? ceil  (rcrop_) : rcrop_); break; // right
  default: jnlst_->Printf(J_WARNING, J_BRANCHING, "Couenne: branching on nonsense way %d\n", way);
  }

  // TODO: apply bound tightening

  if (jnlst_->ProduceOutput(J_DETAILED, J_BRANCHING)) {
    switch (way) {
    case -1: jnlst_->Printf(J_DETAILED, J_BRANCHING,
			    "#3# Branch: x%d <= %g\n",               index, lcrop_); break; // left
    case  0: jnlst_->Printf(J_DETAILED, J_BRANCHING,
			    "#3# Branch: %g <= x%d <= %g\n", lcrop_, index, rcrop_); break; // center
    case  1: jnlst_->Printf(J_DETAILED, J_BRANCHING,
			    "#3# Branch: x%d >= %g\n",               index, rcrop_); break; // right
    default: jnlst_->Printf(J_DETAILED, J_BRANCHING, "Couenne: branching on nonsense way %d\n", way);
    }
  }

  branchIndex_++;

  return 0.; // estimated change in objective function
}
