/* $Id$
 *
 * Name:    CouenneBTPerfIndicatorConstr.cpp
 * Author:  Pietro Belotti
 * Purpose: Measures performance of BT in terms of shrunken bounds -- constructors
 *
 * (C) Pietro Belotti, 2011.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneBTPerfIndicator.hpp"
#include "CouenneProblem.hpp"
#include "IpSmartPtr.hpp"


using namespace Couenne;

///
CouenneBTPerfIndicator::CouenneBTPerfIndicator (CouenneProblem *p, const std::string &name):

  name_            (name),
  nFixed_          (0.),
  boundRatio_      (0.),
  shrunkInf_       (0.),
  shrunkDoubleInf_ (0.),
  nProvedInfeas_   (0.),
  weightSum_       (0.),
  oldLB_           (NULL),
  oldUB_           (NULL),
  totalTime_       (0.),
  nRuns_           (0),
  problem_         (p),
  stats_           ((p != NULL) && 
		    (GetRawPtr (p -> Jnlst ()) != NULL) && 
		    (p -> Jnlst () -> ProduceOutput (Ipopt::J_ERROR, J_COUENNE))) {}


///
CouenneBTPerfIndicator::~CouenneBTPerfIndicator () {

  if (totalTime_ > 0. &&
      nRuns_ && 
      problem_)

    if (stats_)
      problem_->Jnlst()->Printf(Ipopt::J_ERROR, J_COUENNE, "Performance of %30s:\t %10gs, %8d runs. fix: %10g shrnk: %10g ubd: %10g 2ubd: %10g infeas: %10g\n", 
	      name_.c_str (),
	      totalTime_, 
	      nRuns_,
	      nFixed_, boundRatio_, shrunkInf_, shrunkDoubleInf_, nProvedInfeas_);

  //weightSum_ * nFixed_, weightSum_ * boundRatio_, weightSum_ * shrunkInf_, weightSum_ * shrunkDoubleInf_, weightSum_ * nProvedInfeas_);

  if (oldLB_) delete [] oldLB_;
  if (oldUB_) delete [] oldUB_;
}


///
CouenneBTPerfIndicator::CouenneBTPerfIndicator (const CouenneBTPerfIndicator &rhs):

  name_            (rhs.name_),
  nFixed_          (rhs.nFixed_),
  boundRatio_      (rhs.boundRatio_),
  shrunkInf_       (rhs.shrunkInf_),
  shrunkDoubleInf_ (rhs.shrunkDoubleInf_),
  nProvedInfeas_   (rhs.nProvedInfeas_),
  weightSum_       (rhs.weightSum_),
  oldLB_           (!rhs.problem_ || rhs.oldLB_ ? NULL : CoinCopyOfArray (rhs.oldLB_, rhs.problem_ -> nVars ())),
  oldUB_           (!rhs.problem_ || rhs.oldUB_ ? NULL : CoinCopyOfArray (rhs.oldUB_, rhs.problem_ -> nVars ())),
  totalTime_       (rhs.totalTime_),
  nRuns_           (rhs.nRuns_),
  problem_         (rhs.problem_),
  stats_           (rhs.stats_) {}


///
CouenneBTPerfIndicator &CouenneBTPerfIndicator::operator= (const CouenneBTPerfIndicator &rhs) {

  name_            = rhs.name_;
  nFixed_          = rhs.nFixed_;
  boundRatio_      = rhs.boundRatio_;
  shrunkInf_       = rhs.shrunkInf_;
  shrunkDoubleInf_ = rhs.shrunkDoubleInf_;
  nProvedInfeas_   = rhs.nProvedInfeas_;
  weightSum_       = rhs.weightSum_;
  oldLB_           = !rhs.problem_ || !rhs.oldLB_ ? NULL : CoinCopyOfArray (rhs.oldLB_, rhs.problem_ -> nVars ());
  oldUB_           = !rhs.problem_ || !rhs.oldUB_ ? NULL : CoinCopyOfArray (rhs.oldUB_, rhs.problem_ -> nVars ());
  totalTime_       = rhs.totalTime_; 
  nRuns_           = rhs.nRuns_;
  problem_         = rhs.problem_;
  stats_           = rhs.stats_;

  return *this;
}


///
void CouenneBTPerfIndicator::setOldBounds (const CouNumber *lb, const CouNumber *ub) const {

  if (problem_) {

    oldLB_ = CoinCopyOfArray (lb, problem_ -> nVars ());
    oldUB_ = CoinCopyOfArray (ub, problem_ -> nVars ());

  } else {

    printf ("CouenneBTPerfIndicator::setOldBounds(): no problem information, exiting\n");
    exit (-1);
  }
}


/// add to timer
void CouenneBTPerfIndicator::addToTimer (double time) const 
{totalTime_ += time;}
