/* $Id$
 *
 * Name:    CouenneAggrProbing.cpp
 * Author:  Giacomo Nannicini
 * Purpose: Aggressive probing
 *
 * (C) Giacomo Nannicini, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneAggrProbing.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprOpp.hpp"
//#include "BonCbc.hpp"
#include "CouenneBab.hpp"
#include "CouenneCutGenerator.hpp"
#include <string>

#define COUENNE_AGGR_PROBING_FINITE_BOUND 1.0e+10
#define COUENNE_AGGR_PROBING_MIN_INTERVAL 1.0e-2
#define COUENNE_AGGR_PROBING_BND_RELAX COUENNE_EPS

using namespace Couenne;

CouenneAggrProbing::CouenneAggrProbing(CouenneSetup *setup,
				       const Ipopt::SmartPtr<Ipopt::OptionsList> options)
{
  couenne_ = setup;
  numCols_ = couenne_->couennePtr()->Problem()->nVars();
  maxTime_ = COIN_DBL_MAX;
  maxFailedSteps_ = 10;
  maxNodes_ = 1000;
  initCutoff_ = COUENNE_INFINITY;
  restoreCutoff_ = false;

}

CouenneAggrProbing::CouenneAggrProbing(const CouenneAggrProbing &rhs){
  couenne_ = new CouenneSetup(*rhs.couenne_);
  numCols_ = rhs.numCols_;
  maxTime_ = rhs.maxTime_;
  maxFailedSteps_ = rhs.maxFailedSteps_;
  maxNodes_ = rhs.maxNodes_;
  initCutoff_ = rhs.initCutoff_;
  restoreCutoff_ = rhs.restoreCutoff_;
}

CouenneAggrProbing::~CouenneAggrProbing(){
}

void CouenneAggrProbing::registerOptions(Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {
  // Nothing for the moment, but will be added later as needed
}

double CouenneAggrProbing::getMaxTime() const {
  return maxTime_;
}

void CouenneAggrProbing::setMaxTime(double value){
  maxTime_ = value;
}

int CouenneAggrProbing::getMaxFailedSteps() const {
  return maxFailedSteps_;
}

void CouenneAggrProbing::setMaxFailedSteps(int value){
  maxFailedSteps_ = value;
}

int CouenneAggrProbing::getMaxNodes() const {
  return maxNodes_;
}

void CouenneAggrProbing::setMaxNodes(int value){
  maxNodes_ = value;
}

bool CouenneAggrProbing::getRestoreCutoff() const {
  return restoreCutoff_;
}

void CouenneAggrProbing::setRestoreCutoff(bool value){
  restoreCutoff_ = value;
}

double CouenneAggrProbing::probeVariable(int index, bool probeLower){

  // Useful objects for easy access
  OsiSolverInterface* nlp = couenne_->nonlinearSolver();
  OsiSolverInterface* lp = couenne_->continuousSolver();
  CouenneProblem* problem = couenne_->couennePtr()->Problem();

  // Save initial bounds
  double initUpper = lp->getColUpper()[index];
  double initLower = lp->getColLower()[index];

  double* initLowerLp = new double[numCols_];
  double* initUpperLp = new double[numCols_];

  memcpy(initLowerLp, lp->getColLower(), numCols_*sizeof(double));
  memcpy(initUpperLp, lp->getColUpper(), numCols_*sizeof(double));

  if (initUpper < initLower + COUENNE_EPS){
    // Variable is fixed, so we can't tighten
    return ((probeLower) ? initLower : initUpper);
  }

  // Index of the aux variable representing the objective function
  int indobj = problem->Obj(0)->Body()->Index();

  // Initial cutoff value
  double initCutoff = problem->Ub()[indobj];

  double* initCutoffSol = NULL; 

  if (restoreCutoff_ && problem->getCutOff() < COUENNE_INFINITY){
    initCutoffSol = new double[numCols_];
    memcpy(initCutoffSol, problem->getCutOffSol(), numCols_*sizeof(double));
  }

  // Save parameters
  Bonmin::BabSetupBase::NodeComparison initNodeComparison = 
    couenne_->nodeComparisonMethod();
  int initMaxNodes = couenne_->getIntParameter(Bonmin::BabSetupBase::MaxNodes);
  double initMaxTime = couenne_->getDoubleParameter(Bonmin::BabSetupBase::MaxTime);
  int initMaxSol = couenne_->getIntParameter(Bonmin::BabSetupBase::MaxSolutions);
  couenne_->setNodeComparisonMethod(Bonmin::BabSetupBase::bestBound);
  //couenne_->nodeComparisonMethod() = Bonmin::BabSetupBase::bestBound;
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxNodes, maxNodes_);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxSolutions, COIN_INT_MAX);
  problem->setCheckAuxBounds(true);

  /// First store, then disable all heuristics.
  Bonmin::BabSetupBase::HeuristicMethods heuristics = couenne_->heuristics();
  couenne_->heuristics().clear();

  double currentBound = (probeLower) ? initLower : initUpper;
  double startTime = CoinCpuTime();
  int failedSteps = 0;
  double intervalSize = 0.0;
  double tryBound = 0.0;

  int iter = 0;

  if (probeLower)
    std::cout << "Probing lower on var " << index << std::endl;
  else
    std::cout << "Probing upper on var " << index << std::endl;

  if ((fabs(currentBound) > COUENNE_AGGR_PROBING_FINITE_BOUND) &&
      ((probeLower && initUpper > -COUENNE_AGGR_PROBING_FINITE_BOUND) ||
       (!probeLower && initLower < COUENNE_AGGR_PROBING_FINITE_BOUND))){
    // The bound is too large to apply the standard probing method;
    // try to reduce it to a finite value. We only do this if we want
    // to probe a variable on an infinite (or close to) bound, and the
    // other bound of the variable is sufficiently far away
    if (probeLower){
      tryBound = -COUENNE_AGGR_PROBING_FINITE_BOUND;
      lp->setColLower(index, currentBound);
      problem->Lb()[index] = currentBound;
      lp->setColUpper(index, tryBound);
      problem->Ub()[index] = tryBound;
      if (index < problem->nOrigVars()){
	nlp->setColLower(index, currentBound);
	nlp->setColUpper(index, tryBound);
      }
    }
    else{
      tryBound = COUENNE_AGGR_PROBING_FINITE_BOUND;
      lp->setColLower(index, tryBound);
      problem->Lb()[index] = tryBound;
      lp->setColUpper(index, currentBound);
      problem->Ub()[index] = currentBound;
      if (index < problem->nOrigVars()){
	nlp->setColLower(index, tryBound);
	nlp->setColUpper(index, currentBound);
      }
    }

    /// Setup Branch-and-Bound limits
    couenne_->setDoubleParameter(Bonmin::BabSetupBase::MaxTime, 
				 CoinMin(maxTime_-(CoinCpuTime()-startTime),
					 maxTime_*0.5));

    if (restoreCutoff_){
      problem->resetCutOff(initCutoff);
      problem->Ub()[indobj] = initCutoff;
      problem->installCutOff();
    }

    std::cout << "Iteration " << iter << ", current bound " << currentBound
	      << ", try bound " << tryBound << std::endl;

    /// Now do Branch-and-Bound and see if probing succeeded
    // Bonmin::Bab bb;
    // bb.setUsingCouenne(true);

    CouenneBab bb;

    bb(couenne_);
    if (bb.model().isProvenInfeasible()){
      /// Problem is infeasible; therefore, probing was successful.
      currentBound = tryBound;
      std::cout << "Probing succeeded; brought to finite" << std::endl;
    }
    else{
      /// Problem is not infeasible; we failed
      std::cout << "Probing failed; still infinity, exit" << std::endl;
    }    
    iter++;
  }

  // Now that we have a finite bound, pick size of the probing interval

  // Override (for testing - will be chosen automatically in final
  // implementation)
  intervalSize = 0.1;

  if (intervalSize < COUENNE_AGGR_PROBING_MIN_INTERVAL){
    intervalSize = COUENNE_AGGR_PROBING_MIN_INTERVAL;
  }
  
  while ((fabs(currentBound) <= COUENNE_AGGR_PROBING_FINITE_BOUND) &&
	 ((CoinCpuTime() - startTime) < maxTime_) &&
	 (failedSteps < maxFailedSteps_) &&
	 (intervalSize >= COUENNE_AGGR_PROBING_MIN_INTERVAL) && 
	 iter < 100){

    /// Set the bound that we want to try
    if (probeLower){
      tryBound = currentBound + intervalSize;
      if (tryBound > initUpper){
	// It does not make sense to use bounds larger than the initial
	// ones
	tryBound = initUpper;
      }
      if (lp->isInteger(index)){
	tryBound = floor(tryBound);
      }
      // Relax bounds a little bit
      lp->setColLower(index, currentBound - COUENNE_AGGR_PROBING_BND_RELAX);
      problem->Lb()[index] = currentBound - COUENNE_AGGR_PROBING_BND_RELAX;
      lp->setColUpper(index, tryBound + COUENNE_AGGR_PROBING_BND_RELAX);
      problem->Ub()[index] = tryBound + COUENNE_AGGR_PROBING_BND_RELAX;
      if (index < problem->nOrigVars()){
	nlp->setColLower(index, currentBound - COUENNE_AGGR_PROBING_BND_RELAX);
	nlp->setColUpper(index, tryBound + COUENNE_AGGR_PROBING_BND_RELAX);
      }
    }
    else{
      tryBound = currentBound - intervalSize;
      if (tryBound < initLower){
	// It does not make sense to use bounds larger than the initial
	// ones
	tryBound = initLower;
      }
      if (lp->isInteger(index)){
	tryBound = ceil(tryBound);
      }
      // Relax bounds a little bit
      lp->setColLower(index, tryBound - COUENNE_AGGR_PROBING_BND_RELAX);
      problem->Lb()[index] = tryBound - COUENNE_AGGR_PROBING_BND_RELAX;
      lp->setColUpper(index, currentBound + COUENNE_AGGR_PROBING_BND_RELAX);
      problem->Ub()[index] = currentBound + COUENNE_AGGR_PROBING_BND_RELAX;
      if (index < problem->nOrigVars()){
	nlp->setColLower(index, tryBound - COUENNE_AGGR_PROBING_BND_RELAX);
	nlp->setColUpper(index, currentBound + COUENNE_AGGR_PROBING_BND_RELAX);
      }
    }

    lp->resolve();
    problem->domain()->push(numCols_, lp->getColSolution(),
			    lp->getColLower(), lp->getColUpper());

    /// Setup Branch-and-Bound limits
    couenne_->setDoubleParameter(Bonmin::BabSetupBase::MaxTime, 
    				 CoinMin(maxTime_-(CoinCpuTime()-startTime),
    					 maxTime_*0.5));

    if (restoreCutoff_){
      problem->Ub()[indobj] = initCutoff;
      problem->resetCutOff(initCutoff);
      problem->installCutOff();
    }

    std::cout << "Iteration " << iter << ", current bound " << currentBound
	      << ", try bound " << tryBound << std::endl;

    /// Now do Branch-and-Bound and see if probing succeeded
    // Bonmin::Bab bb;
    // bb.setUsingCouenne(true);

    CouenneBab bb;
    bb(couenne_);

    problem->domain()->pop();

    double obj = 0.0;
    /// Is the search in the current interval complete?
    bool intervalSearched = (bb.model().isProvenOptimal() || 
			     bb.model().isProvenInfeasible());

    if ((!intervalSearched) || // If the search is not complete
	(restoreCutoff_ && // or we don't want to accept new solutions
	 problem->getCutOffSol() && // and we have a new feasible solution
	 problem->checkNLP(problem->getCutOffSol(), obj, true))){
      /// Try again in a smaller interval
      if (lp->isInteger(index) && fabs(tryBound-currentBound) < 0.5){
	/// There is no smaller interval that we can try; bail out
	failedSteps = maxFailedSteps_;
      }
      else{
	intervalSize /= 2;
      }
      failedSteps++;
      std::cout << "Probing failed; shrinking interval" << std::endl;
    }
    else{
      /// We fully explored the current interval, and there is no
      /// feasible solution, or there is a solution and we have
      /// already updated the cutoff. So, we can eliminate the current
      /// interval. We also double the size of the search interval.
      if (lp->isInteger(index) && fabs(tryBound-currentBound) < 0.5){
	/// Make sure we increase by at least one if it is an integer
	/// variable
	intervalSize = 1.0;
      }
      else{
	intervalSize *= 2;
      }
      currentBound = tryBound;
      if (lp->isInteger(index)){
	if (probeLower){
	  currentBound += 1.0;
	}
	else {
	  currentBound -= 1.0;
	}
      }
      failedSteps = 0;
      std::cout << "Probing succeeded; enlarging interval" << std::endl;
    }

    // Check early termination condition: if we manage to fix the
    // variable (unlikely), there is nothing more we can do
    if ((probeLower && fabs(currentBound-initUpper) < COUENNE_EPS) ||
	(!probeLower && fabs(currentBound-initLower) < COUENNE_EPS)){
      failedSteps = maxFailedSteps_;
    }

    // Reset cutoff
    if (restoreCutoff_){
      problem->Ub()[indobj] = initCutoff;
      problem->resetCutOff(initCutoff);
      problem->installCutOff();
    }

    problem->domain()->pop();

    iter++;
  }

  /// Restore initial bounds (we do not want to modify the
  /// CouenneSetup object: the caller will do that, if needed)
  lp->setColLower(initLowerLp);
  lp->setColUpper(initUpperLp);
  nlp->setColLower(initLowerLp);
  nlp->setColUpper(initUpperLp);
  memcpy(problem->Lb(), initLowerLp, numCols_*sizeof(double));
  memcpy(problem->Ub(), initUpperLp, numCols_*sizeof(double));

  /// Restore parameters and heuristics
  problem->setCheckAuxBounds(false);
  //couenne_->nodeComparisonMethod() = initNodeComparison;
  couenne_->setNodeComparisonMethod(initNodeComparison);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxSolutions, initMaxSol);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxNodes, initMaxNodes);
  couenne_->setDoubleParameter(Bonmin::BabSetupBase::MaxTime, initMaxTime);
  couenne_->heuristics() = heuristics;

  /// Restore cutoff
  if (restoreCutoff_){
    problem->resetCutOff();
    problem->setCutOff(initCutoff, initCutoffSol);
    if (initCutoffSol){
      delete[] initCutoffSol;
    }
  }
  
  delete[] initLowerLp;
  delete[] initUpperLp;

  /// We are done; return best bound found.
  return currentBound;

}

double CouenneAggrProbing::probeVariable2(int index, bool probeLower){
  // Does not work! It's impossible to get Maximization problems working...
  // Adding extra variables doesn't seem to do the trick

  // Useful objects for easy access
  OsiSolverInterface* lp = couenne_->continuousSolver();
  CouenneProblem* problem = couenne_->couennePtr()->Problem();

  // Save initial bounds
  double initUpper = lp->getColUpper()[index];
  double initLower = lp->getColLower()[index];

  if (initUpper < initLower + COUENNE_EPS){
    // Variable is fixed, so we can't probe anything
    return ((probeLower) ? initLower : initUpper);
  }

  /// Modify the CouenneSetup object to use our options.
  /// We store the initial values of all parameters that we modify,
  /// so that we can restore them when we are done.
  Bonmin::BabSetupBase::NodeComparison initNodeComparison = 
    couenne_->nodeComparisonMethod();
  int initMaxNodes = couenne_->getIntParameter(Bonmin::BabSetupBase::MaxNodes);
  double initMaxTime = couenne_->getDoubleParameter(Bonmin::BabSetupBase::MaxTime);
  int initMaxSol = couenne_->getIntParameter(Bonmin::BabSetupBase::MaxSolutions);
  couenne_->setNodeComparisonMethod (Bonmin::BabSetupBase::bestBound);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxNodes, maxNodes_);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxSolutions, COIN_INT_MAX);

  /// First store, then disable all heuristics - we do not need upper
  /// bounds, plus they probably use a NLP object which we do not know
  /// how to modify
  Bonmin::BabSetupBase::HeuristicMethods heuristics = couenne_->heuristics();
  couenne_->heuristics().clear();

  /// Now, store and modify objective function in the CouenneProblem object
  double* initLpObj = new double[numCols_];
  memcpy(initLpObj, lp->getObjCoefficients(), numCols_*sizeof(double));
  expression* initProbObj = problem->Obj(0)->Body();

  double* newLpObj = new double[numCols_];
  memset(newLpObj, 0, numCols_*sizeof(double));

  //expression* exprObj = NULL; // PB: unused
  expression* extraVar = NULL;

  lp->writeLp("before");

  if (probeLower){
    std::cout << "Probing LOWER" << std::endl;
    // Set LP objective function
    newLpObj[index] = 1.0;
    lp->setObjective(newLpObj);

    lp->writeLp("lower");

    // Set CouenneProblem objective
    problem->Obj(0)->Body(problem->Variables()[index]);
    // couenne_->setDoubleParameter(Bonmin::BabSetupBase::Cutoff, COIN_DBL_MAX);
    // problem->setCutOff(initUpper, lp->getColUpper());
    // problem->installCutOff();

  }
  else{
    // We cannot maximize an objective function in Couenne; so, we
    // have to introduce an additional variable, equal to the opposite
    // of the variable that we want to maximize, and minimize that.

    // Add one column and one row to the LP
    int extraCol = numCols_;
    lp->setObjective(newLpObj);
    lp->addCol(0, NULL, NULL, -initUpper, -initLower, 1.0);

    // Create the row x_{extraCol} = -x_{index}
    int rowIndices[2] = {index, extraCol};
    double rowElements[2] = {1.0, 1.0};
    lp->addRow(2, rowIndices, rowElements, 0.0, 0.0);
    lp->resolve();

    // Create the expression x_{extraCol} = -x_{index} in CouenneProblem
    extraVar = problem->addVariable(lp->isInteger(index), NULL);
    // exprObj = new exprOpp(problem->Variables()[index]->clone());
    // problem->addEQConstraint(extraVar, exprObj);
    problem->Obj(0)->Body(extraVar);

    // couenne_->setDoubleParameter(Bonmin::BabSetupBase::Cutoff, COIN_DBL_MAX);
    // problem->setCutOff(-initLower, lp->getColLower());
    // problem->installCutOff();

    lp->writeLp("upper");
  }

  couenne_->setNodeComparisonMethod (Bonmin::BabSetupBase::bestBound);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxNodes, maxNodes_);
  couenne_->setDoubleParameter(Bonmin::BabSetupBase::MaxTime, 
			       maxTime_);

  // Bonmin::Bab bb;
  // bb.setUsingCouenne(true);

  CouenneBab bb;
  bb(couenne_);

  double bestBound = bb.model().getBestPossibleObjValue();

  std::cout << "Obtained bound: " << bb.model().getBestPossibleObjValue() << std::endl;


  /// Restore parameters
  couenne_->setNodeComparisonMethod (initNodeComparison);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxNodes, initMaxNodes);
  couenne_->setDoubleParameter(Bonmin::BabSetupBase::MaxTime, initMaxTime);
  couenne_->setIntParameter(Bonmin::BabSetupBase::MaxSolutions, initMaxSol);
  couenne_->heuristics() = heuristics;

  if (!probeLower){
    int extra = lp->getNumCols()-1;
    lp->deleteCols(1, &extra);
    extra = lp->getNumRows()-1;
    lp->deleteRows(1, &extra);
    problem->Variables().pop_back();
    delete extraVar;
    // delete exprObj;
  }

  lp->setObjective(initLpObj);
  problem->Obj(0)->Body(initProbObj);

  delete[] initLpObj;
  delete[] newLpObj;

  return ((probeLower) ? bestBound : -bestBound);

}
