/* $Id$
 *
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 *          Pietro Belotti, Lehigh University
 *          Francois Margot, Carnegie Mellon University
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "BonChooseVariable.hpp"

#include "CouenneObject.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneRecordBestSol.hpp"

// The recommended ones:
#define FM_SORT_STRONG
#define FM_SEC_SORT_USEFUL
#define USE_NOT_TRUSTED

//#define TRACE_STRONG
//#define TRACE_STRONG2
//#define FM_ALWAYS_SORT
//#define USE_SMALL_GAP
//#define OLD_STYLE

using namespace Ipopt;
using namespace Couenne;

const CouNumber estProdEps = 1e-6;

#ifdef COIN_HAS_NTY
#include "Nauty.h"
#endif

  /// constructor
  CouenneChooseStrong::CouenneChooseStrong (Bonmin::BabSetupBase &b, CouenneProblem* p, JnlstPtr jnlst) :

  Bonmin::BonChooseVariable (b, b.continuousSolver()),
  problem_          (p),
  jnlst_            (jnlst),
  branchtime_       (0.) {

    std::string s;

    b.options () -> GetStringValue ("pseudocost_mult_lp", s, "couenne.");
    pseudoUpdateLP_ = (s == "yes");      

    b.options () -> GetStringValue ("estimate_select", s, "couenne.");
    estimateProduct_ = (s == "product");

    b.options () -> GetStringValue ("trust_strong", s, "couenne.");

    // trust solution from strong branching to provide right LB

    setTrustStrongForSolution (s == "yes");
    setTrustStrongForBound    (s == "yes");
  }

  /// copy constructor
  CouenneChooseStrong::CouenneChooseStrong (const CouenneChooseStrong& rhs) :
    Bonmin::BonChooseVariable (rhs),
    problem_          (rhs.problem_),
    pseudoUpdateLP_   (rhs.pseudoUpdateLP_),
    estimateProduct_  (rhs.estimateProduct_),
    jnlst_            (rhs.jnlst_),
    branchtime_       (rhs.branchtime_)
  {}

  /// destructor
  CouenneChooseStrong::~CouenneChooseStrong()
  {if (branchtime_ > 1e-9) jnlst_ -> Printf (J_ERROR, J_BRANCHING, "Strong branching: total time %g\n", branchtime_);}

  /// cloning method
  OsiChooseVariable *
  CouenneChooseStrong::clone() const
  {return new CouenneChooseStrong(*this);}

  /// assignment operator
  CouenneChooseStrong&
  CouenneChooseStrong::operator=(const CouenneChooseStrong & rhs)
  {
    if (this != &rhs) {
      Bonmin::BonChooseVariable::operator=(rhs);
      problem_         = rhs.problem_;
      pseudoUpdateLP_  = rhs.pseudoUpdateLP_;
      estimateProduct_ = rhs.estimateProduct_;
      jnlst_           = rhs.jnlst_;
      branchtime_      = rhs.branchtime_;
    }
    return *this;
  }

/***********************************************************************/
int CouenneChooseStrong::goodCandidate(OsiSolverInterface *solver,
				       OsiBranchingInformation *info,
				       OsiObject **object, const int iObject,
				       const double prec) {

  int vInd = object [iObject] -> columnNumber ();

  if (vInd < 0) return 3;  // not a variable object, so deem it good

  bool varIsInt = solver -> isInteger (vInd);

  // all object now are CouenneObjects or derivates (CouenneVarObject,
  // CouenneSOSObject, etc.)

  //CouenneObject    *co    = dynamic_cast <CouenneObject    *> (object [iObject]);
  //OsiSimpleInteger *simpl = dynamic_cast <OsiSimpleInteger *> (object [iObject]);
  
  // int vInd = -1;
  // bool varIsInt = false;
  // if(co) {
  //   vInd = co->Reference()->Index();
  //   if(vInd >= 0) {
  //     varIsInt = co->Reference()->isInteger();
  //   }
  // }
  // else {
  //   if(simpl) {
  //     vInd = object[iObject]->columnNumber();
  //     varIsInt = true;
  //   }
  //   else {

  //     // printf("CouenneChooseStrong::goodCandidate: ### ERROR: unknown object\n");
  //     // exit(1);

  //     // this is probably a SOS object, and anyhow we don't want to
  //     // exit on it.

  //     return 3;
  //   }
  // }
  
  int goodCand = 3; // good candidate
 
  // Must not call branch() for integer variable vInd with
  // upper == lower or for OsiSimpleInteger with
  // info->solution[vInd] not between lower and upper
  if((vInd >= 0) && varIsInt) {
    double vUp = solver->getColUpper()[vInd];
    double vLow = solver->getColLower()[vInd];
    double infoVal = info->solution_[vInd];
    double distToInt = fabs(infoVal - floor(infoVal + 0.5));
    if(distToInt > 0.5) {
      distToInt = 1 - distToInt;
    }
    // if(simpl) {
    //   goodCand = 0; // bad candidate
    //   if((distToInt > info->integerTolerance_) &&
    //      (vUp > vLow + prec)) {
    //     goodCand = 3; // good candidate
    //     if((vUp + prec < infoVal) || (infoVal < vLow - prec)) {
    //       goodCand = 1; // bad candidate
    //     }
    //   }
    // }
    //if(co) {
      goodCand = 2;
      if(vUp > vLow + prec) {
        goodCand = 3; // good candidate
      }
      //}
  }
    
  return(goodCand);
} /* goodCandidate */

/***********************************************************************/
bool CouenneChooseStrong::saveBestCand(OsiObject **object, const int iObject, 
				       const double value, 
				       const double upEstimate, 
				       const double downEstimate,
				       double &bestVal1, 
				       double &bestVal2, int &bestIndex,
				       int &bestWay) {
  bool retval = false;
  if(value > bestVal1) {
    retval = true;
    bestVal1 = value;
    bestIndex = iObject;
    bestWay = upEstimate > downEstimate ? 0 : 1;
    // but override if there is a preferred way
    const OsiObject * obj = object[iObject];
    if (obj->preferredWay() >= 0 && obj->infeasibility()) {
      bestWay = obj->preferredWay();
    }
  }
  return(retval);
} /* saveBestCand */

/*******************************************************************/
  /* Choose a variable. Returns:

    -1  Node is infeasible
     0  Normal termination - we have a candidate
     1  All look satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates

     We can pick up branch from whichObject() and whichWay()
     We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
  */
  int CouenneChooseStrong::chooseVariable (OsiSolverInterface * solver,
					   OsiBranchingInformation *info,
					   bool fixVariables) {

    /// Note: had to copy code from
    /// BonChooseVariable::chooseVariable() in order to test product
    /// thing

    problem_ -> domain () -> push
      (problem_ -> nVars (),
       info -> solution_,
       info -> lower_,
       info -> upper_);

#ifdef TRACE_STRONG2
    int pv = -1;
    if(problem_->doPrint_) {
      if(pv > -1) {
	printf("CCS: beg: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("CCS: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("CCS: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }
    }
 #endif                                                                         

    int retval;
    const double prec = problem_->getFeasTol();

    //int retval = BonChooseVariable::chooseVariable (solver, info, fixVariables);

    // COPY of Bonmin starts here ////////////////////////////////////////////

    // We assume here that chooseVariable is called once at the very
    // beginning with fixVariables set to true.  This is then the root
    // node.

    bool isRoot = isRootNode(info);
    int numberStrong = numberStrong_;

    if (isRoot) {
      numberStrong = CoinMax(numberStrong_, numberStrongRoot_);
    }

    if (numberUnsatisfied_) {
      int cardIndForPseudo = 0, 
	*indForPseudo = new int[numberUnsatisfied_];
      OsiObject ** object = solver->objects();
      const double* upTotalChange = pseudoCosts_.upTotalChange();
      const double* downTotalChange = pseudoCosts_.downTotalChange();
      const int* upNumber = pseudoCosts_.upNumber();
      const int* downNumber = pseudoCosts_.downNumber();
      int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();

      // number of objects to be chosen by doStrongBranching()
      int numberLeft = CoinMin (numberStrong - numberStrongDone_, numberUnsatisfied_);

      results_.clear();
      int returnCode = 0;
      int returnCodeSB = 0;
      bestObjectIndex_ = -1;
      bestWhichWay_ = -1;
      firstForcedObjectIndex_ = -1;
      firstForcedWhichWay_ =-1;
      double bestTrustedVal1 = -COIN_DBL_MAX;
      double bestTrustedVal2 = -COIN_DBL_MAX;

      bool smallGap = false;
      bool sbObjPosImp = false; // true if an object on which strong branching
                                // was performed has a positive improvement 
                                // in both branches; used only when gap is
                                // deemed small

#ifdef USE_SMALL_GAP
      int objInd = problem_ -> Obj (0) -> Body () -> Index ();
      double lbGap = objInd >= 0 ? info -> lower_ [objInd] : problem_ -> Obj (0) -> Body () -> Value ();
      double ubGap = problem_ -> getCutOff ();
      double currentGap = 
	(ubGap >  COUENNE_INFINITY    / 10 ||
	 lbGap < -Couenne_large_bound / 10) ? 1e3 : 
	fabs (ubGap - lbGap) / (1.e-3 + CoinMin (fabs (ubGap), fabs (lbGap)));

      if(currentGap < 1e-3) {
	smallGap = true;
      }
#endif

#ifdef TRACE_STRONG
      if((problem_->doPrint_) && (number_not_trusted_ > 0)) {
	printf("number_not_trusted: %d\n", number_not_trusted_);
      }
#endif

      for (int i=0;i<numberLeft;i++) {
        int iObject = list_[i];
        if (numberBeforeTrusted == 0||
            i < minNumberStrongBranch_ ||
            (
              only_pseudo_when_trusted_ && number_not_trusted_>0 ) ||
              ( !isRoot && (upNumber[iObject]<numberBeforeTrusted ||
                          downNumber[iObject]<numberBeforeTrusted ))||
              ( isRoot && (!upNumber[iObject] && !downNumber[iObject])) ) {
         
#ifdef TRACE_STRONG
	  if(problem_->doPrint_) {
	    printf("Push object %d for strong branch\n", iObject);
	  }
#endif
	  results_.push_back (Bonmin::HotInfo (solver, info, object, iObject));
        }
        else {

#ifdef TRACE_STRONG
	  if(problem_->doPrint_) {
	    printf("Use pseudo cost for object %d\n", iObject);
	  }
#endif
	  indForPseudo[cardIndForPseudo] = iObject;
	  cardIndForPseudo++;
	}
      }

      int numberFixed=0;

      if (results_.size() > 0) {

	//
	// do strong branching
	//

        returnCodeSB = doStrongBranching (solver, info, results_.size(), 1);

        if (bb_log_level_>=3) {
          const char* stat_msg[] = {"NOTDON", "FEAS", "INFEAS", "NOFINI"};
          message(SB_HEADER)<<CoinMessageEol;
          for (unsigned int i = 0; i< results_.size(); i++) {
            double up_change = results_[i].upChange();
            double down_change = results_[i].downChange();
            int up_status = results_[i].upStatus();
            int down_status = results_[i].downStatus();
            message(SB_RES)<<(int) i<<stat_msg[down_status+1]<<down_change
            <<stat_msg[up_status+1]<< up_change<< CoinMessageEol;
          }
        }

        if (returnCodeSB >= 0 && returnCodeSB <= 2) { // 1, 2: some can be fixed
	  if(returnCodeSB > 0) {
	    returnCode = 4; // no bestObject yet
	  }
	  else {
	    returnCode = 0;
	  }
          for (unsigned int i=0;i < results_.size (); i++) {

	    if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0))
	      continue;

            // if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0)) {
	    //   continue;
	    // }

            int iObject = results_[i].whichObject();
	    const OsiObject * obj = object[iObject];
	    int needBranch = goodCandidate(solver, info, object, iObject, prec);

	    ///////////////////////////////////////////////////////////////////

            double upEstimate;

            if (results_[i].upStatus()!=1) {
	      assert (results_[i].upStatus()>=0);
              upEstimate = results_[i].upChange();
            }
            else {
              // infeasible - just say expensive
              if (info->cutoff_<1.0e50)
                upEstimate = 2.0*(info->cutoff_-info->objectiveValue_);
              else
                upEstimate = 2.0*fabs(info->objectiveValue_);
              if (firstForcedObjectIndex_ <0) {
                // first fixed variable
                firstForcedObjectIndex_ = iObject;
                firstForcedWhichWay_ =0;
              }

              numberFixed++;
              if (fixVariables) {
		if(needBranch >= 2) { // for OsiSimpleInteger: do not branch
                  // if upper == lower or if info value is outside bounds
                  // for other objects: branch
		  OsiBranchingObject * branch = 
                                        obj->createBranch(solver, info, 0);
		  branch -> branch (solver);
		  delete branch;
		}
	      }
            }

	    ///////////////////////////////////////////////////////////////////

            double downEstimate;

            if (results_[i].downStatus()!=1) {
	      assert (results_[i].downStatus()>=0);
	      downEstimate = results_[i].downChange();
            }
            else {
              // infeasible - just say expensive
              if (info->cutoff_<1.0e50)
                downEstimate = 2.0*(info->cutoff_-info->objectiveValue_);
              else
                downEstimate = 2.0*fabs(info->objectiveValue_);
              if (firstForcedObjectIndex_ <0) {
                firstForcedObjectIndex_ = iObject;
                firstForcedWhichWay_ =1;
              }
              numberFixed++;
              if (fixVariables) {
		if(needBranch >= 2) { // for OsiSimpleInteger: do not branch
                  // if upper == lower or if info value is outside bounds
                  // for other objects: branch
		  OsiBranchingObject * branch = 
                                         obj->createBranch(solver, info, 1);
		  branch -> branch (solver);
		  delete branch;
		}
	      }
            }
	  
            double
	      MAXMIN_CRITERION = maxminCrit (info),
	      minVal, maxVal, value;

	    if (downEstimate < upEstimate) {
	      minVal = downEstimate;
	      maxVal = upEstimate;
	    } else {
	      minVal = upEstimate;
	      maxVal = downEstimate;
	    }

	    if(smallGap) { // use change in objective value
	      value = minVal;
	    }
	    else {
	      value = 
		estimateProduct_ ? 
		((estProdEps + minVal) * maxVal) :
		(       MAXMIN_CRITERION  * minVal + 
			(1.0 - MAXMIN_CRITERION) * maxVal);
	    }

	    if((needBranch == 3) &&
	       saveBestCand(object, iObject, value, upEstimate, downEstimate,
			    bestTrustedVal1, 
			    bestTrustedVal2, bestObjectIndex_, bestWhichWay_)) {
	      if(returnCodeSB) { // 1 or 2
		returnCode = 2; 
	      }

#ifdef USE_SMALL_GAP
	      if(smallGap && (minVal > 1e-3)) {
		sbObjPosImp = true;
	      }
#endif

	    }
          }
        }
        else { // returnCodeSB is -1 or 3
	  if (returnCodeSB == 3) { // max time - just choose one
	    if(bestObjectIndex_ < 0) {
	      // should not select an integer var with fixed bounds
	      // taken care of below
	      bestObjectIndex_ = list_[0];
	      bestWhichWay_ = 0;
	      returnCode = 0;
	    }
	  }
	  else {
	    returnCode = -1;
	  }
        }
#ifdef OLD_STYLE
	if(bestObjectIndex_ < 0) {
	  bestObjectIndex_ = list_[0];
	}
	cardIndForPseudo = 0; // do not scan other vars using pseudo costs
#endif
      }

      if((returnCodeSB != -1) && 
	 ((returnCode != 0) || (!sbObjPosImp))) {  
	          // if returnCodeSB == -1 (i.e. problem is infeasible)
	          // no need to scan objects with pseudocosts
	          // if returnCode == 0 and sbObjPOsImp is true
	          // we want to branch on that object
                  // to be sure to improve the lower bound
	
	// to keep best object skipped on basis of bounds and branching value
	int bestObjectIndex2 = -1;
	int bestWhichWay2 = 0;
	double bestTrusted2Val1 = -COIN_DBL_MAX;
	double bestTrusted2Val2 = -COIN_DBL_MAX;
	
	for(int ips=0; ips<cardIndForPseudo; ips++) {
	  int iObject = indForPseudo[ips];
	  const OsiObject * obj = object[iObject];
	  int needBranch = goodCandidate(solver, info, object, iObject, prec);

	  double
	    upEstimate       = (upTotalChange   [iObject] * obj -> upEstimate   ()) / upNumber   [iObject],
	    downEstimate     = (downTotalChange [iObject] * obj -> downEstimate ()) / downNumber [iObject],
	    MAXMIN_CRITERION = maxminCrit (info),
	    minVal, maxVal, value;
	  
	  if (downEstimate < upEstimate) {
	    minVal = downEstimate;
	    maxVal = upEstimate;
	  } else {
	    minVal = upEstimate;
	    maxVal = downEstimate;
	  }
	  
	  value = 
	    estimateProduct_ ? 
	    ((estProdEps + minVal) * maxVal) :
	    (       MAXMIN_CRITERION  * minVal + 
		    (1.0 - MAXMIN_CRITERION) * maxVal);
	  
	  
	  // store bad candidates in secondary best
	  if(needBranch != 3) {
	    if(saveBestCand(object, iObject, value, 
			    upEstimate, downEstimate,
			    bestTrusted2Val1, 
			    bestTrusted2Val2, bestObjectIndex2, 
			    bestWhichWay2)) {
	      // no returnCode change
	    }

#ifdef TRACE_STRONG
            if(problem_->doPrint_) {
              printf("Object %d skip pseudocost\n", iObject);
            }
#endif
	  }
	  else {
	    if(saveBestCand(object, iObject, value, 
			    upEstimate, downEstimate,
			    bestTrustedVal1, 
			    bestTrustedVal2, bestObjectIndex_, bestWhichWay_)) {
	      returnCode = (returnCode ? 3 : 0); // if returnCode was 2 or 3
                                                 // it becomes 3
	    }

#ifdef TRACE_STRONG
            if(problem_->doPrint_) {
              printf("Object %d use pseudocost\n", iObject);
            }
#endif
	  }
	}
    
	if((bestObjectIndex_ < 0) && (bestObjectIndex2 >= 0)) {
	  bestObjectIndex_ = bestObjectIndex2;
	  bestWhichWay_ = bestWhichWay2;
	  bestTrustedVal1 = bestTrusted2Val1;
	  returnCode = 4;
	}
      }

      int objectVarInd = -1;
      if(bestObjectIndex_ >= 0) {
        OsiObject * obj = object[bestObjectIndex_];
        obj->setWhichWay(bestWhichWay_);
	// CouenneObject *co =  dynamic_cast <CouenneObject *>(object[bestObjectIndex_]);
	// if(co) {
	//   objectVarInd = co->Reference()->Index();
	// }
	// else {
	objectVarInd = obj->columnNumber();
	//}

        message(BRANCH_VAR)<<bestObjectIndex_<< bestWhichWay_ <<CoinMessageEol;

#ifdef TRACE_STRONG
	if(problem_->doPrint_) {
	  if(objectVarInd >= 0) {
            double vUb = solver->getColUpper()[objectVarInd];
            double vLb = solver->getColLower()[objectVarInd];
            double vSolI = info->solution_[objectVarInd];
            double vSolS = solver->getColSolution()[objectVarInd];
            printf("Branch on object %d (var: %d): solInfo: %10.4f  SolSolver: %10.4f low: %10.4f  up: %10.4f\n",
                   bestObjectIndex_, objectVarInd, vSolI, vSolS, vLb, vUb);
	  }
	  else {
	    printf("Branch on object %d (var: -1)\n", bestObjectIndex_);
	  }
	}
#endif
      }
      message(CHOSEN_VAR)<<bestObjectIndex_<<CoinMessageEol;
    
      if (numberFixed==numberUnsatisfied_&&numberFixed)
        returnCode=4;

      if((returnCode == 2) || (returnCode == 3)) {
        if((objectVarInd > -1) && 
	   (goodCandidate(solver, info, object, bestObjectIndex_, prec) != 2)) {
          // Can occur: two objects for same var, first scanned object
          // has both branches feasible and is saved as bestObjectIndex_,
          // second scanned object fixes the variable
          returnCode = 4;
        }
      }
      retval = returnCode;

      delete[] indForPseudo;
    }
    else {
      retval = 1;
    }

    // COPY of Bonmin ends here //////////////////////////////////////////////


#ifdef TRACE_STRONG2
    if(problem_->doPrint_) {
      if(pv > -1) {
	printf("CCS: end: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("CCS: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("CCS: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }
    }
#endif

#ifdef TRACE_STRONG
    if(problem_->doPrint_) {
      printf("CouenneChooseStrong::ChooseVariable(): retval: %d\n", retval);
    }
#endif
  
    problem_ -> domain () -> pop ();

    return retval;
  }

void eliminateIntegerObjects (OsiSolverInterface *model);
void eliminateIntegerObjects (CbcModel           *model);

/*******************************************************************/
  // Sets up strong list and clears all if initialize is true.
  // Returns number of infeasibilities.
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {

    static bool
      firstCall = true,
      warned    = false;

    if (firstCall) {

      eliminateIntegerObjects (const_cast <OsiSolverInterface *> (solver_));
      eliminateIntegerObjects (const_cast <OsiSolverInterface *> (info -> solver_));

      firstCall = false;
    }

    initialize = true; // to avoid failed assert in BonChooseVariable::setupList()

    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       info -> solution_, 
       info -> lower_, 
       info -> upper_); // have to alloc+copy

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "----------------- (strong) setup list\n");

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
	printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
		info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
    }

#ifdef TRACE_STRONG
    OsiObject ** object = info->solver_->objects();
    if(problem_->doPrint_) {
      printObjViol(info);
    }
#endif

    // int way;
    // for (int i=0; i<info->solver_->numberObjects(); ++i)
    //   printf ("[%d:%d,%g] ", 
    // 	      info -> solver_ -> objects () [i] -> columnNumber (), 
    // 	      info -> solver_ -> objects () [i] -> priority (), 
    // 	      info -> solver_ -> objects () [i] -> infeasibility (info,way));
    // printf ("\n");

    //
    // Real list setup
    //

    int retval = gutsOfSetupList (info, initialize);

    if (retval == 0) { // No branching is possible

#ifdef FM_CHECKNLP2
      if(!(problem_->checkNLP2(info->solution_, 
			       info->objectiveValue_, true, // care about obj
			       false, // do not stop at first viol 
			       true, // checkAll
			       problem_->getFeasTol()))) {
                                // false for NOT stopping at first violation
	if (!warned) {
	  printf("CouenneChooseStrong::setupList(): ### WARNING: checkNLP2() returns infeasible, no branching object selected\n");
          warned = true;
	}
      }
#else /* not FM_CHECKNLP2 */
      double ckObj = info->objectiveValue_;
      if(!(problem_->checkNLP(info->solution_, ckObj, true))) {
	if (!warned) {
	  printf("CouenneChooseStrong::setupList(): ### WARNING: checkNLP() returns infeasible, no branching object selected\n");
          warned = true;
	}
      }
#endif /* not FM_CHECKNLP2 */
    	
#ifdef FM_TRACE_OPTSOL
#ifdef FM_CHECKNLP2
      problem_->getRecordBestSol()->update();
#else /* not FM_CHECKNLP2 */
      problem_->getRecordBestSol()->update(info->solution_, problem_->nVars(),
					   ckObj, problem_->getFeasTol());
#endif /* not FM_CHECKNLP2 */
#endif
    }

#ifdef TRACE_STRONG
    if(problem_->doPrint_) {
      printf("Strong list: (obj_ind var_ind priority useful)\n");
      printf("numberStrong: %d  numberStrongRoot: %d  retval: %d\n", 
	     numberStrong_, numberStrongRoot_, retval);
      for(int i=0; i<retval; i++) {
	// CouenneObject *co =  dynamic_cast <CouenneObject *>(object[list_[i]]);
	int objectInd = -1;
	// if(co) {
	//   objectInd = co->Reference()->Index();
	// }
	// else {
	  objectInd = object[list_[i]]->columnNumber();
	// }
	printf(" (%d %d %d %6.4f)", list_[i], objectInd, 
	       object[list_[i]]->priority(), useful_[i]);
      }
      printf("\n");
    }
#endif

    // for (int i=0; i < (numberStrong_ ? CoinMin (numberStrong_, solver_ -> numberObjects ()) : 1); i++) {
    //   printf ("list %3d: %3d ", i, list_ [i]);
    //   if (!((i+1) % 12)) printf ("\n");
    // }

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "----------------- (strong) setup list done - %d infeasibilities\n", retval);

    problem_ -> domain () -> pop ();
    return retval;
  }

/****************************************************************************/
  /// Add list of options to be read from file ////////////////////////////////////////
  void CouenneChooseStrong::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

    roptions -> AddStringOption6
      ("pseudocost_mult",
       "Multipliers of pseudocosts for estimating and update estimation of bound",
       "interval_br_rev",

       "infeasibility", "infeasibility returned by object",

       "projectDist",   "distance between current LP point and resulting branches' LP points",

       "interval_lp",   "width of the interval between bound and current lp point",
       "interval_lp_rev",   "similar to interval_lp, reversed",

       "interval_br",   "width of the interval between bound and branching point",
       "interval_br_rev",   "similar to interval_br, reversed");

    roptions -> AddStringOption2
      ("pseudocost_mult_lp",
       "Use distance between LP points to update multipliers of pseudocosts "  
       "after simulating branching",
       "no",
       "yes", "",
       "no",  "");

    roptions -> AddStringOption2
      ("estimate_select",
       "How the min/max estimates of the subproblems' bounds are used in strong branching",
       "normal",
       "normal",   "as usual in literature",
       "product",  "use their product");

    roptions -> AddStringOption2
      ("trust_strong",
       "Fathom strong branching LPs when their bound is above the cutoff",
       "yes",
       "yes", "",
       "no",  "");
  }


  // Returns true if solution looks feasible against given objects
  bool CouenneChooseStrong::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects) {

#ifdef FM_CHECKNLP2
    return problem_ -> checkNLP2 (solution, 0, false, true, true, 
				  problem_->getFeasTol());
#else
    int indobj = problem_ -> Obj (0) -> Body () -> Index ();
    return problem_ -> checkNLP (solution, indobj >= 0 ? solution [indobj] : problem_ -> Obj (0) -> Body () -> Value ());
#endif
  }

/****************************************************************************/
  void CouenneChooseStrong::printObjViol(OsiBranchingInformation *info) {

    OsiObject ** object = info->solver_->objects();
    int numberObjects = info->solver_->numberObjects();

    printf("CouenneChooseStrong::printObjViol(): Object violations: (obj_ind  var_ind  violation)");
    double maxViol = 0;
    double minPosViol = 1e50;
    for(int i=0; i<numberObjects; i++) {
      //CouenneObject *co =  dynamic_cast <CouenneObject *>(object[i]);
      int indVar = -1;
      // if(co) {
      // 	indVar = co->Reference()->Index();
      // }
      // else {
      indVar = object[i]->columnNumber();
      // }
      int way;
      double value = object[i]->infeasibility(info,way);
      //double value = object[i]->checkInfeasibility(info);
      maxViol = (value > maxViol ? value : maxViol);
      if(value > 0.0) {
	printf("(%d %d %f)", i, indVar, value);
	minPosViol = (value < minPosViol ? value : minPosViol);
      }
    }
    printf("\nmaxViol: %g  minPosViol: %g\n", maxViol, minPosViol);

  } /* printObjViol */

//}
