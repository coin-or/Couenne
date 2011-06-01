/* $Id$
 *
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 *          Pietro Belotti, Lehigh University
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneObject.hpp"
#include "BonChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneBranchingObject.hpp"

#include "CouenneRecordBestSol.hpp"

//#define TRACE_STRONG
//#define TRACE_STRONG2
//#define FM_MOD
//#define FM_SORT_STRONG
//#define FM_ALWAYS_SORT
//#define OLD_STYLE

using namespace Couenne;

const CouNumber estProdEps = 1e-6;


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

    minDepthPrint_ = -1;
  }

  /// copy constructor
  CouenneChooseStrong::CouenneChooseStrong (const CouenneChooseStrong& rhs) :
    Bonmin::BonChooseVariable (rhs),
    problem_          (rhs.problem_),
    pseudoUpdateLP_   (rhs.pseudoUpdateLP_),
    estimateProduct_  (rhs.estimateProduct_),
    jnlst_            (rhs.jnlst_),
    branchtime_       (rhs.branchtime_),
    minDepthPrint_    (rhs.minDepthPrint_)
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
      minDepthPrint_   = rhs.minDepthPrint_;
    }
    return *this;
  }


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
#ifdef FM_MOD
/*******************************************************************/
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
    if(info->depth_ > minDepthPrint_) {
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
      int cardIndForPseudo = 0, *indForPseudo = new int[numberUnsatisfied_];
      OsiObject ** object = solver->objects();
      const double* upTotalChange = pseudoCosts_.upTotalChange();
      const double* downTotalChange = pseudoCosts_.downTotalChange();
      const int* upNumber = pseudoCosts_.upNumber();
      const int* downNumber = pseudoCosts_.downNumber();
      int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
      int numberLeft = CoinMin(numberStrong -numberStrongDone_,numberUnsatisfied_);
      results_.clear();
      int returnCode = 0;
      int returnCodeSB = 0;
      bestObjectIndex_ = -1;
      bestWhichWay_ = -1;
      firstForcedObjectIndex_ = -1;
      firstForcedWhichWay_ =-1;
      double bestTrusted=-COIN_DBL_MAX;
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
	  if(info->depth_ > minDepthPrint_) {
	    printf("Push object %d for strong branch\n", iObject);
	  }
#endif
	  results_.push_back(Bonmin::HotInfo(solver, info,
                                object, iObject));
        }
        else {

#ifdef TRACE_STRONG
	  if(info->depth_ > minDepthPrint_) {
	    printf("Use pseudo cost for object %d\n", iObject);
	  }
#endif
	  indForPseudo[cardIndForPseudo] = iObject;
	  cardIndForPseudo++;
	}
      }

      int numberFixed=0;

      if (results_.size() > 0) {

	// do strong branching //////////////////////////////////////////////////
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
	  const double prec = problem_->getFeasTol();

          for (unsigned int i=0;i < results_.size (); i++) {

	    if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0))
	      continue;


            if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0)) {
	      continue;
	    }

            int iObject = results_[i].whichObject();

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
		bool needBranch = true;
		const OsiObject * obj = object[iObject];
		CouenneObject *co = dynamic_cast <CouenneObject *>(object[iObject]);
		OsiSimpleInteger *simpl = dynamic_cast <OsiSimpleInteger *>(object[iObject]);
               int vInd = -1;
                if(co) {
                  vInd = co->Reference()->Index();
                }
                else {
                  vInd = obj->columnNumber();
                }
                if((vInd >= 0) && (simpl)) {
                  needBranch = false;
                  double nearest = floor(info->solution_[vInd] + 0.5);
                  if(nearest > 0.5) {
                    nearest = 1 - nearest;
                  }
                                                                                
                  if((nearest > info->integerTolerance_) &&
                     (solver->getColUpper()[vInd] > solver->getColLower()[vInd]\
 + prec)) {
                    needBranch = true;
                  }
                }
                if(needBranch) {
		  OsiBranchingObject * branch = obj -> createBranch (solver, info, 0);
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
		bool needBranch = true;
		const OsiObject * obj = object[iObject];
		CouenneObject *co = dynamic_cast <CouenneObject *>(object[iObject]);
		OsiSimpleInteger *simpl = dynamic_cast <OsiSimpleInteger *>(object[iObject]);
		int vInd = -1;
                if(co) {
                  vInd = co->Reference()->Index();
                }
                else {
                  vInd = obj->columnNumber();
                }
                if((vInd >= 0) && (simpl)) {
                  needBranch = false;
                  double nearest = floor(info->solution_[vInd] + 0.5);
                  if(nearest > 0.5) {
                    nearest = 1 - nearest;
                  }
                                                                                
                  if((nearest > info->integerTolerance_) &&
                     (solver->getColUpper()[vInd] > solver->getColLower()[vInd]\
 + prec)) {
                    needBranch = true;
                  }
                }
                if(needBranch) {
		  OsiBranchingObject * branch = obj -> createBranch (solver, info, 1);
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

	    value = 
	      estimateProduct_ ? 
	      ((estProdEps + minVal) * maxVal) :
	      (       MAXMIN_CRITERION  * minVal + 
	       (1.0 - MAXMIN_CRITERION) * maxVal);

	    if (value>bestTrusted) {
              bestTrusted = value;
              bestObjectIndex_ = iObject;
              bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
              // but override if there is a preferred way
              const OsiObject * obj = object[iObject];
              if (obj->preferredWay()>=0&&obj->infeasibility())
                bestWhichWay_ = obj->preferredWay();
              if(returnCodeSB) { // 1 or 2
                returnCode = 2; 
	      }
	    }
          }
        }
        else { // returnCodeSB is -1 or 3
	  if (returnCodeSB == 3) { // max time - just choose one
	    if(bestObjectIndex_ < 0) {
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

      if(returnCodeSB != -1) {  
	          // if returnCodeSB == -1 (i.e. problem is infeasible)
	          // no need to scan objects with pseudocosts
	const double *solverUb = solver->getColUpper();
	const double *solverLb = solver->getColLower();
	const double prec = problem_->getFeasTol();
	
	// to keep best object skipped on basis of bounds and branching value
	int bestObjectIndex2 = -1;
	int bestWhichWay2 = 0;
	double bestTrusted2 = -COIN_DBL_MAX;
	
	for(int ips=0; ips<cardIndForPseudo; ips++) {
	  int iObject = indForPseudo[ips];
	  const OsiObject * obj = object[iObject];
	  double
	    upEstimate       = (upTotalChange[iObject]*obj->upEstimate())/upNumber[iObject],
	    downEstimate     = (downTotalChange[iObject]*obj->downEstimate())/downNumber[iObject],
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
	  
	  
	  // skip OsiSimpleInteger objects with variable fixed or 
	  // branching value outside bounds
	  bool skipIt = false;
	  CouenneObject *co = dynamic_cast <CouenneObject *>(object[iObject]);
	  OsiSimpleInteger *simpl = dynamic_cast <OsiSimpleInteger *>(object[iObject]);
	  int vInd = -1;
	  if(co) {
	    vInd = co->Reference()->Index();
	  }
	  else {
	    vInd = obj->columnNumber();
	  }
	  if(simpl && (vInd >= 0)) {
	    double vUb = solverUb[vInd];
	    double vLb = solverLb[vInd];
	    double vSol = info->solution_[vInd];
	    if((vSol < vLb + prec) || (vSol > vUb - prec) || (vUb-vLb < prec)) {
	      skipIt = true;
	      numberFixed++;
	    }
	  }
	  
	  if(skipIt) {
	    if (value > bestTrusted2) {
	      bestObjectIndex2 = iObject;
	      bestWhichWay2 = upEstimate>downEstimate ? 0 : 1;
	      bestTrusted2 = value;
	    }
#ifdef TRACE_STRONG
            if(info->depth_ > minDepthPrint_) {
              printf("Object %d skip pseudocost\n", iObject);
            }
#endif
	  }
	  else {
	    if (value > bestTrusted) {
	      bestObjectIndex_ = iObject;
	      bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
	      bestTrusted = value;
	      returnCode = (returnCode ? 3 : 0); // if returnCode was 2 or 3
                                                 // it becomes 3
	    }
#ifdef TRACE_STRONG
            if(info->depth_ > minDepthPrint_) {
              printf("Object %d use pseudocost\n", iObject);
            }
#endif
	  }
	}
    
	delete[] indForPseudo;

	if((bestObjectIndex_ < 0) && (bestObjectIndex2 >= 0)) {
	  bestObjectIndex_ = bestObjectIndex2;
	  bestWhichWay_ = bestWhichWay2;
	  bestTrusted = bestTrusted2;
	  returnCode = 4;
	}
      }

      int objectVarInd = -1;
      if(bestObjectIndex_ >= 0) {
        OsiObject * obj = object[bestObjectIndex_];
        obj->setWhichWay(bestWhichWay_);
	CouenneObject *co =  dynamic_cast <CouenneObject *>(object[bestObjectIndex_]);
	if(co) {
	  objectVarInd = co->Reference()->Index();
	}
	else {
	  objectVarInd = obj->columnNumber();
	}

        message(BRANCH_VAR)<<bestObjectIndex_<< bestWhichWay_ <<CoinMessageEol;

#ifdef TRACE_STRONG
	if(info->depth_ > minDepthPrint_) {
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
	   (solver->getColUpper()[objectVarInd] < solver->getColLower()[objectVarInd] + problem_->getFeasTol())) {
          // Can occur: two objects for same var, first scanned object
          // has both branches feasible and is saved as bestObjectIndex_,
          // second scanned object fixes the variable
          returnCode = 4;
        }
      }
      retval = returnCode;
    }
    else {
      retval = 1;
    }

    // COPY of Bonmin ends here //////////////////////////////////////////////


#ifdef TRACE_STRONG2
    if(info->depth_ > minDepthPrint_) {
      if(pv > -1) {
	printf("CCS: end: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("CCS: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("CCS: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }
      printf("retval: %d\n", retval);
    }
#endif
  
    problem_ -> domain () -> pop ();

    return retval;
  }
/*******************************************************************/
#else /* not FM_MOD */
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

    int retval;

#ifdef TRACE_STRONG2
    int pv = -1;
    if(info->depth_ > minDepthPrint_) {
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
      const double* upTotalChange = pseudoCosts_.upTotalChange();
      const double* downTotalChange = pseudoCosts_.downTotalChange();
      const int* upNumber = pseudoCosts_.upNumber();
      const int* downNumber = pseudoCosts_.downNumber();
      int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
      int numberLeft = CoinMin(numberStrong -numberStrongDone_,numberUnsatisfied_);
      results_.clear();
      int returnCode=0;
      bestObjectIndex_ = -1;
      bestWhichWay_ = -1;
      firstForcedObjectIndex_ = -1;
      firstForcedWhichWay_ =-1;
      double bestTrusted=-COIN_DBL_MAX;
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
	  if(info->depth_ > minDepthPrint_) {
	    printf("Push object %d for strong branch\n", iObject);
	  }
#endif
	  results_.push_back(Bonmin::HotInfo(solver, info,
					     solver->objects(), iObject));
        }
        else {

#ifdef TRACE_STRONG
	  if(info->depth_ > minDepthPrint_) {
	    printf("Use pseudo cost for object %d\n", iObject);
	  }
#endif

          const OsiObject * obj = solver->object(iObject);
          double
	    upEstimate       = (upTotalChange[iObject]*obj->upEstimate())/upNumber[iObject],
	    downEstimate     = (downTotalChange[iObject]*obj->downEstimate())/downNumber[iObject],
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

          if (value > bestTrusted) {
            bestObjectIndex_=iObject;
            bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
            bestTrusted = value;
          }
        }
      }

      int numberFixed=0;

      if (results_.size() > 0) {

	// do strong branching //////////////////////////////////////////////////
        returnCode = doStrongBranching (solver, info, results_.size(), 1);

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

        if (returnCode >= 0 && 
	    returnCode <= 2) {

          if (returnCode)
            returnCode = (bestObjectIndex_>=0) ? 3 : 4;

          for (unsigned int i=0;i < results_.size (); i++) {

	    if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0))
	      continue;


            if((results_[i].upStatus() < 0) || (results_[i].downStatus() < 0)) {
	      continue;
	    }

            int iObject = results_[i].whichObject();
	    
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

                const OsiObject * obj = solver->objects()[iObject];
		OsiBranchingObject * branch = obj -> createBranch (solver, info, 0);
		branch -> branch (solver);
		delete branch;
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

                const OsiObject * obj = solver->objects()[iObject];
		OsiBranchingObject * branch = obj -> createBranch (solver, info, 1);
		branch -> branch (solver);
		delete branch;
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

	    value = 
	      estimateProduct_ ? 
	      ((estProdEps + minVal) * maxVal) :
	      (       MAXMIN_CRITERION  * minVal + 
	       (1.0 - MAXMIN_CRITERION) * maxVal);

	    if (value>bestTrusted) {
              bestTrusted = value;
              bestObjectIndex_ = iObject;
              bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
              // but override if there is a preferred way
              const OsiObject * obj = solver->object(iObject);
              if (obj->preferredWay()>=0&&obj->infeasibility())
                bestWhichWay_ = obj->preferredWay();
              if (returnCode)
                returnCode = 2;
            }
          }
        }
        else if (returnCode==3) {
          // max time - just choose one
          bestObjectIndex_ = list_[0];
          bestWhichWay_ = 0;
          returnCode=0;
        }
      }
      else {
        bestObjectIndex_=list_[0];
      }
    
      if(bestObjectIndex_ >=0) {
        OsiObject * obj = solver->objects()[bestObjectIndex_];
        obj->setWhichWay(bestWhichWay_);
        message(BRANCH_VAR)<<bestObjectIndex_<< bestWhichWay_ <<CoinMessageEol;
	CouenneObject *co =  dynamic_cast <CouenneObject *>(solver->objects()[bestObjectIndex_]);

	int objectInd = -1;
	if(co) {
	  objectInd = co->Reference()->Index();
	}
	else {
	  objectInd = obj->columnNumber();
	}

#ifdef TRACE_STRONG
	if(info->depth_ > minDepthPrint_) {
	  if(objectInd >= 0) {
            double vUb = solver->getColUpper()[objectInd];                      
            double vLb = solver->getColLower()[objectInd];                      
            double vSolI = info->solution_[objectInd];                          
            double vSolS = solver->getColSolution()[objectInd];                 
            printf("Branch on object %d (var: %d): solInfo: %10.4f  SolSolver: \
%10.4f low: %10.4f  up: %10.4f\n",                                              
                   bestObjectIndex_, objectInd, vSolI, vSolS, vLb, vUb);        
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
      retval = returnCode;
    }
    else {
      retval = 1;
    }

    // COPY of Bonmin ends here //////////////////////////////////////////////

#ifdef TRACE_STRONG2
    if(info->depth_ > minDepthPrint_) {
      if(pv > -1) {
	printf("CCS: end: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("CCS: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("CCS: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }                                                   
      printf("retval: %d\n", retval);
    }
#endif
    problem_ -> domain () -> pop ();

    return retval;
    }
#endif /* not FM_MOD */

  // Sets up strong list and clears all if initialize is true.
  // Returns number of infeasibilities.
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {
    static bool warned = false;

    initialize = true; // to avoid failed assert in BonChooseVariable::setupList()

    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       info -> solution_, 
       info -> lower_, 
       info -> upper_); // have to alloc+copy

#ifdef COIN_HAS_NTY

    if (problem_ -> orbitalBranching ()) {

      problem_ -> ChangeBounds (info -> lower_,
				info -> upper_, 
				problem_ -> nVars ());
    
      problem_ -> Compute_Symmetry();
    }

#endif


    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "----------------- (strong) setup list\n");

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
	printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
		info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
    }

    OsiObject ** object = info->solver_->objects();
    int numberObjects = info->solver_->numberObjects();

#ifdef TRACE_STRONG
    if(info->depth_ > minDepthPrint_) {
      printObjViol(info);
    }
#endif

    int retval = gutsOfSetupList(info, initialize);

    if(retval == 0) { // No branching is possible
      double ckObj = info->objectiveValue_;

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
    if(info->depth_ > minDepthPrint_) {
      printf("Strong list: (obj_ind var_ind priority useful)\n");
      printf("numberStrong: %d  numberStrongRoot: %d  retval: %d\n", 
	     numberStrong_, numberStrongRoot_, retval);
      for(int i=0; i<retval; i++) {
	CouenneObject *co =  dynamic_cast <CouenneObject *>(object[list_[i]]);
	int objectInd = -1;
	if(co) {
	  objectInd = co->Reference()->Index();
	}
	else {
	  objectInd = object[list_[i]]->columnNumber();
	}
	printf(" (%d %d %d %6.4f)", list_[i], objectInd, 
	       object[list_[i]]->priority(), useful_[i]);
      }
      printf("\n");
    }
#endif
  

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "----------------- (strong) setup list done - %d infeasibilities\n", retval);

    problem_ -> domain () -> pop ();
    return retval;
  }

/****************************************************************************/
// Copied from BonChooseVariable.cpp and modified slightly
  int CouenneChooseStrong::gutsOfSetupList(OsiBranchingInformation *info, 
				      bool initialize)
  {
    if (numberBeforeTrustedList_ < 0) {
      number_not_trusted_ = 1;
      printf("CouenneChooseStrong::gutsOfSetupList(): Did not think we were using this; Please double check ...\n");
      exit(1);
      return OsiChooseVariable::setupList(info, initialize);
    }
    if (initialize) {
      status_=-2;
      delete [] goodSolution_;
      bestObjectIndex_=-1;
      numberStrongDone_=0;
      numberStrongIterations_ = 0;
      numberStrongFixed_ = 0;
      goodSolution_ = NULL;
      goodObjectiveValue_ = COIN_DBL_MAX;
      number_not_trusted_=0;
    }
    else {
      throw CoinError(CNAME,"setupList","Should not be called with initialize==false");
    }
    numberOnList_=0;
    numberUnsatisfied_=0;
    int numberObjects = solver_->numberObjects();
    assert (numberObjects);
    if (numberObjects>pseudoCosts_.numberObjects()) {
      //std::cout<<"Number objects "<<numberObjects<<std::endl;
      //AW : How could that ever happen?
      //PB : It happens for instance when SOS constraints are added. They are added after the creation of this.
      //   assert(false && "Right now, all old content is deleted!");
      // redo useful arrays
      int saveNumberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
      pseudoCosts_.initialize(numberObjects);
      pseudoCosts_.setNumberBeforeTrusted(saveNumberBeforeTrusted);
    }
    double check = -COIN_DBL_MAX;
    int checkIndex=0;
    int bestPriority=COIN_INT_MAX;
    int putOther = numberObjects;
    int i;

#ifdef FM_SORT_STRONG
    int numStr = numberStrong_;
    if(isRootNode(info)) {
      numStr = numberStrongRoot_;
    }
    int maximumStrong = CoinMin(numStr, numberObjects) ;
    int lastPrio = problem_->getLastPrioSort();
    int card_vPriority = 0;
    int posEnd_vPriority = numberObjects;
    double *vPriority = new double[numberObjects];
#else /* not FM_SORT_STRONG */
    int maximumStrong = CoinMin(CoinMax(numberStrong_,numberStrongRoot_),
        numberObjects) ;
    for (i=0;i<numberObjects;i++) {
      list_[i]=-1;
      useful_[i]=0.0;
    }
    // We make a second list for most fractional variables
    int* list2 = NULL;
    double* useful2 = NULL;
    double check2 = -COIN_DBL_MAX;
    int checkIndex2=0;
    int max_most_fra = setup_pseudo_frac_ > 0. ? (int)floor(setup_pseudo_frac_*(double)maximumStrong): 0;
    if (setup_pseudo_frac_ > 0.) {
      max_most_fra = CoinMax(1, max_most_fra);
    }
    if (max_most_fra) {
      list2 = new int[max_most_fra];
      useful2 = new double[max_most_fra];
      for (i=0;i<max_most_fra;i++) {
        list2[i]=-1;
        useful2[i]=0.0;
      }
    }
#endif /* not FM_SORT_STRONG */

#ifdef FM_CHECK
    const double* upTotalChange = pseudoCosts_.upTotalChange();
    const double* downTotalChange = pseudoCosts_.downTotalChange();
    int pseudoNum = pseudoCosts_.numberObjects();
    for(i=0; i<pseudoNum; i++) {
      if(isnan(upTotalChange[i]) || isinf(upTotalChange[i])) {
	printf("CouenneChooseStrong::gutsOfSetupList(): upTotalChange[%d]: not a number or infinite\n", i);
	exit(1);
      }
      if(isnan(downTotalChange[i]) || isinf(downTotalChange[i])) {
	printf("CouenneChooseStrong::gutsOfSetupList(): downTotalChange[%d]: not a number or infinite\n", i);
	exit(1);
      }
    }
#endif

    OsiObject ** object = info->solver_->objects();
    double upMultiplier, downMultiplier;
    computeMultipliers(upMultiplier, downMultiplier);

    // Say feasible
    bool feasible = true;
    const double MAXMIN_CRITERION = maxminCrit(info);

    bool firstPass = false; // not important; useful for making two
                            // passes, picking different objects
    while(numberOnList_ == 0) {
      for(i=0;i<numberObjects;i++) {
	int way;
	double value = object[i]->infeasibility(info, way);
	double lbForInfeas = 0.0;
	if(value > lbForInfeas) {
	  numberUnsatisfied_++;
	  if(value >= 1e50) {
	    // infeasible
	    feasible=false;
	    break;
	  }
	  int priorityLevel = object[i]->priority();

#ifdef FM_SORT_STRONG
	  if(priorityLevel > lastPrio) {
	    posEnd_vPriority--;
	    vPriority[posEnd_vPriority] = priorityLevel;
	    list_[posEnd_vPriority] = i;
	  }
	  else {
	    vPriority[card_vPriority] = priorityLevel;
	    list_[card_vPriority] = i;
	    card_vPriority++;
	  }
#else /* not FM_SORT_STRONG */
	  // Better priority? Flush choices.
	  if(priorityLevel < bestPriority) {
	    for (int j=maximumStrong-1; j>=0; j--) {
	      if(list_[j] >= 0) {
		int iObject = list_[j];
		list_[j]=-1;
		useful_[j]=0.0;
		list_[--putOther]=iObject;
	      }
	    }
	    maximumStrong = CoinMin(maximumStrong,putOther);
	    bestPriority = priorityLevel;
	    check=-COIN_DBL_MAX;
	    checkIndex=0;
	    check2=-COIN_DBL_MAX;
	    checkIndex2=0;
	    number_not_trusted_=0;
	    if(max_most_fra > 0) {
	      for(int j=0; j<max_most_fra; j++) {
		list2[j]=-1;
		useful2[j]=0.0;
	      }
	    }
	  }
	  if(priorityLevel == bestPriority) {
	    // Modify value
	    double value2;
	    value = computeUsefulness(MAXMIN_CRITERION,
				      upMultiplier, downMultiplier, value,
				      object[i], i, value2);
	    if(value > check) {
	      //add to list
	      int iObject = list_[checkIndex];
	      if(iObject >= 0) {
		assert (list_[putOther-1]<0);
		list_[--putOther]=iObject;  // to end
	      }
	      list_[checkIndex]=i;
	      assert (checkIndex<putOther);
	      useful_[checkIndex]=value;
	      // find worst
	      check=COIN_DBL_MAX;
	      maximumStrong = CoinMin(maximumStrong,putOther);
	      for (int j=0; j<maximumStrong; j++) {
		if(list_[j]>=0) {
		  if (useful_[j]<check) {
		    check=useful_[j];
		    checkIndex=j;
		  }
		}
		else {
		  check=0.0;
		  checkIndex = j;
		  break;
		}
	      }
	    }
	    else {
	      // to end
	      assert (list_[putOther-1]<0);
	      list_[--putOther]=i;
	      maximumStrong = CoinMin(maximumStrong,putOther);
	    }
	    if(max_most_fra > 0 && value2 > check2) {
	      // add to list of integer infeasibilities
	      number_not_trusted_++;
	      list2[checkIndex2]=i;
	      useful2[checkIndex2]=value2;
	      // find worst
	      check2=COIN_DBL_MAX;
	      for(int j=0; j<max_most_fra; j++) {
		if(list2[j] >= 0) {
		  if(useful2[j] < check2) {
		    check2=useful2[j];
		    checkIndex2=j;
		  }
		}
		else {
		  check2=0.0;
		  checkIndex2 = j;
		  break;
		}
	      }
	    }
	  }
	  else {
	    // worse priority
	    // to end
	    assert (list_[putOther-1]<0);
	    list_[--putOther]=i;
	    maximumStrong = CoinMin(maximumStrong,putOther);
	  }
#endif /* not FM_SORT_STRONG */
	}
      }

#ifdef FM_SORT_STRONG

#ifdef FM_CHECK
      if(card_vPriority - posEnd_vPriority + numberObjects != numberUnsatisfied_) {
	printf("CouenneChooseStrong::gutsOfSetupList(): ### ERROR: card_vPriority: %d  posEnd_vPriority: %d  numberUnsatisfied: %d numberObjects: %d\n",
	       card_vPriority, posEnd_vPriority, numberUnsatisfied_, numberObjects);
	exit(1);
      }
#endif

      numberOnList_ = 0;
      if(feasible) {
	int card_smallerThanPrio = card_vPriority;
	if(posEnd_vPriority > card_vPriority) {
	  for(i=posEnd_vPriority; i<numberObjects; i++) {
	    list_[card_vPriority] = list_[i];
	    list_[i] = -1;
	    vPriority[card_vPriority] = vPriority[i]; 
	    card_vPriority++;
	  }
	}
	else {
	  card_vPriority = numberUnsatisfied_;
	}
	// correct bounds if card_smallThanPrio >= maximumStrong
	int sortFrom = 0;
	int sortUpTo = card_smallerThanPrio;
	if(card_smallerThanPrio < maximumStrong) {
	  sortFrom = card_smallerThanPrio;
	  sortUpTo = card_vPriority;
	}
	if(card_vPriority > 0) {
	  numberOnList_ = (card_vPriority < maximumStrong ? card_vPriority : maximumStrong);

#ifdef FM_ALWAYS_SORT
	  bool alwaysSort = true;
#else	  
	  bool alwaysSort = false;
#endif
	  if(alwaysSort) {
	    sortFrom = 0;
	    sortUpTo = card_vPriority;
	  }
	  if((sortUpTo > maximumStrong) || alwaysSort){
	    // sort list_[card_sortFrom..card_sortUpTo-1] according to priority
	    CoinSort_2(vPriority + sortFrom, vPriority + sortUpTo, 
		       list_ + sortFrom);
	  }
	  for(i=0; i<card_vPriority; i++) {
	    int indObj = list_[i];
	    double value, value2;
	    value = computeUsefulness(MAXMIN_CRITERION,
				      upMultiplier, downMultiplier, value,
				      object[indObj], indObj, value2);

#ifdef OLD_USEFULLNESS
	    useful_[i] = -value;
#else
	    if ((sortCrit_ & 1) == 0) {
	      useful_[i] = -value;
	    }
	    else {
	      useful_[i] = value;
	    }
#endif
	  }
	
	  if(sortUpTo > maximumStrong) {
	    // compute from, upto such that 
	    // vPriority[k] == vPriority[maximumStrong] for k in [from..upto-1]
	    int from = maximumStrong-1, upto = maximumStrong;
	    int msPrio = vPriority[maximumStrong-1];
	    problem_->setLastPrioSort(msPrio);
	    while((from > -1) && (vPriority[from] == msPrio)) {
	      from--;
	    }
	    from++;
	    while((upto < sortUpTo) && (vPriority[upto] == msPrio)) {
	      upto++;
	    }
	    // sort list[from]..list[upto-1] according to 
	    // useful_[from]..useful_[upto-1]
	    CoinSort_2(useful_+from, useful_+upto, list_+from);
	  }
	}

#ifdef FM_CHECK
	// priority of last selected object
	double ckPrio = (card_vPriority < numberUnsatisfied_ ?
			 vPriority[card_vPriority] : 100000);
	double ckUse = (card_vPriority < numberUnsatisfied_ ?
		      useful_[card_vPriority] : 100000);
	for(i=0; i<card_vPriority; i++) {
	  int indObj = list_[i];
	  if(object[indObj]->priority() > ckPrio + 1e-3) {
	    printf("CouenneChooseStrong::gutsOfSetupList(): ### ERROR: object[%d]->priority(): %d  > ckPrio: %d\n", 
		   indObj, object[indObj]->priority(), ckPrio);
	    exit(1);
	  }
	  if(fabs(object[indObj]->priority() - ckPrio) < 1e-3) {
	    if(useful_[i] > ckUse + 1e-3) {
	      printf("CouenneChooseStrong::gutsOfSetupList(): ### ERROR: object[%d]->useful: %f  > ckUse: %d\n", 
		   indObj, useful_[i], ckUse);
	      exit(1);
	    }
	  }
	}
	for(i=card_vPriority; i<numberUnsatisfied_; i++) {
	  int indObj = list_[i];
	  if(object[indObj]->priority() < ckPrio - 1e-3) {
	    printf("CouenneChooseStrong::gutsOfSetupList(): ### ERROR: object[%d]->priority(): %d  < ckPrio: %d\n", 
		   indObj, object[indObj]->priority(), ckPrio);
	    exit(1);
	  }
	  if(fabs(object[indObj]->priority() - ckPrio) < 1e-3) {
	    if(useful_[i] < ckUse - 1e-3) {
	      printf("CouenneChooseStrong::gutsOfSetupList(): ### ERROR: object[%d]->useful: %f  < ckUse: %d\n", 
		   indObj, useful_[i], ckUse);
	      exit(1);
	    }
	  }
	}
#endif
	

      }
      else {
	numberUnsatisfied_ = -1;
      }
#else /* not FM_SORT_STRONG */
      // Get list
      numberOnList_=0;
      if (feasible) {
	maximumStrong = CoinMin(maximumStrong,putOther);
	for (i=0;i<maximumStrong;i++) {
	  if (list_[i]>=0) {
#ifdef OLD_USEFULLNESS
	    list_[numberOnList_]=list_[i];
	    useful_[numberOnList_++]=-useful_[i];
	    
#else
	    list_[numberOnList_]=list_[i];
	    if ((sortCrit_ & 1) == 0) {
	      useful_[numberOnList_++]=-useful_[i];
	    }
	    else useful_[numberOnList_++] = useful_[i];
#endif
	    message(CANDIDATE_LIST2)<<numberOnList_-1
				    <<list_[numberOnList_-1]<<numberOnList_-1<<useful_[numberOnList_-1]
				    <<CoinMessageEol;
	  }
	}
	if (numberOnList_) {
	  int tmp_on_list = 0;
	  if (max_most_fra > 0 && numberOnList_ >= maximumStrong) {
	    // If we want to force non-trusted in the list, give them huge
	    // weight here
	    number_not_trusted_=0;
	    for (i=0;i<max_most_fra;i++) {
	      if (list2[i]>=0) {
		list2[number_not_trusted_] = list2[i];
		useful2[number_not_trusted_++] = useful2[i];
		message(CANDIDATE_LIST3)<<number_not_trusted_-1
					<<list2[number_not_trusted_-1]<<number_not_trusted_-1
					<<useful2[number_not_trusted_-1]<<CoinMessageEol;
	      }
	    }
	    if (number_not_trusted_) {
	      CoinSort_2(list_,list_+numberOnList_,useful_);
	      CoinSort_2(list2,list2+number_not_trusted_,useful2);
	      int i1=0;
	      int i2=0;
	      for (i=0; i<numberObjects; i++) {
		bool found1 = (list_[i1]==i);
		bool found2 = (list2[i2]==i);
		if (found1 && found2) {
		  useful_[i1] = -1e150*(1.+useful2[i2]);
		  list2[i2] = -1;
		}
		if (found1) i1++;
		if (found2) i2++;
		if (i2==max_most_fra) break;
	      }
	      for (i=0; i<number_not_trusted_; i++) {
		if (list2[i] >= 0) {
		  list_[numberOnList_+tmp_on_list] = list2[i];
		  useful_[numberOnList_+tmp_on_list] = -1e150*(1.+useful2[i]);
		  tmp_on_list++;
		}
	      }
	    }
	  }
	  // Sort
	  CoinSort_2(useful_,useful_+numberOnList_+tmp_on_list,list_);
	  // move others
	  i = numberOnList_;
	  for (;putOther<numberObjects;putOther++)
	    list_[i++]=list_[putOther];
	  assert (i==numberUnsatisfied_);
	  if (!CoinMax(numberStrong_,numberStrongRoot_))
	    numberOnList_=0;
	}
      }
      else {
	// not feasible
	numberUnsatisfied_=-1;
      }
#endif /* not FM_SORT_STRONG */

      if(!firstPass) {
	break;
      }
      firstPass = false;
    } /* while(numberOnList_ == 0) */

#ifdef TRACE_STRONG
      printf("numberStrong_: %d   maximumStrong: %d\n", 
	     numberStrong_, maximumStrong);
#endif

#ifdef FM_SORT_STRONG
      delete [] vPriority;
#else  /* not FM_SORT_STRONG */
    delete [] list2;
    delete [] useful2;
#endif /* not FM_SORT_STRONG */

    // Get rid of any shadow prices info
    info->defaultDual_ = -1.0; // switch off
    delete [] info->usefulRegion_;
    delete [] info->indexRegion_;

    int way;
    if (bb_log_level_>3) {
      //for (int i=0; i<Min(numberUnsatisfied_,numberStrong_); i++)
      for (int i=0; i<numberOnList_; i++){
        message(CANDIDATE_LIST)<<i<< list_[i]<< i<< useful_[i]
        <<object[list_[i]]->infeasibility(info,way)
        <<CoinMessageEol;
      }
    }
    // return -1 if infeasible to differentiate with numberOnList_==0
    // when feasible
    if(numberUnsatisfied_ == -1) {
      return(-1);
    }
    return numberOnList_;
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

    double obj = solution [problem_ -> Obj (0) -> Body () -> Index ()];

#ifdef FM_CHECKNLP2
    bool isFeas = problem_->checkNLP2(solution, 0, false, true, true, 
				      problem_->getFeasTol());
    return isFeas;
#else
    return problem_ -> checkNLP (solution, obj);
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
      CouenneObject *co =  dynamic_cast <CouenneObject *>(object[i]);
      int indVar = -1;
      if(co) {
	indVar = co->Reference()->Index();
      }
      else {
	indVar = object[i]->columnNumber();
      }
      int way;
      double value = object[i]->infeasibility(info,way);
      maxViol = (value > maxViol ? value : maxViol);
      if(value > 0.0) {
	printf("(%d %d %f)", i, indVar, value);
	minPosViol = (value < minPosViol ? value : minPosViol);
      }
    }
    printf("\nmaxViol: %g  minPosViol: %g\n", maxViol, minPosViol);

  } /* printObjViol */


//}
