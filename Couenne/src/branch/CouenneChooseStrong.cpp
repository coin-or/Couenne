/* $Id$
 *
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 *          Pietro Belotti, Lehigh University
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "BonChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneBranchingObject.hpp"

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
         
	  results_.push_back(Bonmin::HotInfo(solver, info,
                                solver->objects(), iObject));
        }
        else {
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
          for (unsigned int i=0;i < results_.size();i++) {
            int iObject = results_[i].whichObject();
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
                const OsiObject * obj = solver->object(iObject);
                OsiBranchingObject * branch = obj->createBranch(solver,info,0);
                branch->branch(solver);
                delete branch;
              }
            }
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
                const OsiObject * obj = solver->object(iObject);
                OsiBranchingObject * branch = obj->createBranch(solver,info,1);
                branch->branch(solver);
                delete branch;
              }
            }

            double
	      MAXMIN_CRITERION = maxminCrit(info),	      
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

      if ( bestObjectIndex_ >=0 ) {
        OsiObject * obj = solver->objects()[bestObjectIndex_];
        obj->setWhichWay(bestWhichWay_);
        message(BRANCH_VAR)<<obj->columnNumber()<< bestWhichWay_ <<CoinMessageEol;
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

    problem_ -> domain () -> pop ();

    return retval;
  }


  // Sets up strong list and clears all if initialize is true.
  // Returns number of infeasibilities.
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {

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

    // call Bonmin's setuplist
    int retval = Bonmin::BonChooseVariable::setupList (info, initialize);

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "----------------- (strong) setup list done - %d infeasibilities\n", retval);

    problem_ -> domain () -> pop ();
    return retval;
  }


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
  }


  // Returns true if solution looks feasible against given objects
  bool CouenneChooseStrong::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects) {

    double obj = solution [problem_ -> Obj (0) -> Body () -> Index ()];
    return problem_ -> checkNLP (solution, obj);
  }
//}
