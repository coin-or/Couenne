/* $Id$
 *
 * Name:    doStrongBranching.cpp
 * Authors: Andreas Waechter, IBM Corp.
 *          Pietro Belotti, CMU
 * Purpose: actual strong branching method
 *
 * (C) Carnegie-Mellon University, 2008-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinTime.hpp"
#include "BonChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

//#define TRACE_STRONG
//#define TRACE_STRONG2

using namespace Ipopt;
using namespace Couenne;

/// Called from simulateBranch and from disjunctive cut generators
/// when object is not CouenneObject and therefore needs explicit FBBT
bool BranchingFBBT (CouenneProblem *problem,
		    OsiObject *Object,
		    OsiSolverInterface *solver);

/// compute Euclidean distance between two points (most likely LP solutions)
/// l_2 norm by default, but can change it by fourth parameter
double distance (const double *p1, const double *p2, int size, double k=2.) {

  double 
    result = 0.,
    element;

  if (k == 2.) // a bit faster, probably

    while (size--) {
      element = *p1++ - *p2++;
      result += element * element;
    }

  else

    while (size--) {
      element = *p1++ - *p2++;
      result += pow (element, k);
    }

  return pow (result, 1. / k);
}


/**  This is a utility function which does strong branching on
     a list of objects and stores the results in OsiHotInfo.objects.
     On entry the object sequence is stored in the OsiHotInfo object
     and maybe more.
     It returns -

     -1 - one branch was infeasible both ways
     0 - all inspected - nothing can be fixed
     1 - all inspected - some can be fixed (returnCriterion==0)
     2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
     3 - returning because max time
*/
  int CouenneChooseStrong::doStrongBranching (OsiSolverInterface *solver, 
					      OsiBranchingInformation *info,
					      int numberToDo, int returnCriterion) {

#ifdef TRACE_STRONG2
    int pv = -1;
    if(problem_->doPrint_) {
      if(pv > -1) {
	printf("doSB: beg: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("doSB: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("doSB: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }
    }
#endif

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "\n-\n------- CCS: trying %d objects:\n", numberToDo);

    //solver -> doingResolve () = false; // turns off setCutoff and restoreUnused

    int numberColumns = solver -> getNumCols ();

    solver -> markHotStart (); // save current LP point

    // save initial bounds
    const double
      *initLower = info -> lower_,
      *initUpper = info -> upper_;

    // save intersection of bounds obtained by branching on all objects
    double 
      *saveLower = CoinCopyOfArray (info -> lower_, numberColumns),
      *saveUpper = CoinCopyOfArray (info -> upper_, numberColumns),

      *Lower0 = NULL,
      *Upper0 = NULL,

      // save union of bounds on both branches of one object;
      // reset to (current) saveLower when branching on a new object
      *unionLower  = new double [numberColumns],
      *unionUpper  = new double [numberColumns],

      *lpSol     = NULL, 
       timeStart = CoinCpuTime ();

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      Lower0 = CoinCopyOfArray (info -> lower_, numberColumns); // delete afterwards
      Upper0 = CoinCopyOfArray (info -> upper_, numberColumns);
    }

    // LP solution for distance
    if (pseudoUpdateLP_) 
      lpSol = CoinCopyOfArray (info -> solution_, numberColumns);

    // provide Couenne problem with point/bounds contained in info
    // problem_ -> domain () -> push
    //   (problem_ -> nVars (),
    //    info -> solution_,
    //    info -> lower_,
    //    info -> upper_);

    //Bonmin::HotInfo * results = results_ ();

    int returnCode = 0, iDo = 0;

    for (iDo = 0; iDo < numberToDo; iDo++) {

      Bonmin::HotInfo * result = results_ () + iDo; // retrieve i-th object to test

      OsiObject *Object = solver_ -> objects () [result -> whichObject ()];

      // TODO: apply isCuttable()     

      // TODO: set a cutoff for dual bound in dual simplex
      //       do the same for primal based on SB's alpha

      // For now just 2 way
      OsiBranchingObject * branch = result -> branchingObject ();
      assert (branch->numberBranches()==2);

      CouenneBranchingObject *cb = dynamic_cast <CouenneBranchingObject *> (branch);

      if (cb) cb -> setSimulate (true);

      int 
	status0 = -1, 
	status1 = -1;

      ///////////////////////////////////////////////////////////////////////////

      /* Try the first direction.  Each subsequent call to branch()
	 performs the specified branch and advances the branch object
	 state to the next branch alternative. */

      bool isInf0 = false;
      bool isInf1 = false;

      double indUb = 0, indLb = 0;
      CouenneObject *CouObj = dynamic_cast <CouenneObject *> (Object);
      OsiSimpleInteger *simpl = dynamic_cast <OsiSimpleInteger *>(solver_->objects()[result->whichObject ()]);
      // if OsiSimpleInteger Object with branching point outside 
      // current solver bound interval, one branch must
      // be set as infeasible, otherwise bounds are enlarged 
      // in one branch
      int objVarIndex = -1;
      if(CouObj) {
	objVarIndex = CouObj->Reference()->Index();
      }
      else {
	objVarIndex = Object->columnNumber();
      }
      if(simpl) { 
	if(objVarIndex >= 0) {
	  indUb = solver->getColUpper()[objVarIndex];
	  indLb = solver->getColLower()[objVarIndex];
	  if(info->solution_[objVarIndex] < indLb) {
	    isInf0 = true;
	  }
	  if(info->solution_[objVarIndex] > indUb) {
	    isInf1 = true;
	  }
	}
      }

      status0 = simulateBranch (Object, info, branch, solver, result, -1);
      if(isInf0) {
	status0 = 1; // branch was known to be infeasible
	result->setDownStatus(1);
      }

      // save current bounds as tightened by the down branch; will be           
      // used below to update global bounding box in solver                     
      // if status0 == 1 unionLower will be ignored below                         
      CoinCopyN (solver->getColLower(), numberColumns, unionLower);
      CoinCopyN (solver->getColUpper(), numberColumns, unionUpper);

      // Restore pre-left-branch bounds in solver and problem
      for (int j=0; j<numberColumns; j++) {
	if(problem_->Lb(j) > unionLower[j]) {
	  unionLower[j] = problem_->Lb(j);
	}
	if(problem_->Ub(j) < unionLower[j]) {
	  unionLower[j] = problem_->Ub(j);
	}

        solver->setColLower(j, saveLower [j]);
        solver->setColUpper (j, saveUpper [j]);
	problem_ -> Lb (j) = saveLower [j];
	problem_ -> Ub (j) = saveUpper [j];
      }

      /* second direction */

      status1 = simulateBranch (Object, info, branch, solver, result, +1);
      if(isInf1) {
	status1 = 1; // branch was known to be infeasible
	result->setUpStatus(1);
      }

#ifdef TRACE_STRONG
      if(problem_->doPrint_) {
	printf("Strong on object %d: status0: %d  status1: %d\n", 
	       result->whichObject(), status0, status1);
      }
#endif

      ///////////////////////////////////////////////////////////////////////////

      jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "-------\n");

      if (cb) 
	cb -> setSimulate (false);

      /////////////////////////////////////////////////////////////////////////////

      bool tightened = false;

      t_chg_bounds *chg_bds = new t_chg_bounds [numberColumns];

      const double *sLb = solver->getColLower();
      const double *sUb = solver->getColUpper();

      if(status1 != 1) { // feasible
        if(status0 != 1) { // feasible; take union of both branches
	  for (int j=0; j<numberColumns; j++) {
	    double maxLb = (problem_->Lb(j) < sLb[j] ? sLb[j] : problem_->Lb(j));
	    double minUb = (problem_->Ub(j) < sUb[j] ? problem_->Ub(j) : sUb[j]);
	    problem_->Lb(j) = (unionLower[j] > maxLb ? maxLb : unionLower [j]);
	    problem_->Ub(j) = (unionUpper[j] < minUb ? minUb : unionUpper [j]);
	  }
	}
	else { // keep current bounds; best of problem_ and solver
	  for (int j=0; j<numberColumns; j++) {
	    problem_->Lb(j) = (problem_->Lb(j) < sLb[j] ? sLb[j] : problem_->Lb(j));
	    problem_->Ub(j) = (problem_->Ub(j) < sUb[j] ? problem_->Ub(j) : sUb[j]);
	  }
	}
      }
      else { // branch 1 infeasible
	if(status0 != 1) { // feasible; otherwise both branches are infeasible  
                           // keep current inconsistant bounds in solver        
          for (int j=0; j<numberColumns; j++) {                                 
            problem_->Lb (j) = unionLower [j];                                  
            problem_->Ub (j) = unionUpper [j];                                  
          }                                                                     
        }                                                                       
      }                                                                         

      if((status0 == 1) && (status1 == 1)) {
	tightened = false;

	// make sure that bounds in solver proves problem is
	// infeasible
	double lbVar0 = solver->getColLower()[0];
	if(lbVar0 < 1) {
	  solver->setColLower(0, 1);
	  solver->setColUpper(0, 0);
	}
	else {
	  solver->setColUpper(0, lbVar0-1);
	}
      }
      else {
	for (int j=0; j<numberColumns; j++) {
	  if (problem_ -> Lb (j) > initLower [j] + COUENNE_EPS) {
	    chg_bds [j].setLower (t_chg_bounds::CHANGED);
	    tightened = true;
	  }
	  
	  if (problem_ -> Ub (j) < initUpper [j] - COUENNE_EPS) {
	    chg_bds [j].setUpper (t_chg_bounds::CHANGED);
	    tightened = true;
	  }
	}
      }
      if (tightened &&                     // have tighter bounds
	  (problem_ -> doFBBT ()) &&       // selected FBBT
	  !(problem_ -> btCore (chg_bds))) // tighten again on root

	status0 = status1 = 1;	           // if returns false, problem is infeasible

      delete [] chg_bds;


      if((status0 != 1) || (status1 != 1)) {

	// set new bounding box as the possibly tightened one (a subset
	// of the initial)
	for (int j=0; j<numberColumns; j++) {
	  solver -> setColLower (j, saveLower [j] = problem_ -> Lb (j));
	  solver -> setColUpper (j, saveUpper [j] = problem_ -> Ub (j));
	}
      }

      /*
        End of evaluation for this candidate object. Possibilities are:

        * Both sides below cutoff; this variable is a candidate for
          branching.

        * Both sides infeasible or above the objective cutoff: no
          further action here. Break from the evaluation loop and
          assume the node will be purged by the caller.

        * One side feasible and below cutoff: Install the branch
          (i.e., fix the variable). Possibly break from the evaluation
          loop and assume the node will be reoptimised by the caller.
      */

      if (status0 == 1 && 
	  status1 == 1) { // infeasible
        returnCode=-1;
        break; // exit loop
      } else if (status0==1 || status1==1) {
        numberStrongFixed_++;
        if (!returnCriterion) {
	  returnCode=1;
        } else {
	  returnCode=2;
	  break;
        }
      }

      bool hitMaxTime = ( CoinCpuTime()-timeStart > info->timeRemaining_);
      if (hitMaxTime) {
        returnCode=3;
        break;
      }
    } // end loop /***********************************/

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("strong branching: tightened bounds. ");
      // create union of bounding box from both branching directions
      for (int j=0; j<numberColumns; j++) {
      
	if (problem_ -> Lb (j) > Lower0 [j]) printf ("l%d (%g-->%g) ", j,Lower0[j], problem_->Lb (j));
	if (problem_ -> Ub (j) < Upper0 [j]) printf ("u%d (%g-->%g) ", j,Upper0[j], problem_->Ub (j));
      }

      delete [] Lower0;
      delete [] Upper0;
    }

    //problem_ -> domain () -> pop (); // discard current point/bounds from problem

    delete [] lpSol;

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "----------------------done\n\n\n");

#ifdef TRACE_STRONG2
    if(problem_->doPrint_) {
      if(pv > -1) {
	printf("doSB: beg: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, solver->getColSolution()[pv], solver->getColLower()[pv], solver->getColUpper()[pv]);
	printf("doSB: info: x[%d]: %10.4f  lb: %10.4f  ub: %10.4f\n",
	       pv, info->solution_[pv], info->lower_[pv], info->upper_[pv]);
	printf("doSB: problem: lb: %10.4f  ub: %10.4f\n",
	       problem_->Lb(pv), problem_->Ub(pv));
      }
    }
#endif

    if (iDo < numberToDo) iDo++; // exited due to infeasibility
    assert (iDo <= (int) results_.size());
    results_.resize (iDo);

    delete [] unionLower;
    delete [] unionUpper;

    delete [] saveLower;
    delete [] saveUpper;

    solver -> unmarkHotStart ();     // Delete the snapshot

    //solver -> doingResolve () = true;
    branchtime_ += CoinCpuTime () - timeStart;

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "Done doStrongBranching\n");

    return returnCode;
  }

// Do one side of strong branching
int CouenneChooseStrong::simulateBranch (OsiObject *Object,
					 OsiBranchingInformation *info,
					 OsiBranchingObject *branch,
					 OsiSolverInterface *solver,
					 Bonmin::HotInfo * result,
					 int direction) {

  bool boundBranch = branch -> boundBranch ();

  int status = -1;

  OsiSolverInterface *thisSolver = 
    boundBranch ? solver : solver -> clone ();

  CouenneObject *CouObj = dynamic_cast <CouenneObject *> (Object);

  if ((branch -> branch (thisSolver) > COUENNE_INFINITY) || // branch is infeasible
      // Bound tightening if not a CouenneObject -- explicit since
      // FBBT is done at ::branch() for CouenneObjects
      (!CouObj && !BranchingFBBT (problem_, Object, thisSolver))) {

    status = 1;

    if (direction < 0) result -> setDownStatus (1);
    else               result -> setUpStatus   (1);

  } else {

    if (boundBranch) // branching rule is a variable bound, can use hotstart

      thisSolver -> solveFromHotStart ();

    else { // branching rule is more complicated, need a resolve

      int limit;
      thisSolver -> getIntParam (OsiMaxNumIterationHotStart, limit);
      thisSolver -> setIntParam (OsiMaxNumIteration,         limit); 

      thisSolver -> resolve ();
    }

    if (pseudoUpdateLP_ && CouObj && thisSolver -> isProvenOptimal ()) {
      CouNumber dist = distance (info -> solution_, thisSolver -> getColSolution (), 
				 problem_ -> nVars ());
      if (dist > COUENNE_EPS)
	CouObj -> setEstimate (dist, direction < 0 ? 0 : 1);
    }
  }

  // Can check if we got solution
  // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution

  // only update information if this branch is feasible
  if (status < 0)
    status = result -> updateInformation (thisSolver, info, this);

  numberStrongIterations_ += thisSolver -> getIterationCount ();

  if ((status == 3) && (trustStrongForSolution_)) {
    // new solution already saved
    info -> cutoff_ = goodObjectiveValue_;
    //problem_ -> setCutOff (goodObjectiveValue_);
    status = 0;
  }

  if (solver != thisSolver)
    delete thisSolver;

  return status;
}

/// Called from simulateBranch when object is not CouenneObject and
/// therefore needs explicit FBBT
bool BranchingFBBT (CouenneProblem *problem,
		    OsiObject *Object,
		    OsiSolverInterface *solver) {

  bool feasible = true;

  if (problem -> doFBBT ()) {

    problem -> domain () -> push (solver);

    int 
      indVar = Object  -> columnNumber (),
      nvars  = problem -> nVars ();

    // do not perform this if object is not a variable object

    if (indVar >= 0) {

      // Tell Couenne one bound has changed

      t_chg_bounds *chg_bds = new t_chg_bounds [nvars];
      chg_bds [indVar].setUpper (t_chg_bounds::CHANGED);
      chg_bds [indVar].setLower (t_chg_bounds::CHANGED);

      problem -> installCutOff ();

      if ((feasible = problem -> btCore (chg_bds))) {

	const double
	  *lb  = solver -> getColLower (),
	  *ub  = solver -> getColUpper (),
	  *nLB = problem -> Lb (),
	  *nUB = problem -> Ub ();

	for (int i=0; i<nvars; i++) {
	  if (nLB [i] > lb [i]) solver -> setColLower (i, nLB [i]);
	  if (nUB [i] < ub [i]) solver -> setColUpper (i, nUB [i]);
	}
      }

      delete [] chg_bds;
    }

    problem -> domain () -> pop ();
  }

  return feasible;
}
