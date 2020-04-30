/* $Id$
 *
 * Name:    CouenneChooseVariable.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiSolverInterface.hpp"

#include "CbcBranchActual.hpp"
#include "CbcModel.hpp"

#include "CouenneChooseVariable.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneObject.hpp"

#ifdef COUENNE_HAS_NAUTY
#include "CouenneNauty.hpp"
#endif

struct objPri {
  int objIndex_;
  int priority_;
};

bool compPri (struct objPri *one, struct objPri *two)  {
  return (one -> priority_ < 
	  two -> priority_);
}

using namespace Couenne;

void eliminateIntegerObjects (OsiSolverInterface *model);
void eliminateIntegerObjects (CbcModel           *model);

/// Default Constructor 
CouenneChooseVariable::CouenneChooseVariable (): 
  OsiChooseVariable (),
  problem_ (NULL) {}


/// Constructor from solver (so we can set up arrays etc)
CouenneChooseVariable::CouenneChooseVariable (const OsiSolverInterface *si,
					      CouenneProblem *p,
					      JnlstPtr jnlst):
  OsiChooseVariable (si),
  problem_ (p),
  jnlst_   (jnlst) {}


/// Copy constructor 
CouenneChooseVariable::CouenneChooseVariable (const CouenneChooseVariable &source):
  OsiChooseVariable (source),
  problem_ (source.problem_),
  jnlst_   (source.jnlst_) {}


/// Assignment operator 
CouenneChooseVariable & CouenneChooseVariable::operator= (const CouenneChooseVariable& rhs) {
  problem_ = rhs.problem_; 
  jnlst_   = rhs.jnlst_;
  return *this;
}


/// Sets up strong list and clears all if initialize is true.
/// Returns number of infeasibilities. 
/// If returns -1 then node is found infeasible
int CouenneChooseVariable::setupList (OsiBranchingInformation *info, bool initialize) {

  static bool firstCall = true;

  int n = problem_ -> nVars ();

  problem_ -> domain () -> push 
    (n,
     info -> solution_, 
     info -> lower_, 
     info -> upper_);

  jnlst_ -> Printf (Ipopt::J_ITERSUMMARY, J_BRANCHING, "----------------- setup list\n");
  if (jnlst_ -> ProduceOutput (Ipopt::J_DETAILED, J_BRANCHING)) {
    printf ("----------------- setup list\n");
    for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++) 
      if (problem_ -> Var (i) -> Multiplicity () > 0) {
	printf ("%4d %20.4g [%20.4g %20.4g]", i, info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
	if (problem_ -> Var (i) -> Type () == AUX) {
	  printf (" expr. %20.4g [%+e] ", (*(problem_ -> Var (i) -> Image ())) (), (*(problem_ -> Var (i) -> Image ())) () - info -> solution_ [i]);
	  problem_ -> Var (i) -> Image () -> print ();
	}
	printf ("\n");
      }
  }

  int retval;

  // Make it stable, in OsiChooseVariable::setupList() numberObjects must be 0.
  //  retval = (solver_ -> numberObjects ()) ? 
  //    OsiChooseVariable::setupList (info, initialize) : 0;

  // Copied OsiChooseVariable::setupList to adjust it to Orbital Branching

  {
    if (initialize) {
      status_=-2;
      delete [] goodSolution_;
      bestObjectIndex_=-1;
      numberStrongDone_=0;
      numberStrongIterations_ = 0;
      numberStrongFixed_ = 0;
      goodSolution_ = NULL;
      goodObjectiveValue_ = COIN_DBL_MAX;
    }

    numberOnList_=0;
    numberUnsatisfied_=0;

    if (firstCall) {

      eliminateIntegerObjects (const_cast <OsiSolverInterface *> (solver_));
      eliminateIntegerObjects (const_cast <OsiSolverInterface *> (info -> solver_));

      firstCall = false;
    }

    int numberObjects = solver_ -> numberObjects();

    assert (numberObjects);

    OsiObject ** object = info -> solver_ -> objects ();

    // false when problem found infeasible
    bool feasible = true;

    // CouenneChooseVariable has numberStrong_ set once, to one, so
    // CouenneChooseVariable::chooseVariable() simply picks the first
    // in the list. There is no point in scanning the whole list of
    // objects, given that infeasibility() and checkInfeasibility()
    // are expensive.
    //
    // Subdivide them by priority, then find the one with minimum
    // priority (primary key) and maximum infeasibility, and store it
    // at list_ [0]. No need for code on orbital branching here, the
    // only change will be during the branching proper.

#define NEW_SETUPLIST

#ifdef NEW_SETUPLIST

    int way;

    std::vector <struct objPri *> listPri;

    for (int i=0; i<numberObjects; i++) {

      struct objPri *singleton = new struct objPri;
      singleton -> objIndex_ = i;
      singleton -> priority_ = object [i] -> priority ();

      listPri.push_back (singleton);
    }

    // for (std::vector <struct objPri *>::iterator i=listPri.begin (); i != listPri.end (); ++i)
    //   printf ("[%d:%d,%e] ", (*i) -> objIndex_, (*i) -> priority_, object [(*i) -> objIndex_] -> infeasibility (info,way));

    // printf (" ---- ");

    std::sort (listPri.begin (), listPri.end (), compPri);

    // for (std::vector <struct objPri *>::iterator i=listPri.begin (); i != listPri.end (); ++i)
    //   printf ("[%d:%d,%e] ", (*i) -> objIndex_, (*i) -> priority_, object [(*i) -> objIndex_] -> infeasibility (info,way));
    // printf ("\n");

    int minPriority = -1;

    double maxInfeas = 0.;

    for (int i=0; i<numberObjects; ++i) {

      int 
	currIndex = listPri [i] -> objIndex_,
	priority  = listPri [i] -> priority_;

      if ((minPriority >= 0) && 
	  (priority > minPriority))
	break;

      //double infeas = object [currIndex] -> infeasibility (info,way); // MOST EXPENSIVE PART
      double infeas = object [currIndex] -> checkInfeasibility (info); // Less expensive

      //printf ("<%d:%d,%e> ", currIndex, priority, infeas);

      if (((minPriority < 0) || (priority == minPriority)) &&
	  (infeas > maxInfeas)) { 

	// printf (" bingo! ");

	if (minPriority < 0)
	  minPriority = priority;

	maxInfeas = infeas;

	++numberUnsatisfied_;

	if (infeas == COIN_DBL_MAX) {

	  feasible = false;
	  break;

	} else {

	  list_   [0] = currIndex;
	  useful_ [0] = infeas;
	}
      }
    }

    // if (feasible) 
    //   if (!numberUnsatisfied_)
    // 	printf ("no violations ");
    //   else
    // 	printf ("Selected %d (%e)", list_ [0], useful_ [0]);
    // else printf ("infeasible ");

    for (std::vector <struct objPri *>::iterator i=listPri.begin (); i != listPri.end (); ++i)
      delete (*i);

#else

    int maximumStrong = numberStrong_ ? CoinMin (numberStrong_, numberObjects) : 1;

    double check = 0.0;

    int
      checkIndex    = 0,
      bestPriority  = COIN_INT_MAX,
      // pretend one strong even if none
      putOther      = numberObjects; // counts downward from end of list

    // init member list_

    for (int i=0; i < maximumStrong; i++) {
      list_   [i] = -1;
      useful_ [i] = 0.;
    }

    // printf ("numObj: %d, maxStrong: %d\n", numberObjects, maximumStrong);
    // for (int i=0; i<numberObjects;) {
    //   printf ("(%d,%d,%d,%g)", i, object [i] -> columnNumber (), object [i] -> priority (), object [i] -> checkInfeasibility (info));
    //   if (!(++i % 10) || i>=numberObjects) printf ("\n");
    // }

    for (int i=0; i<numberObjects; i++) {

      int way;
      //double value = object [i] -> infeasibility (info,way); // TODO: use checkInfeasibility instead
      double value = object [i] -> checkInfeasibility (info); // 

      if (value > 0.) {

	numberUnsatisfied_++;

	if (value == COIN_DBL_MAX) { // infeasible
	  feasible = false;
	  break;
	}

	int priorityLevel = object [i] -> priority ();

	// Better priority? Flush choices
	if (priorityLevel < bestPriority) {

	  for (int j=0; j<maximumStrong; j++) {
	    
	    if (list_ [j] >= 0) {

	      int iObject = list_[j];

	      list_   [j]          = -1;
	      useful_ [j]          = 0.;
	      list_   [--putOther] = iObject;
	    }
	  }

	  bestPriority = priorityLevel;
	  check = 0.;
	} 

	if ((priorityLevel == bestPriority) // only consider those with equal priority
	    && (value > check)

	    ) {

	  int iObject = list_ [checkIndex];
	  if (iObject >= 0)
	    list_ [--putOther] = iObject;  // put former best to end

	  // add to list
	  list_   [checkIndex] = i;
	  useful_ [checkIndex] = value;

	  // Now find next element to (possibly) replace, i.e., find worst
	  check=COIN_DBL_MAX;
	  for (int j=0;j<maximumStrong;j++) {
	    if (list_[j]>=0) {
	      if (useful_[j]<check) {
		check=useful_[j];
		checkIndex=j;
	      }
	    } else {
	      check=0.0;
	      checkIndex = j;
	      break;
	    }
	  }
	} else list_ [--putOther] = i; // greater (weaker) priority or
				       // not most fractional
      }
    }
#endif

    // Get list
    numberOnList_=0;

    if (feasible) {
#ifndef NEW_SETUPLIST
      for (int i=0;i<maximumStrong;i++) {
      	if (list_[i]>=0) {
      	  list_[numberOnList_]=list_[i];
      	  useful_[numberOnList_++]=-useful_[i];
      	}
      }
      if (numberOnList_) {
      	// Sort 
      	CoinSort_2(useful_,useful_+numberOnList_,list_);
      	// move others
      	int i = numberOnList_;
      	for (;putOther<numberObjects;putOther++) 
      	  list_[i++]=list_[putOther];
      	assert (i==numberUnsatisfied_);
      	if (!numberStrong_)
      	  numberOnList_=0;
      } 
#endif
    } else {
      // infeasible
      numberUnsatisfied_ = -1;
    }

    retval = numberUnsatisfied_;
  }

  ////////////////////////////////////////////////////////
  //
  // End copy
  //
  ////////////////////////////////////////////////////////

  problem_ -> domain () -> pop ();

  jnlst_ -> Printf (Ipopt::J_ITERSUMMARY, J_BRANCHING, "----------------- setup list done, %d objects\n", retval);

  // for (int i=0, k = (numberStrong_ ? CoinMin (numberStrong_, solver_ -> numberObjects ()) : 1); i<k; i++) {
  //   printf ("list %3d: %3d ", i, list_ [i]);
  //   if (!((i+1) % 12) || i == k-1) printf ("\n");
  // }

  // printf ("returning %d\n", retval);

  return retval;
}


/// choose object to branch based on earlier setup
// int CouenneChooseVariable::chooseVariable (OsiSolverInterface * solver, 
// 					   OsiBranchingInformation *info, 
// 					   bool fixVariables) {

//   // !!!should go -- just choose the first element
//   problem_ -> domain () -> push 
//     (problem_ -> nVars (),
//      info -> solution_, 
//      info -> lower_, 
//      info -> upper_);

//   int retval = OsiChooseVariable::chooseVariable (solver, info, fixVariables);
//   problem_ -> domain () -> pop ();
//   return retval;
// }


// Returns true if solution looks feasible against given objects
bool CouenneChooseVariable::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects) {

#ifdef FM_CHECKNLP2
  bool isFeas = problem_->checkNLP2(solution,
				    0, 
				    false, // do not care about obj
				    true,  // stopAtFirstViol
				    true,  // checkAll
				    problem_ -> getFeasTol());
  
  return isFeas;
#else
  int indobj = problem_ -> Obj (0) -> Body () -> Index ();
  double obj = indobj >= 0 ? solution [indobj] : problem_ -> Obj (0) -> Body () -> Value ();
  return problem_ -> checkNLP (solution, obj);
#endif
}


/// Add list of options to be read from file
void CouenneChooseVariable::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2 
    ("enable_sos",
     "Use Special Ordered Sets (SOS) as indicated in the MINLP model",
     "no",
     "no","",
     "yes","");

  roptions -> AddStringOption2
    ("branch_fbbt",
     "Apply bound tightening before branching",
     "yes",
     "no","",
     "yes","",
     "After applying a branching rule and before re-solving the subproblem, apply Bound Tightening.");

  roptions -> AddStringOption2
    ("branch_conv_cuts",
     "Apply convexification cuts before branching (for now only within strong branching)",
     "yes",
     "no","",
     "yes","",
     "After applying a branching rule and before resolving the subproblem, generate a round of linearization cuts with the new bounds enforced by the rule."
    );

  roptions -> AddStringOption6
    ("branch_pt_select",
     "Chooses branching point selection strategy",
     "mid-point",
     "lp-clamped", "LP point clamped in [k,1-k] of the bound intervals (k defined by lp_clamp)",
     "lp-central", "LP point if within [k,1-k] of the bound intervals, middle point otherwise" 
     "(k defined by branch_lp_clamp)",
     "balanced", "minimizes max distance from curve to convexification",
     "min-area", "minimizes total area of the two convexifications",
     "mid-point", "convex combination of current point and mid point",
     "no-branch", "do not branch, return null infeasibility; for testing purposes only",
     "");

  std::string br_ops [] = {"prod", "div", "exp", "log", "trig", 
			   "pow",  "negpow", "sqr", "cube", ""};

  for (int i=0; br_ops [i] != ""; i++) {

    char optname [40], optname2 [40], description [90];
    sprintf (optname,  "branch_pt_select_%s", br_ops [i].c_str ());
    sprintf (optname2, "branch_lp_clamp_%s",  br_ops [i].c_str ());
    sprintf (description, "Chooses branching point selection strategy for operator %s.",
	     br_ops [i].c_str ());

    roptions -> AddStringOption7
      (optname,
       description,
       "common",
       "common",    "use strategy defined for generic operators",
       "lp-clamped", "LP point clamped in [k,1-k] of the bound intervals "
       "(k defined by lp_clamp_${this operator}$)",
       "lp-central", "LP point if within [k,1-k] of the bound intervals, middle point otherwise" 
       "(k defined by branch_lp_clamp_${this operator}$)",
       "balanced",  "minimizes max distance from curve to convexification",
       "min-area",  "minimizes total area of the two convexifications",
       "mid-point", "convex combination of current point and mid point",
       "no-branch", "do not branch, return null infeasibility; for testing purposes only",
       "");

    roptions -> AddBoundedNumberOption
      (optname2,
       "Defines safe interval percentage [0,0.5] for using LP point as a branching point.",
       0.,false,
       0.5,false,
       0.2);
  }

  roptions -> AddBoundedNumberOption
    ("branch_midpoint_alpha",
     "Defines convex combination of mid point and current LP point: "
     "b = alpha x_lp + (1-alpha) (lb+ub)/2.",
     0.,false,
     1.,false,
     default_alpha);

  roptions -> AddBoundedNumberOption
    ("branch_lp_clamp",
     "Defines safe interval percentage for using LP point as a branching point.",
     0.,false,
     1.,false,
     0.2);

  // Setting priorities slightly below CbcBranchingObjects' priority,
  // so that Couenne's integer branching is used

  roptions -> AddLowerBoundedIntegerOption
    ("cont_var_priority",
     "Priority of continuous variable branching",
     1, 99,
     "When branching, this is compared to the priority of integer variables, whose priority is given by int_var_priority, and SOS, whose priority is 10. "
     "Higher values mean smaller priority."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("int_var_priority",
     "Priority of integer variable branching",
     1, 98,
     "When branching, this is compared to the priority of continuous variables, whose priority is given by cont_var_priority, and SOS, whose priority is 10. "
     "Higher values mean smaller priority."
    );

  roptions -> AddStringOption2
    ("red_cost_branching",
     "Apply Reduced Cost Branching (instead of the Violation Transfer) -- MUST have vt_obj enabled",
     "no",
     "no", "Use Violation Transfer with $\\sum |\\pi_i a_{ij}|$",
     "yes","Use Reduced cost branching with $|\\sum \\pi_i a_{ij}|$");

  roptions -> AddStringOption2 
    ("orbital_branching",
     "detect symmetries and apply orbital branching",
     "no",
     "yes", "",
     "no", "");

  roptions -> AddLowerBoundedIntegerOption
    ("orbital_branching_depth",
     "Maximum depth at which the symmetry group is computed",
     -1, 10,
     "Select -1 if you want to compute the symmetry group at all nodes");
}


// Cbc adds automatically objects for integer variables, but we did
// not ask for it. In order to overcome this, we add objects the usual
// way and then eliminate the CbcSimpleInteger objects one by one.

int gutsofEIO (OsiObject **objects, int nco) {

  int 
    nRealObj, 
    currObj  = 0;

  for (; currObj < nco; ++currObj) 

    if ((NULL != dynamic_cast <CbcSimpleInteger *> (objects [currObj])) ||
	(NULL != dynamic_cast <OsiSimpleInteger *> (objects [currObj]))) {

      // At [currObj] is a Cbc integer variable object. Kill it! Kill it with fire!
      delete objects [currObj];
      objects [currObj] = NULL;
    }

  // squeeze the sparse vector into a dense one with only non-NULL entries

  for (nRealObj = 0, currObj = -1; nRealObj < nco; ++nRealObj)

    if (NULL == objects [nRealObj]) {

      if (currObj < 0) 
	currObj = nRealObj + 1;

      while ((currObj < nco) && 
	     (NULL == objects [currObj]))
	++currObj;

      if (currObj >= nco)
	break;

      objects [nRealObj] =
      objects [currObj];

      objects [currObj] = NULL;
    }

  //printf ("%d real objects out of %d (s.co %d)\n", nRealObj, currObj, s.continuousSolver () -> numberObjects ());

  return nRealObj;
}

void eliminateIntegerObjects (OsiSolverInterface *model) {model -> setNumberObjects (gutsofEIO (model -> objects (), model -> numberObjects ()));}
void eliminateIntegerObjects (CbcModel           *model) {model -> setNumberObjects (gutsofEIO (model -> objects (), model -> numberObjects ()));}
