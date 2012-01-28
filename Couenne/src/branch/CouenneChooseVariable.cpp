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

#include "CouenneChooseVariable.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneObject.hpp"

#ifdef COIN_HAS_NTY
#include "Nauty.h"
#endif

using namespace Couenne;

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

  int n = problem_ -> nVars ();

  problem_ -> domain () -> push 
    (n,
     info -> solution_, 
     info -> lower_, 
     info -> upper_);

#ifdef COIN_HAS_NTY

  bool useOrbitalBranching = problem_ -> orbitalBranching ();

  int  *whichOrbit    = NULL;
  bool *orbitIncluded = NULL; // true if one variable from orbit already candidate

  if (useOrbitalBranching) {

    problem_ -> ChangeBounds (info -> lower_,
			      info -> upper_, n);
    
    problem_ -> Compute_Symmetry ();

    whichOrbit    = new int  [n];
    orbitIncluded = new bool [n]; // true if one variable from orbit already candidate

    CoinFillN (whichOrbit,    n, -1);
    CoinFillN (orbitIncluded, n, false);

    std::vector<std::vector<int> > *orbits = problem_ -> getNtyInfo () -> getOrbits ();

    int orbIndex = 0;

    // Objects in orbits aren't just variables, but also constants.
    // Variable indices are therefore mixed with others, we need to
    // consider only variables.

    for   (std::vector <std::vector<int> >::iterator i = orbits -> begin (); i != orbits -> end (); ++i, ++orbIndex)
      for (             std::vector<int>  ::iterator j = i      -> begin (); j != i      -> end (); ++j)
	if ((*j < n) && (*j >= 0)) // only save it if it is a variable's index
	  whichOrbit [*j] = orbIndex;

    delete orbits;
  }

#endif

  jnlst_ -> Printf (Ipopt::J_ITERSUMMARY, J_BRANCHING, "----------------- setup list\n");

  if (jnlst_ -> ProduceOutput (Ipopt::J_DETAILED, J_BRANCHING)) {

    printf ("----------------- setup list\n");

    for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++) 

      if (problem_ -> Var (i) -> Multiplicity () > 0) {
	printf ("%4d %20.4g [%20.4g %20.4g]", i,
		info -> solution_ [i],
		info -> lower_ [i],
		info -> upper_ [i]);

	if (problem_ -> Var (i) -> Type () == AUX) {
	  printf (" expr. %20.4g [%+e] ", 
		  (*(problem_ -> Var (i) -> Image ())) (), 
		  (*(problem_ -> Var (i) -> Image ())) () - info -> solution_ [i]);
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

    int numberObjects = solver_ -> numberObjects();

    assert (numberObjects);

    double check = 0.0;

    int
      checkIndex    = 0,
      bestPriority  = COIN_INT_MAX,
      // pretend one strong even if none
      maximumStrong = numberStrong_ ? CoinMin (numberStrong_, numberObjects) : 1,
      putOther      = numberObjects; // counts downward from end of list

    // init member list_

    for (int i=0; i < maximumStrong; i++) {
      list_   [i] = -1;
      useful_ [i] = 0.;
    }

    OsiObject ** object = info -> solver_ -> objects ();

    // false when problem found infeasible
    bool feasible = true;

    for (int i=0; i<numberObjects; i++) {

      int way;
      double value = object [i] -> infeasibility (info, way); // TODO: use checkInfeasibility instead

      if (value > 0.) {

	numberUnsatisfied_++;

	if (value == COIN_DBL_MAX) { // infeasible
	  feasible = false;
	  break;
	}

	int priorityLevel = object [i] -> priority ();

	// Better priority? Flush choices
	if (priorityLevel < bestPriority) {

	  // TODO: isn't it "j<putOther"? This overwrites list_
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

#ifdef COIN_HAS_NTY

	    // Add further check if Orbital Branching is used: is
	    // another variable from the same orbit already in the
	    // list? If so, no need to do add object
	    && (!useOrbitalBranching                              || 
		(object [i] -> columnNumber () < 0)              ||
		(whichOrbit [object [i] -> columnNumber ()] < 0) || 
		!(orbitIncluded [whichOrbit [object [i] -> columnNumber ()]]))
#endif
	    ) {

	  int iObject = list_ [checkIndex];
	  if (iObject >= 0)
	    list_ [--putOther] = iObject;  // put former best to end

	  // add to list
	  list_   [checkIndex] = i;
	  useful_ [checkIndex] = value;

#ifdef COIN_HAS_NTY
	  if (useOrbitalBranching && 
	      (object [i] -> columnNumber () >= 0) && 
	      (whichOrbit [object [i] -> columnNumber ()] >= 0))
	    orbitIncluded [whichOrbit [object [i] -> columnNumber ()]] = true;
#endif

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

    // Get list
    numberOnList_=0;

    if (feasible) {
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

  jnlst_ -> Printf (Ipopt::J_ITERSUMMARY, J_BRANCHING, "----------------- setup list done, %d objects\n", 
		    retval);

#ifdef COIN_HAS_NTY
  if (useOrbitalBranching) {
    delete [] whichOrbit;
    delete [] orbitIncluded;
  }
#endif

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

  double obj = solution [problem_ -> Obj (0) -> Body () -> Index ()];

#ifdef FM_CHECKNLP2
  bool isFeas = problem_->checkNLP2(solution,
				    0, false, // do not care about obj
				    true, // stopAtFirstViol
				    true, // checkAll
				    problem_->getFeasTol());
  
  return isFeas;
#else
  return problem_ -> checkNLP(solution, obj);
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
    sprintf (description, "Chooses branching point selection strategy for operator %s", 
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
       "Defines safe interval percentage [0,0.5] for using LP point as a branching point",
       0.,false,
       0.5,false,
       0.2,
       "Default value is 0.2.");
  }

  roptions -> AddBoundedNumberOption
    ("branch_midpoint_alpha",
     "Defines convex combination of mid point and current LP point: "
     "b = alpha x_lp + (1-alpha) (lb+ub)/2.",
     0.,false,
     1.,false,
     default_alpha,
     "Default value is 0.25.");

  roptions -> AddBoundedNumberOption
    ("branch_lp_clamp",
     "Defines safe interval percentage for using LP point as a branching point",
     0.,false,
     1.,false,
     0.2,
     "Default value is 0.2.");

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

    roptions -> AddStringOption2 (
      "orbital_branching",
      "detect symmetries and apply orbital branching",
      "no",
      "yes", "",
      "no", "");
}
