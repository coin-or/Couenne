/* $Id$
 *
 * Name:    StrongBranchingSetupList.cpp
 * Authors: Andreas Waechter
 *          Pietro Belotti
 *          Francois Margot, Carnegie Mellon University
 * Purpose: Guts of setup list (including orbital branching)
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneObject.hpp"
#include "BonChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
//#include "CouenneProblemElem.hpp"
//#include "CouenneBranchingObject.hpp"

//#include "CouenneRecordBestSol.hpp"

//The recommended ones:
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

// service types for orbital branching

struct objStrongPri {

  int objIndex_;

  int priority_;
  double value_;
};

inline bool compStrongPri (struct objStrongPri *one, struct objStrongPri *two)  {

  return (one   -> priority_  <  two -> priority_ || 
	  ((one -> priority_  == two -> priority_) &&
	   (one -> value_     >  two -> value_)));
}

/****************************************************************************/
// Copied from BonChooseVariable.cpp and modified slightly
//
// If FM_SORT_STRONG is used:
// Select unsatisfied objects first on priority, then usefulness. 
//
//   If USE_NOT_TRUSTED is also defined, modify usefulness of a fraction
//   of unsatisfied objects with minimum priority (most fractional first)
//   to make them more attractive.  Number of such objects is
//   number_not_trusted_ on return.

// Additional options working with FM_SORT_STRONG (at most one of the two):
//   if FM_ALWAYS_SORT is also defined exact sorting is done. Otherwise, 
//   on return all objects in list[0..numberOnList_] have either smaller 
//   priority or equal priority and usefulness not worse than other 
//   unsatisfied objects.
//
//   if FM_SEC_SORT_USEFUL is defined, objects are selected by priority
//   then usefulness and the final list is sorted according to usefulness.
//
//
// If FM_SORT_STRONG is not used:
// Select unsatisfied objects first on priority, then usefulness
// but in a weird way: Only guarantee is that objects with minimum
// priority will be first in the list, objects with higher priority
// appearing after them in no predictable order. Then
// objects (with priority matching the smallest priority value among 
// unsatisfied objects, most fractional first) have their usefulness 
// modified to make them more attractive. Number of such objects is
// number_not_trusted_ on return. List is then sorted according to usefulness. 
//
// Recommended settings: one of 
// i) define FM_SORT_STRONG USE_NOT_TRUSTED FM_SEC_SORT_USEFUL
// ii) no flags defined (default)

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
    int numberObjects = solver_->numberObjects(); // FIXME: why not info -> solver_ ?
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

    int bestPriority = COIN_INT_MAX;
    //int i;

#ifdef FM_SORT_STRONG
    int numStr = numberStrong_;
    if(isRootNode(info)) {
      numStr = numberStrongRoot_;
    }
    int maximumStrong = CoinMin(numStr, numberObjects) ;
    int lastPrio = problem_->getLastPrioSort();
    int card_vPriority = 0;
    int posEnd_vPriority = numberObjects;
    double *vPriority = new double [numberObjects];
    double *infeasVal = new double [numberObjects];
    int max_most_fra = setup_pseudo_frac_ > 0. ? (int)floor(setup_pseudo_frac_*(double)maximumStrong): 0;
    if (setup_pseudo_frac_ > 0.) {
      max_most_fra = CoinMax(1, max_most_fra);
    }
#else /* not FM_SORT_STRONG */

    int putOther = numberObjects;
    double check = -COIN_DBL_MAX;
    int checkIndex = 0;
    int maximumStrong = CoinMin(CoinMax(numberStrong_,numberStrongRoot_),
        numberObjects) ;
    for (int i=0;i<numberObjects;i++) {
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
      for (int i=0;i<max_most_fra;i++) {
        list2[i]=-1;
        useful2[i]=0.0;
      }
    }
#endif /* not FM_SORT_STRONG */

#ifdef FM_CHECK
    const double* upTotalChange = pseudoCosts_.upTotalChange();
    const double* downTotalChange = pseudoCosts_.downTotalChange();
    int pseudoNum = pseudoCosts_.numberObjects();
    for(int i=0; i<pseudoNum; i++) {
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

    double upMultiplier, downMultiplier;
    computeMultipliers (upMultiplier, downMultiplier);

    // Say feasible
    bool feasible = true;
    const double MAXMIN_CRITERION = maxminCrit(info);

    OsiObject **objectOrig = info -> solver_ -> objects ();

    std::vector <int> objectInd;

    numberObjects = info -> solver_ -> numberObjects ();

    //
    // Strong Orbital branching
    //

    ///////////////////////////////////////////////////////////////////
    //
    // Possibly reduce the set of objects to take into account orbital
    // branching.
    //
    // * Identify objects related to a variable instead of e.g. SOS
    // * Identify orbits
    // * For each orbit:
    //     Keep only one object referring to that orbit
    //
    // No more objects than #orbits should be included in order to
    // save SB time. Solving strong branching LPs for two variables in
    // the same orbit is a waste of time given that, in theory, the
    // branching rules on these two variables are equivalent.

#ifdef COIN_HAS_NTY

    bool useOrbitalBranching = problem_ -> orbitalBranching ();

    if (useOrbitalBranching) {

      int n = problem_ -> nVars ();

      problem_ -> ChangeBounds (info -> lower_, info -> upper_, n);
      problem_ -> Compute_Symmetry ();

      //problem_ -> Print_Orbits ();

      std::vector<std::vector<int> > *orbits = problem_ -> getNtyInfo () -> getOrbits ();

      // bail out if there are only trivial (size 1) orbits

      bool nonTrivialOrbits = false;

      for (std::vector <std::vector<int> >::iterator i = orbits -> begin (); 
	   i != orbits -> end (); ++i)

	if (i -> size () >= 2) {
	  nonTrivialOrbits = true;
	  break;
	}

      //
      // Only do this if there is at least one nontrivial orbit
      //

      if (nonTrivialOrbits) { 

	//printf ("-----------------------------------------------------\n");

	// build map variable -> object (i.e. the inverse of object
	// [i] -> columnNumber ())

	int *varObj = new int [n];

	CoinFillN (varObj, n, -1);

	for (unsigned int i = 0; i < numberObjects; i++) {

	  int indVar = objectOrig [i] -> columnNumber ();

	  if ((indVar >= 0) && 
	      (indVar <  n))
	    varObj [indVar] = i;

	  //if ((indVar >= 0) && (indVar <  n)) printf ("varObj [%d] = %d\n", indVar, i);
	}

	// Objects in orbits aren't just variables, but also
	// constants.  Variable indices are therefore mixed with
	// others, we need to consider only variables.

	for (std::vector <std::vector<int> >::iterator i = orbits -> begin (); 
	     i != orbits -> end (); ++i) {

	  int orbSize = i -> size ();

	  if (orbSize <= 0)
	    continue;

	  if (orbSize == 1) {

	    int orbVar = (*i) [0];

	    if ((orbVar <  n) && 
		(orbVar >= 0)) { // single-variable orbit, not much to do

	      if ((varObj [orbVar] >= 0) &&
		  (objectOrig [varObj [orbVar]] -> checkInfeasibility (info) > COUENNE_EPS))
		objectInd.push_back (varObj [orbVar]); // variable associated with object

	      //if ((varObj [orbVar] >= 0) && (objectOrig [varObj [orbVar]] -> checkInfeasibility (info) > COUENNE_EPS))
	      //printf ("added trivial orbit object %d [x%d]\n", varObj [orbVar], orbVar);

	      continue;
	    }
	  }

	  // Orbit is nontrivial. Sort its variables' objects
	  // according to (non-decreasing) priority.

	  std::vector <struct objStrongPri *> orbitObj; // pre-size vector

	  int
	    minPri =  COIN_INT_MAX,
	    maxPri = -COIN_INT_MAX;

	  double 
	    minValue =  COIN_DBL_MAX,
	    maxValue = -COIN_DBL_MAX;

	  for (std::vector<int>::iterator j = i -> begin (); 
	       j != i -> end (); ++j) 

	    if ((*j < n) && (*j >= 0)) {

	      int objInd = varObj [*j];

	      if (objInd < 0) 
		continue;

	      int pri = objectOrig [objInd] -> priority ();

	      if (pri < minPri) minPri = pri;
	      if (pri > maxPri) maxPri = pri;
		
	      double 
		newValue,
		infeas = objectOrig [objInd] -> checkInfeasibility (info),
		//infeas = objectOrig [(*j) -> objIndex_] -> infeasibility (info, way),
		value  = computeUsefulness (MAXMIN_CRITERION,
					    upMultiplier, downMultiplier, infeas,
					    objectOrig [objInd], objInd, 
					    newValue); // output parameter (ignored)

	      if (value < minValue) minValue = infeas;
	      if (value > maxValue) maxValue = infeas;

	      //printf ("%d(x%d,%d,%g(%g)) ", objInd, *j, pri, infeas, value); fflush (stdout);

	      struct objStrongPri *obj = new struct objStrongPri;

	      obj -> objIndex_ = objInd;
	      obj -> priority_ = pri;
	      obj -> value_    = infeas;

	      orbitObj. push_back (obj);
	    }

	  // If minPri < maxPri (unlikely in orbits), sort objects so
	  // scan below is straightforward

	  if (minPri < maxPri)
	    std::sort (orbitObj.begin (), orbitObj.end (), compStrongPri);

	  // Scan objects in non-decreasing order of priority. Store
	  // in objects (the array of real objects) the one with null
	  // infeasibility and minimum priority.

	  for (std::vector <struct objStrongPri *>::iterator j = orbitObj. begin ();
	       j != orbitObj. end (); ++j) {

	    //printf ("checking object (%d,%d) --> %e\n", (*j) -> objIndex_, (*j) -> priority_, (*j) -> value_);

	    if ((*j) -> value_ > COUENNE_EPS) {

	      // Found object of this orbit that is good enough (and
	      // has highest priority)

	      objectInd.push_back ((*j) -> objIndex_);
	      break;
	    }
	  }

	  for (std::vector <struct objStrongPri *>::iterator j = orbitObj. begin ();
	       j != orbitObj. end (); ++j)
	    delete (*j);
	}

	numberObjects = objectInd.size ();

	// printf ("there are in total %d objects\n", numberObjects);

	// for (int i=0; i<numberObjects; i++) {

	//   printf ("obj %d [%d]: ", i, objectInd [i]); fflush (stdout);

	//   OsiObject *obj = info -> solver_ -> objects () [objectInd [i]];
	//   //int way;

	//   printf ("x%d [%g,%g] %g\n",
	// 	  obj -> columnNumber (),
	// 	  info -> lower_ [obj -> columnNumber ()],
	// 	  info -> upper_ [obj -> columnNumber ()],
	// 	  obj -> checkInfeasibility (info));
	// 	  //obj -> infeasibility (info,way));
	// }

	delete [] varObj;
      }
    }
#endif

    bool firstPass = false; // not important; useful for making two
                            // passes, picking different objects

    OsiObject **object = info -> solver_ -> objects ();

    while (numberOnList_ == 0) { // FIXME: add iteration limit

      for (unsigned int i = 0; i < numberObjects; i++) {

	int indexObj = ((objectInd.size () > 0) && (i < objectInd.size ())) ? objectInd [i] : i;

	int way;
	//double value = object [indexObj] -> checkInfeasibility (info);
	double value = object [indexObj] -> infeasibility (info, way);

	//printf ("object %d[%d]: %g\n", i, indexObj, value);

#ifdef FM_SORT_STRONG
	infeasVal[i] = value;
#endif

	double lbForInfeas = 0.0;
	if(value > lbForInfeas) {
	  numberUnsatisfied_++;
	  if(value >= 1e50) {
	    // infeasible
	    feasible=false;
	    break;
	  }
	  int priorityLevel = object[indexObj]->priority();

#ifdef FM_SORT_STRONG
	  if(priorityLevel < bestPriority) {
	    bestPriority = priorityLevel;	    
	  }
	  if(priorityLevel > lastPrio) {
	    posEnd_vPriority--;
	    vPriority[posEnd_vPriority] = priorityLevel;
	    list_[posEnd_vPriority] = indexObj;
	  }
	  else {
	    vPriority[card_vPriority] = priorityLevel;
	    list_[card_vPriority] = indexObj;
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
	    check = -COIN_DBL_MAX;
	    checkIndex = 0;
	    check2 = -COIN_DBL_MAX;
	    checkIndex2 = 0;
	    number_not_trusted_ = 0;
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
				      object[indexObj], indexObj, value2);
	    if(value > check) {
	      //add to list
	      int iObject = list_[checkIndex];
	      if(iObject >= 0) {
		assert (list_[putOther-1]<0);
		list_[--putOther]=iObject;  // to end
	      }
	      list_[checkIndex]= indexObj;
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
	      list_[--putOther] = indexObj;
	      maximumStrong = CoinMin(maximumStrong,putOther);
	    }

	    if((max_most_fra > 0) && (value2 > check2)) {
	      // add to list of integer infeasibilities
	      number_not_trusted_++;
	      list2[checkIndex2] = indexObj;
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
	    list_[--putOther]= indexObj;
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
	  for(int i=posEnd_vPriority; i<numberObjects; i++) {
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

#ifdef FM_ALWAYS_SORT /* FM_SORT_STRONG */
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
	  for(int i=0; i<card_vPriority; i++) {
	    int indObj = list_[i];
	    double value = 0, value2;
	    value = computeUsefulness(MAXMIN_CRITERION,
				      upMultiplier, downMultiplier, value,
				      object[indObj], indObj, value2);

#ifdef OLD_USEFULLNESS /* FM_SORT_STRONG */
	    useful_[i] = -value; //printf ("useful_ [%d] <-- %g\n", i, -value);
#else
	    if ((sortCrit_ & 1) == 0) {
	      useful_[i] = -value; //printf ("useful_ [%d] <-- %g\n", i, -value);
	    }
	    else {
	      useful_[i] = value; //printf ("useful_ [%d] <-- %g\n", i, value);
	    }
#endif

#ifdef USE_NOT_TRUSTED
	    if(value2 < -COIN_DBL_MAX / 10) { // object with trusted pseudo-cost
	      infeasVal[i] = -COIN_DBL_MAX;
	    }
#endif
	  }

#ifdef USE_NOT_TRUSTED /* FM_SORT_STRONG */
	  // adjust usefulness of objects having priority bestPriority
	  if((card_vPriority > maximumStrong) &&
	     (vPriority[maximumStrong] < bestPriority + COUENNE_EPS)) { 
	    // not all objects with bestPriority will be selected

	    int cardFrac = 0;
	    int *fracInd = new int[card_vPriority]; // holds position 
	                                            // in list_
	    double *fracVal = new double[card_vPriority];

	    // Find all untrusted objects with min priority
	    for(int i=0; i<card_vPriority; i++) {
	      //int indObj = list_[i]; // FIXME: why not used?
	      if(vPriority[i] < bestPriority + COUENNE_EPS) {
		if(infeasVal[i] > -COIN_DBL_MAX/10) {
		  fracInd[cardFrac] = i;
		  fracVal[cardFrac] = -infeasVal[i];
		  cardFrac++;
		}
	      }
	    }

	    if(max_most_fra > 0) {

	    // if more than max_most_fra, sort to pick the ones with max viol
	      if(cardFrac > max_most_fra) {
		CoinSort_2(fracVal, fracVal+cardFrac, fracInd);
	      }
	      for(int i=0; i<cardFrac; i++) {
		useful_[fracInd[i]] = 
		  -1e150*(1. + infeasVal[fracInd[i]]); 
		//-1e150*(1. + infeasVal[list_[fracInd[i]]]);  // FIXME: check if uncommented is correct

		// printf ("useful_ [fracInd[%d]", i);
		// printf (" = %d] <--", fracInd [i]);
		// printf (" %g", useful_ [fracInd [i]]);
		// //printf (" [list[%d]=%d]", fracInd [i], list_ [fracInd [i]]);
		// //printf (" inf=%g\n", infeasVal [list_ [fracInd [i]]]);
		// printf (" [%d]", fracInd [i]);
		// printf (" inf=%g\n", infeasVal [fracInd [i]]);

		number_not_trusted_++;
		if(i == max_most_fra - 1) {
		  break;
		}
	      }
	    }
	    delete[] fracInd;
	    delete[] fracVal;
	  }
#endif /* USE_NOT_TRUSTED and FM_SORT_STRONG */
	}

	//printf ("card_vprio = %d\n", card_vPriority);

#ifdef FM_SEC_SORT_USEFUL
	CoinSort_2(useful_, useful_ + card_vPriority, list_);
#else /* FM_SORT_STRONG not FM_SEC_SORT_USEFUL */
#ifdef FM_ALWAYS_SORT /* FM_SORT_STRONG */
	  int from = 0, upto = 1;
	  while(upto < card_vPriority) {
	    while(vPriority[upto] == vPriority[from]) {
	      upto++;
	      if(upto == card_vPriority) {
		break;
	      }
	    }
	    CoinSort_2(useful_+from, useful_+upto, list_+from);
	    from = upto;
	    upto = from+1;
	  }
#else /* FM_SORT_STRONG not FM_ALWAYS_SORT */
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

#endif /* FM_SORT_STRONG not FM_ALWAYS_SORT */

#ifdef FM_CHECK
	  // priority of last selected object
	  double ckPrio = (card_vPriority < numberUnsatisfied_ ?
			   vPriority[card_vPriority] : 100000);
	  double ckUse = (card_vPriority < numberUnsatisfied_ ?
			  useful_[card_vPriority] : 100000);
	  for(int i=0; i<card_vPriority; i++) {
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
	  for(int i=card_vPriority; i<numberUnsatisfied_; i++) {
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
#endif /* CHECK */
#endif /* FM_SORT_STRONG not FM_ALWAYS_SORT */
      }
      else {
	numberUnsatisfied_ = -1;
      }
#else /* not FM_SORT_STRONG */
      // Get list
      numberOnList_=0;
      if (feasible) {
	maximumStrong = CoinMin(maximumStrong,putOther);
	for (int i=0;i<maximumStrong;i++) {
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
	    for (int i=0;i<max_most_fra;i++) {
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
	      for (int i=0; i<numberObjects; i++) {
		bool found1 = (list_[i1]==i);
		bool found2 = (list2[i2]==i);
		if (found1 && found2) {

#ifdef OLD_USEFULLNESS
		  useful_[i1] = 1e150*(1.+useful2[i2]);
#else
		  useful_[i1] = -1e150*(1.+useful2[i2]);
#endif
		  list2[i2] = -1;
		}
		if (found1) i1++;
		if (found2) i2++;
		if (i2==max_most_fra) break;
	      }
	      for (int i=0; i<number_not_trusted_; i++) {
		if (list2[i] >= 0) {
		  list_[numberOnList_+tmp_on_list] = list2[i];

#ifdef OLD_USEFULLNESS
		  useful_[numberOnList_+tmp_on_list] = 1e150*(1.+useful2[i]);
#else
		  useful_[numberOnList_+tmp_on_list] = -1e150*(1.+useful2[i]);
#endif
		  tmp_on_list++;
		}
	      }
	    }
	  }
	  // Sort
	  CoinSort_2(useful_,useful_+numberOnList_+tmp_on_list,list_);
	  // move others
	  int i = numberOnList_;
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

	//printf ("bailed out of while numberOnList_ == 0\n");
	break;
      }
      firstPass = false;
    } /* while(numberOnList_ == 0) */

#ifdef TRACE_STRONG
    if(problem_->doPrint_) {
      printf("numberStrong_: %d   maximumStrong: %d\n", 
	     numberStrong_, maximumStrong);
    }
#endif

#ifdef FM_SORT_STRONG
      delete [] vPriority;
      delete [] infeasVal;
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
	  //<<object[list_[i]]->checkInfeasibility(info)
        <<CoinMessageEol;
      }
    }

#ifdef COIN_HAS_NTY
    // if (useOrbitalBranching &&
    // 	(objectOrig != object))
    //   delete [] object;
#endif

    // return -1 if infeasible to differentiate with numberOnList_==0
    // when feasible
    if(numberUnsatisfied_ == -1) {
      return(-1);
    }
    return numberOnList_;
  }

