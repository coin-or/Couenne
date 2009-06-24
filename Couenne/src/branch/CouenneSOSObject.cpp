/* $Id: CouenneSOSObject.cpp 141 2009-06-03 04:19:19Z pbelotti $ */
/*
 * Name:    CouenneSOSObject.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Branching SOS object 
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "OsiRowCut.hpp"

#include "CouenneSolverInterface.hpp"
#include "CouenneSOSObject.hpp"
#include "CouenneProblem.hpp"

// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);


/** \brief Execute the actions required to branch, as specified by the
 *	   current state of the branching object, and advance the
 *         object's state.
 *
 *         Returns change in guessed objective on next branch
 */

double CouenneSOSBranchingObject::branch (OsiSolverInterface * solver) {

  // Apply SOS branching first
  double retval = OsiSOSBranchingObject::branch (solver);

  CouenneSolverInterface *couenneSolver = dynamic_cast <CouenneSolverInterface *> (solver);
  CouenneProblem *p = couenneSolver -> CutGen () -> Problem ();

  int 
    nvars  = p -> nVars (),
    objind = p -> Obj (0) -> Body () -> Index ();

  bool infeasible = false;

  p -> domain () -> push (nvars,
			  solver -> getColSolution (), 
			  solver -> getColLower    (), 
			  solver -> getColUpper    ()); // have to alloc+copy

  int       nMembers = dynamic_cast <const OsiSOS *> (originalObject ()) -> numberMembers ();
  const int *Members = dynamic_cast <const OsiSOS *> (originalObject ()) -> members ();

  //CouNumber &estimate = way ? upEstimate_ : downEstimate_;
  CouNumber estimate = 0.;//way ? upEstimate_ : downEstimate_;

  t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

  while (nMembers--) {
    chg_bds [*Members]  .setUpper (t_chg_bounds::CHANGED);
    chg_bds [*Members++].setLower (t_chg_bounds::CHANGED);
  }

  if (     doFBBT_ &&           // this branching object should do FBBT, and
      p -> doFBBT ()) {         // problem allowed to do FBBT

    p -> installCutOff ();

    if (!p -> btCore (chg_bds)) // done FBBT and this branch is infeasible
      infeasible = true;        // ==> report it
    else {

      const double
	*lb = solver -> getColLower (),
	*ub = solver -> getColUpper ();

      //CouNumber newEst = p -> Lb (objind) - lb [objind];
      estimate = CoinMax (0., p -> Lb (objind) - lb [objind]);

      //if (newEst > estimate) 
      //estimate = newEst;

      for (int i=0; i<nvars; i++) {
	if (p -> Lb (i) > lb [i]) solver -> setColLower (i, p -> Lb (i));
	if (p -> Ub (i) < ub [i]) solver -> setColUpper (i, p -> Ub (i));
      }
    }
  }

  /*if (!infeasible && doConvCuts_ && simulate_) { 
    // generate convexification cuts before solving new node's LP

    int nchanged, *changed = NULL;
    OsiCuts cs;

    // sparsify structure with info on changed bounds and get convexification cuts
    sparse2dense (nvars, chg_bds, changed, nchanged);
    couenneSolver -> CutGen () -> genRowCuts (*solver, cs, nchanged, changed, chg_bds);
    free (changed);

    solver -> applyCuts (cs);
    }*/

  delete [] chg_bds;

  p -> domain () -> pop ();

  // next time do other branching
  branchIndex_++;

  // estimated change in objective function
  return (infeasible ? COIN_DBL_MAX : CoinMax (retval, estimate)); 
}


/// create CouenneBranchingObject or CouenneThreeWayBranchObj based
/// on this object
OsiBranchingObject *CouenneSOSObject::createBranch (OsiSolverInterface* si, 
						    const OsiBranchingInformation* info, int way) const
{

  // COPIED FROM OsiBranchingObject.cpp (see code for OsiSOSObject::createBranch)
  int j;
  const double * solution = info->solution_;
  double tolerance = info->primalTolerance_;
  const double * upper = info->upper_;
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    if (upper[iColumn]) {
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (firstNonFixed<0)
	firstNonFixed=j;
      lastNonFixed=j;
      if (value>tolerance) {
	weight += weights_[j]*value;
	if (firstNonZero<0)
	  firstNonZero=j;
	lastNonZero=j;
      }
    }
  }
  assert (lastNonZero-firstNonZero>=sosType_) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  double separator=0.0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  if (sosType_==1) {
    // SOS 1
    separator = 0.5 *(weights_[iWhere]+weights_[iWhere+1]);
  } else {
    // SOS 2
    if (iWhere==lastNonFixed-1)
      iWhere = lastNonFixed-2;
    separator = weights_[iWhere+1];
  }
  // create object

  return new CouenneSOSBranchingObject (si, this, way, separator, jnlst_, doFBBT_, doConvCuts_);
}
