/* $Id$
 *
 * Name:    TwoImpliedGenCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: bound reduction using two inequalities from the LP relaxation
 * 
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <set>
#include <stdlib.h>

#include "BonCbc.hpp"
#include "BonBabInfos.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglCutGenerator.hpp"
#include "CoinTime.hpp"

#include "CouenneTypes.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneTwoImplied.hpp"
#include "CouenneExprVar.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneInfeasCut.hpp"
#include "CouenneJournalist.hpp"

using namespace Ipopt;

// necessary to make updateBranchInfo visible
namespace Couenne {

/// get new bounds from parents' bounds + branching rules
void updateBranchInfo (const OsiSolverInterface &si, CouenneProblem *p, 
		       t_chg_bounds *chg, const CglTreeInfo &info);

// do single pair of inequalities. Return < 0 if infeasible

int combine (CouenneProblem *p,
	     int n1, int n2, 
	     const int *ind1, // indices
	     const int *ind2, 
	     double *sa1, // coeff (sparse array)
	     double *sa2,
	     const double *a1,  // coeff
	     const double *a2, 
	     double *clb, // variable bounds
	     double *cub,
	     double l1, // constraint bounds
	     double l2,  
	     double u1, 
	     double u2, 
	     bool *isInteger,
	     int sign); // invert second constraint? -1: yes, +1: no

/// the main CglCutGenerator
void CouenneTwoImplied::generateCuts (const OsiSolverInterface &si, 
				      OsiCuts &cs, 
				      const CglTreeInfo info)
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
  const
#endif
  {

  // don't perform this is cs has been added an infeasible cut (a
  // result of some bound tightening procedure discovering an
  // infeasible node)

  if (isWiped (cs))
    return;

  double now = CoinCpuTime ();

  // a more elaborate scheme to avoid heavy use of this heavy procedure

  if ((depthStopSeparate_ >= 0 &&           // if -1, there is no limit on depth
       info.level > depthStopSeparate_)     // otherwise, check if too deep for adding these cuts
      ||
      (depthLevelling_ >= 0 &&              // chance to run this procedure
       info.level >= depthLevelling_ &&
       CoinDrand48 () > 1. / (2. + info.level - depthLevelling_)))
    return;

  if (info.level <= 0)
    jnlst_ -> Printf (J_ERROR, J_COUENNE, "TwoImpl-BT: "); fflush (stdout);

  // printf ("probexecute = %g. Level = %d, depthlevelling = %d, depthStop = %d, cond = %d\n", 
  // 	  1. / (2. + info.level - depthLevelling_),
  // 	  info.level,
  // 	  depthLevelling_, 
  // 	  depthStopSeparate_,
  // 	  (depthLevelling_ < 0 || info.level < depthLevelling_));

  // Update CouenneProblem's bounds using si's getCol{Low,Upp}er() and
  // cs's OsiColCuts
  problem_ -> domain () -> push (&si, &cs);

  static int nBadColMatWarnings = 0;

  std::set <std::pair <int, int> > pairs;

  /// In principle, "si. getMatrixByCol ()" should be sufficient.
  /// However, there seems to be a bug (not sure where... Osi? Clp?
  /// Or, most probably, Couenne?) that doesn't return good values
  /// from A (and triggers Valgrind errors and segfaults on ex3_1_1
  /// and others). While waiting for a fix, we'll use the row
  /// representation.

  /// Update: we should probably stick to #defining USE_ROW even when
  /// this is solved, as cuts from cs should also be added.

#define USE_ROW

#ifdef USE_ROW

  const CoinPackedMatrix *mat = si. getMatrixByRow ();

  int 
    m = mat -> getMajorDim (), // # rows
    n = mat -> getMinorDim (); // # cols

#else

  const CoinPackedMatrix *mat = si. getMatrixByCol ();

  int 
    n = mat -> getMajorDim (), // # cols
    m = mat -> getMinorDim (); // # rows

#endif

  const double
    *rlb = si.getRowLower (),
    *rub = si.getRowUpper ();

#ifdef USE_ROW

  /// These will be used

  int 
     nnz   = mat -> getNumElements (), // # nonzeros
     nnzC  = 0,
    *sta   = new int [n+1],
     nCuts = cs.sizeRowCuts ();

  // Count nonzeros in cs

  for (int i=0, ii = cs. sizeRowCuts (); ii--; i++) {

    const OsiRowCut        *cut     = cs. rowCutPtr (i);
    const CoinPackedVector &rowCoe  = cut -> row ();

    nnzC += rowCoe.getNumElements ();
  }

  int    *ind = new int    [nnz + nnzC];
  double *A   = new double [nnz + nnzC];

  /// these are the row-format originals
  {
    const double
      *rA   = mat -> getElements ();

    const int
      *rInd = mat -> getIndices      (),
      *rSta = mat -> getVectorStarts ();

    // copy rA, rInd, rSta into A, ind, sta

    CoinZeroN (sta, n+1);

    /////////////////////////////////////////////////////////

    // pre-fill starting positions with cardinalities of each column

    for (int i=nnz; i--;)
      ++ (sta [1 + *rInd++]);

    // fill sta with nonzeros from cs's OsiRowCuts

    for (int i=0, ii = cs. sizeRowCuts (); ii--; i++) {

      const OsiRowCut        *cut     = cs. rowCutPtr (i);
      const CoinPackedVector &rowCoe  = cut -> row ();
      const int              *indices = rowCoe.getIndices     ();
      int                     nnz     = rowCoe.getNumElements ();
      // Note: nnz redeclared here (no scoping problems)

      for (int i=nnz; i--;)
	++ (sta [1 + *indices++]);
    }

    rInd -= nnz;

    ////////////////////////////////////////////////////////

    // make sta cumulative

    for (int i=1; i<=n; i++)
      sta [i] += sta [i-1];

    // use space marked by sta to fill appropriate
    // indices/coefficients

    for (int i=0; i<m; i++) {

      // filling indices of row i

      int rowStart = rSta [i];

      for (int j = rowStart, jj = rSta [i+1] - rowStart; jj--; j++) {

	int &curSta = sta [rInd [j]];

	ind [curSta]   = i;
	A   [curSta++] = rA [j];
      }
      
      //printf ("\n");
    }

    // Add rowCuts from cs as well.

    for (int i=0, ii = cs. sizeRowCuts (); ii--; i++) {

      const OsiRowCut        *cut      = cs. rowCutPtr (i);
      //printf ("[cut] %4d [%g,%g]: ", m+i, cut -> lb (), cut -> ub ()); 

      const CoinPackedVector &rowCoe   = cut -> row ();
      const int              *indices  = rowCoe.getIndices  ();
      const double           *elements = rowCoe.getElements ();
      int                     nnz      = rowCoe.getNumElements ();
      // Note: nnz redeclared here (no scoping problems)

      for (int j=nnz; j--;) {

	//printf ("%+g x%d ", *elements, *indices); 

	int &curSta = sta [*indices++];

	ind [curSta]   = m+i;
	A   [curSta++] = *elements++;
      }
      //printf ("\n");
    }

    for (int i=n; --i;)
      sta [i] = sta [i-1];

    sta [0] = 0;

    //printf ("%d + %d = %d nonzeros\n", nnz, nnzC, nnz + nnzC);
  }

#else

  const double
    *A   = mat -> getElements ();

  const int
    *ind = mat -> getIndices      (),
    *sta = mat -> getVectorStarts ();

#endif

  /// Prepare vector for integrality test. Since many checks are done
  /// within combine(), it is worth to prepare one here

  bool *isInteger = new bool [n];
  for (int i=0, ii=n; ii--; i++)
    *isInteger++ = problem_ -> Var (i) -> isInteger ();
  isInteger -= n;

  // print out 

  // for (int i=0; i<n; i++) {

  //   printf ("x%04d - %5d -> %5d, %5d elements:", i, sta [i], sta [i+1], sta [i+1] - sta [i]);
  //   fflush (stdout);

  //   for (int j=0; j<sta [i+1] - sta [i]; j++) {
  //     printf ("(%d,%g) ", ind [sta [i] + j], A [sta [i] + j]);
  //     fflush (stdout);
  //   }
  //   printf ("\n");
  // }

  // For every column i, compare pairs of rows j and k with nonzero
  // coefficients.
  //
  // If the coefficients have the same sign and the inequalities are
  // both >= or both <=, skip this pair -- no tightening can come from
  // this pair (this doesn't mean that for another variable the
  // opposite may happen).

  double 
    *sa1 = new double [n], // contains dense representation of a1 i.e. lots of zeros
    *sa2 = new double [n]; //                                  a2

  CoinFillN (sa1, n, 0.);
  CoinFillN (sa2, n, 0.);

  for (int i=0; i<n; i++, sta++) {

    //printf ("x%d:\n", i);

    for   (int jj = *(sta+1) - *sta, j = *sta; jj--; j++)
      for (int kk = jj,              k = j+1;  kk--; k++) {

	if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	  break;

	register int 
	  indj = ind [j],
	  indk = ind [k];

	//printf ("  checking pair %d, %d\n", indj, indk);

	// should never happen, but if it does, just bail out
 	if ((indj >= m + nCuts) || (indj < 0) ||
 	    (indk >= m + nCuts) || (indk < 0)) {

	  if (nBadColMatWarnings++ < 1)
	    //	    jnlst_ -> Printf (J_STRONGWARNING, J_BOUNDTIGHTENING, " 
	    printf ("\
  Couenne: warning, matrix by row has nonsense indices.\n\
  This separator will now return without (column) cuts.\n\
  NOTE: further such inconsistencies won't be reported.\n");

	  delete [] sa1;
	  delete [] sa2;

	  delete [] sta;
	  delete [] ind;
	  delete [] A;

	  delete [] isInteger;

	  problem_ -> domain () -> pop ();

	  totalTime_     += CoinCpuTime () - now;
	  totalInitTime_ += CoinCpuTime () - now;

	  return; 
	}

	double 
	  prod = A [j] * A [k],
	  rlbj, rubj, rlbk, rubk;

	if (indj>=m) {
	  OsiRowCut *cut = cs.rowCutPtr (indj-m);
	  rlbj = cut -> lb ();
	  rubj = cut -> ub ();
	} else {
	  rlbj = rlb [indj];
	  rubj = rub [indj];
	}

	if (indk>=m) {
	  OsiRowCut *cut = cs.rowCutPtr (indk-m);
	  rlbk = cut -> lb ();
	  rubk = cut -> ub ();
	} else {
	  rlbk = rlb [indk];
	  rubk = rub [indk];
	}

	if (prod > 0.) { // same sign -- skip unless finite lb1/ub2 OR
			 // finite ub1/lb2. This is to avoid a situation
			 // in which all coefficients in this pair
			 // have the same sign

	  if (!(((rlbj > -COUENNE_INFINITY/10) && (rubk <  COUENNE_INFINITY/10)) || 
		((rubj <  COUENNE_INFINITY/10) && (rlbk > -COUENNE_INFINITY/10))))

	    continue;

	} else

	if ((prod < 0.) && // opposite sign -- multiply second
			   // inequality by -1 and repeat
	    !(((rlbj > -COUENNE_INFINITY/10) && (rlbk > -COUENNE_INFINITY/10)) || 
	      ((rubj <  COUENNE_INFINITY/10) && (rubk <  COUENNE_INFINITY/10))))

	  continue;

	pairs.insert (std::pair <int, int> (indj, indk));
      }
  }

  // pairs (h,k) are done. Now for each pair set new bounds, if possible ///////////////////////////

#ifdef USE_ROW
  const CoinPackedMatrix *rowA = mat;
#else
  const CoinPackedMatrix *rowA = si. getMatrixByRow ();
#endif

  const double
    *rA  = rowA -> getElements ();

  // TODO: no need for copy, though we need it to compare to old problem's bounds

  double
    *clb = CoinCopyOfArray (problem_ -> Lb (), n),
    *cub = CoinCopyOfArray (problem_ -> Ub (), n),
    *oldLB = CoinCopyOfArray (problem_ -> Lb (), n),
    *oldUB = CoinCopyOfArray (problem_ -> Ub (), n);

  const int
    *rInd = rowA -> getIndices      (),
    *rSta = rowA -> getVectorStarts ();

  int 
    ntightened = 0,
    ntrials    = 0,
    nCurTightened;

  // info about LP problem: upper bound, dual bound

  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (si.getAuxiliaryInfo ());

  int result = 0;

  // data structure for FBBT

  t_chg_bounds *chg_bds = new t_chg_bounds [n];

  // for (int i=0; i<n; i++) 
  //   if (problem_ -> Var (i) -> Multiplicity () <= 0) {
  //     chg_bds [i].setLower (t_chg_bounds::UNCHANGED);
  //     chg_bds [i].setUpper (t_chg_bounds::UNCHANGED);
  //   }

  // fills in chg_bds with changed bounds from previous BB node

  updateBranchInfo (si, problem_, chg_bds, info);

  // Comparing all pairs of inequalities is overkill if the comparison
  // has been done in a previous iteration: a pair of inequalities of
  // the old LP relaxation should be skipped unless some bounds have
  // been changed in at least one of the variables (of either
  // inequality). Moreover, all pairs of inequalities should be
  // considered where the first is from the LP relaxation and the
  // second is from the cuts added so far.

  // Repeat the scan as long as there is tightening and the bounding
  // box is nonempty.

  do {

    nCurTightened = 0;

    // scan all pairs. All are potential pairs of inequalities that
    // can give a better (combined) implied bound

    for (std::set <std::pair <int, int> >:: iterator p = pairs.begin (); p != pairs.end (); ++p) {

      if (CoinCpuTime () > problem_ -> getMaxCpuTime ())
	break;

      // indices of the two inequalities

      int 
	h = p -> first,
	k = p -> second,
	n1, n2;

      double l1, u1, l2, u2;

      const double *a1, *a2;

      const int	*ind1, *ind2;

      // set indices/counters depending on whether h or k are cuts or
      // inequalities

      if (h>=m) {

	const OsiRowCut *cut = cs.rowCutPtr (h-m);
	const CoinPackedVector &rowCoe  = cut -> row ();

	l1 = cut -> lb ();
	u1 = cut -> ub ();
	n1   = rowCoe. getNumElements ();
	ind1 = rowCoe. getIndices     ();
	a1   = rowCoe. getElements    ();

      } else {

	l1 = rlb [h];
	u1 = rub [h];
	n1 = rSta [h+1] - rSta [h];
	ind1 = rInd + rSta [h];
	a1 = rA + rSta [h];
      }

      ////////////////////////////////////////////////////////////////

      if (k>=m) {

	OsiRowCut *cut = cs.rowCutPtr (k-m);
	const CoinPackedVector &rowCoe  = cut -> row ();

	l2 = cut -> lb ();
	u2 = cut -> ub ();
	n2   = rowCoe. getNumElements ();
	ind2 = rowCoe. getIndices     ();
	a2   = rowCoe. getElements    ();

      } else {

	l2 = rlb [k];
	u2 = rub [k];
	n2 = rSta [k+1] - rSta [k];
	ind2 = rInd + rSta [k];
	a2 = rA + rSta [k];
      }

      // If both indices are from the LP relaxation, only check if
      // there is at least one changed bound in the variables of one
      // (or both) of them

      if ((h < m) && (k < m) && !firstCall_) {

	bool new_bound = false;

	for (int i=n1; i--;) {

	  t_chg_bounds &chg = chg_bds [ind1 [i]];

	  if ((chg. lower () != t_chg_bounds::UNCHANGED) ||
	      (chg. upper () != t_chg_bounds::UNCHANGED))

	    new_bound = true;
	}

	if (!new_bound) 
	  for (int i=n2; i--;) {

	    t_chg_bounds &chg = chg_bds [ind2 [i]];

	    if ((chg. lower () != t_chg_bounds::UNCHANGED) ||
		(chg. upper () != t_chg_bounds::UNCHANGED))

	      new_bound = true;
	  }

	if (!new_bound) continue;
      }

      // fill in sa1 and sa2

      for (int i=n1; i--;) sa1 [ind1 [i]] = a1 [i];
      for (int i=n2; i--;) sa2 [ind2 [i]] = a2 [i];

      // A few cases here on the inequalities' bounds, which may
      // determine a change of sign in the second when combining the
      // two: for example, if we have -inf < ax < b and c < dx < +inf,
      // it makes sense to combine -inf < ax < b and -inf < - dx < -c
      // as the upper bound will create a finite convex combination.

      if ((u1 <   COUENNE_INFINITY && u2 <   COUENNE_INFINITY) ||
	  (l1 > - COUENNE_INFINITY && l2 > - COUENNE_INFINITY))

	result = combine (problem_, n1, n2, ind1, ind2, sa1, sa2, a1, a2, clb, cub, l1, l2, u1, u2, isInteger, 1);

      if (result < 0)
	break;

      nCurTightened += result;
      result = 0;

      // fill in sa2 with opposite values
      for (int i=n2; i--;) 
	sa2 [ind2 [i]] = - a2 [i];

      if ((u1 <   COUENNE_INFINITY && l2 > - COUENNE_INFINITY) ||
	  (l1 > - COUENNE_INFINITY && u2 <   COUENNE_INFINITY))

	// do NOT invert l2 and u2, this is done in combine
	result = combine (problem_, n1, n2, ind1, ind2, sa1, sa2, a1, a2, clb, cub, l1, l2, u1, u2, isInteger, -1);

      if (result < 0)
	break;

      nCurTightened += result;

      // clean sa1 and sa2

      for (int i=n1; i--;) sa1 [ind1 [i]] = 0.;
      for (int i=n2; i--;) sa2 [ind2 [i]] = 0.;
    }

    if (result < 0) 
      break;

    int objInd = problem_ -> Obj (0) -> Body () -> Index ();

    if (nCurTightened &&
	(objInd >= 0) && 
	babInfo && 
	babInfo -> babPtr ()) {

#ifdef DEBUG
      printf ("FBBT\n");
#endif

      CouNumber
	UB      = babInfo -> babPtr () -> model (). getObjValue(),
	LB      = babInfo -> babPtr () -> model (). getBestPossibleObjValue (),
	primal0 = problem_ -> Ub (objInd), 
	dual0   = problem_ -> Lb (objInd);

      // Do one round of BT

      if ((UB < COUENNE_INFINITY) && 
	  (UB < primal0 - COUENNE_EPS)) { // update primal bound (MIP)

	problem_ -> Ub (objInd) = UB;
	chg_bds [objInd].setUpper (t_chg_bounds::CHANGED);
      }

      if ((LB > - COUENNE_INFINITY) && 
	  (LB > dual0 + COUENNE_EPS)) { // update dual bound
	problem_ -> Lb (objInd) = LB;
	chg_bds [objInd].setLower (t_chg_bounds::CHANGED);
      }

      // change chg_bds for all new bounds stored in cs

      /*for (int i = cs. sizeColCuts (); i--;) {

	const CoinPackedVector
	  &lbs = cs. colCutPtr (i) -> lbs (),
	  &ubs = cs. colCutPtr (i) -> ubs ();

	// copy lbs

	const int *indices = lbs. getIndices ();

	for (int j = lbs. getNumElements (); j--; indices++)
	  chg_bds [*indices].setLower (t_chg_bounds::CHANGED);

	// copy ubs

	indices = ubs. getIndices ();

	for (int j = ubs. getNumElements (); j--; indices++)
	  chg_bds [*indices].setUpper (t_chg_bounds::CHANGED);
      }
      */

      if (!(problem_ -> btCore (chg_bds))) {

	//problem infeasible, add IIS of size 2

	result = -1;
	break;

      } else {

	// update tightened bounds from problem to clb

	for (int i=0; i<n; i++) {

	  if (problem_ -> Lb (i) > clb [i] + COUENNE_EPS) clb [i] = problem_ -> Lb (i);
	  if (problem_ -> Ub (i) < cub [i] - COUENNE_EPS) cub [i] = problem_ -> Ub (i);
	}
      }
    }

    ntightened += nCurTightened;

  } while (++ntrials < nMaxTrials_ && nCurTightened);

  totalInitTime_ += CoinCpuTime () - now;

  // check if bounds improved, in that case create OsiColCuts

  if (result >= 0 && ntightened) {

    // check old and new bounds

    int 
      *indLB = new int [n],
      *indUB = new int [n],
      ntightenedL = 0,
      ntightenedU = 0;

    double 
      *valLB = new double [n],
      *valUB = new double [n];

    for (int i=0; i<n; i++) {

      if (clb [i] > oldLB [i]) {

	if (problem_ -> bestSol () &&
	    problem_ -> bestSol () [i] > oldLB [i] &&
	    problem_ -> bestSol () [i] < clb   [i] - 1e-12) {

	  jnlst_ -> Printf (J_ERROR, J_COUENNE, 
			    "Warning, twoImplBounds new LB cuts optimal solution: LB x_%d = %g --> %g, opt %g (diff: %g)\n",
			    i, oldLB [i], clb [i], problem_ -> bestSol () [i], clb [i] - problem_ -> bestSol () [i]);

	}

	//printf ("L %4d: %g -> %g\n", i, oldLB [i], clb [i]);
	indLB [ntightenedL]   = i;
	valLB [ntightenedL++] = clb [i];
      }

      if (cub [i] < oldUB [i]) {

	if (problem_ -> bestSol () &&
	    problem_ -> bestSol () [i] < oldUB [i] &&
	    problem_ -> bestSol () [i] > cub   [i] + 1e-12) {

	  jnlst_ -> Printf (J_ERROR, J_COUENNE, 
			    "Warning, twoImplBounds new UB cuts optimal solution: UB x_%d = %g --> %g, opt %g (diff: %g)\n",
			    i, oldUB [i], cub [i], problem_ -> bestSol () [i], problem_ -> bestSol () [i] - cub [i]);

	}

	//printf ("U %4d: %g -> %g\n", i, oldUB [i], cub [i]);
	indUB [ntightenedU]   = i;
	valUB [ntightenedU++] = cub [i];
      }
    }

    if (ntightenedL || ntightenedU) {

      OsiColCut newBound;

      newBound.setLbs (ntightenedL, indLB, valLB);
      newBound.setUbs (ntightenedU, indUB, valUB);

      //newBound.print ();

      cs.insert (newBound);
    }

    delete [] indLB;
    delete [] indUB;
    delete [] valLB;
    delete [] valUB;

    ntightened = ntightenedL + ntightenedU;
  }

  else 

    if (result < 0)
      WipeMakeInfeas (cs);

  delete [] clb;
  delete [] cub;
  delete [] oldLB;
  delete [] oldUB;
  delete [] sa1;
  delete [] sa2; 
  delete [] chg_bds;

  delete [] isInteger;

  problem_ -> domain () -> pop ();

#ifdef USE_ROW
  delete [] A;
  delete [] ind;
  delete [] (sta-n);
#endif

  if (firstCall_)
    firstCall_ = false;

  totalTime_ += CoinCpuTime () - now;

  if (info.level <= 0)
    jnlst_ -> Printf (J_ERROR, J_COUENNE, "%d improved bounds\n", ntightened);
}
}
