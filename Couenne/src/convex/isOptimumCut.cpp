/* $Id$ */
/*
 * Name:    isOptimumCut.cpp
 * Author:  Pietro Belotti
 * Purpose: check if known optimal solution (read from .txt) is
 *          erroneously cut by the row/col cuts we've added
 *
 * (C) Carnegie-Mellon University, 2009.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CglCutGenerator.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"


bool isOptimumCut (const CouNumber *opt, OsiCuts &cs, CouenneProblem *p) {

  bool retval = false;

  // Column cuts ////////////////////////////////////////////////////////////////////////

  if (cs.sizeColCuts ()) {

    //printf ("Checking col cuts:\n");

    for (int i = cs.sizeColCuts (); i--;) {

      // lower bounds

      const CoinPackedVector &lbs = cs.colCutPtr (i) -> lbs ();
      const int    *lindices = lbs.getIndices ();
      const double *lvalues  = lbs.getElements ();

      for (int j = lbs.getNumElements (); j--;) {
	register double lb  = *lvalues++;
	register int    ind = *lindices++;

	if (lb > opt [ind] + COUENNE_EPS) {
	  printf ("################################## new lb [%d] = %g cuts opt %g by %g\n",
		  ind, lb, opt [ind], lb - opt [ind]);
	  retval = true;
	}
      }

      // upper bounds

      const CoinPackedVector &ubs = cs.colCutPtr (i) -> ubs ();
      const int    *uindices = ubs.getIndices ();
      const double *uvalues  = ubs.getElements ();

      for (int j = ubs.getNumElements (); j--;) {
	register double ub  = *uvalues++;
	register int    ind = *uindices++;

	if (ub < opt [ind] - COUENNE_EPS) {
	  printf ("################################## new ub [%d] = %g cuts opt %g by %g\n",
		  ind, ub, opt [ind], opt [ind] - ub);
	  retval = true;
	}
      }
    }
  }

  // Row cuts ///////////////////////////////////////////////////////////////////////////

  if (cs.sizeRowCuts ()) {

    //printf ("Checking row cuts:\n");

    for (int jj=0; jj < cs.sizeRowCuts (); jj++) {

      OsiRowCut        *cut = cs.rowCutPtr (jj);
      CoinPackedVector  row = cut -> row ();

      int           n   = cut -> row (). getNumElements();
      const double *el  = row. getElements ();
      const int    *ind = row. getIndices ();

      double        lb  = cut -> lb ();
      double        ub  = cut -> ub ();

      double lhs = 0;

      while (n--) 
	lhs += el [n] * opt [ind [n]];


      if ((lhs < lb - COUENNE_EPS) || 
	  (lhs > ub + COUENNE_EPS)) {

	printf ("################################## new cut [%d] [%g,%g] cuts opt %g by %g:",
		jj, lb, ub, lhs, CoinMax (lb - lhs, lhs - ub));

	cut -> print ();
	retval = true;
      }
    }
  }

  if (retval) {

    printf ("== genrowcuts on LP =============");

    for (int i = 0; i < p -> nVars (); i++) {
      if (!(i % 3)) 
	printf ("\n");
      if (p -> Var (i) -> Multiplicity () > 0)
	printf ("%3d %+10.3g [%+10.3g,%+10.3g] ", i,
		p -> X  (i),
		p -> Lb (i),
	      p -> Ub (i));
    }

    printf ("\n=============================\n");
  }

  return retval;
}
