/* $Id$
 *
 * Name:    CouenneFPFindSolution.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Find solution by looping through MILP solvers/heuristics
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinTime.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneFPpool.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

#include "cons_rowcuts.h"

#ifdef COIN_HAS_SCIP
#include "scip/scip.h"
#endif

using namespace Couenne;

/// find a feasible or optimal solution of MILP
double CouenneFeasPump::findSolution (double* &sol, int niter, int* nsuciter) {

  /// As found on the notes, these methods can be used, from the most
  /// expensive and accurate (exact) method to a cheap, inexact one:
  ///
  /// 1. Solve a MILP relaxation with Manhattan distance as objective
  /// 2. Partially solve the MILP with emphasis on good solutions
  /// 3. Apply RENS to 1
  /// 4. Use Objective FP 2.0 for MILPs
  /// 5. round-and-propagate
  /// 6. choose from pool, see 4
  /// 7. random perturbation

  // What order should we use? I suggest we use priorities, assigned
  // at the beginning but changeable in the event of multiple failures
  // (or successes) of a given method.
  //
  // Rule of thumb: 
  //
  // 1) Assign all methods i a number p[i] (for instance those in the
  //    list above)
  //
  // 2) Call each in the order defined by p[i], return a solution if
  //    found, otherwise proceed to next method
  //
  // 3) If K consecutive successes at finding new solution (not
  //    necessarily new best feasible), --p[i]
  //
  // 4) if H consecutive failutes, ++p[i]

  double obj;

  /// solve MILP 

#ifdef COIN_HAS_SCIP

  if (useSCIP_ && problem_ -> nIntVars () > 0) { // if LP, use Clp below

    SCIP_RETCODE retcode = ScipSolve (sol, niter, nsuciter, obj);

    if (retcode != SCIP_OKAY) {

      printf ("SCIP did not return successfully\n");
      return COIN_DBL_MAX;
    }
  } else

#endif      
  {

     if (problem_ -> nIntVars () > 0) milp_ -> branchAndBound ();
     else                             milp_ -> initialSolve ();

     if (!sol)
       sol = new CouNumber [problem_ -> nVars ()];

     if (milp_ -> getColSolution ())
       CoinCopyN (milp_ -> getColSolution (), problem_ -> nVars (), sol);
     else {

       if (sol)
	 delete [] sol;
       sol = NULL;
     }

     obj = milp_ -> getObjValue ();
  }

  return obj;
}

/// initialize MILP solvers if needed
void CouenneFeasPump::init_MILP () {}
