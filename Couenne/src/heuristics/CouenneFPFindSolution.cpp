/* $Id$
 *
 * Name:    CouenneFPFindSolution.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Find solution by looping through MILP solvers/heuristics
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneFeasPump.hpp"

using namespace Couenne;

/// find a feasible or optimal solution of MILP
void CouenneFeasPump::findSolution () {

  /// as found on the notes, these methods can be used, from the most
  /// expensive and accurate (exact) method to a cheap, inexact one:
  ///
  /// 1. Solve a MILP relaxation with Manhattan distance as objective
  /// 2. Apply RENS to 1
  /// 3. Use Objective FP 2.0 for MILPs
  /// 4. round-and-propagate
  /// 5. choose from pool, see 4
  /// 6. random perturbation

  // What order should we use? I suggest we use priorities, assigned
  // at the beginning but changeable in the event of multiple failures
  // (or successes) of a given method.
  //
  // Rule of thumb: 
  //
  // 1) Assign all methods i a number p[i] (for instance those in the
  //    list above)
  //
  // 2) Call each in the order define by p[i], return a solution if
  //    found, otherwise proceed to next method
  //
  // 3) If K consecutive successes at finding new solution (not
  //    necessarily new best feasible), --p[i]
  //
  // 4) if H consecutive failutes, ++p[i]

  /// solve MILP 

  milp_ -> branchAndBound ();
}

/// initialize MILP solvers if needed
void CouenneFeasPump::init_MILP () {

}
