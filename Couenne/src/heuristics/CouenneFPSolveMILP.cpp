/* $Id$
 *
 * Name:    CouenneFPSolveMILP.cpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Solve the MILP within the Feasibility Pump 
 * 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CbcModel.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneFeasPump.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"

using namespace Couenne;

/// find integer (possibly NLP-infeasible) point isol closest
/// (according to the l-1 norm of the hessian) to the current
/// NLP-feasible (but fractional) solution nsol
CouNumber CouenneFeasPump::solveMILP (CouNumber *nSol0, CouNumber *&iSol) {

  printf ("FP:solveMILP\n");

  // The problem is of the form
  //
  // min  f(x)
  // s.t. Ax >= b
  //      x_i in Z, i in N
  //
  // where N is the index set of integer variables and f(x) is one of
  // the following:
  //
  // 1) sum {i in Vars} |x_i - x_i^0|
  // 2) sum {i in N}    |x_i - x_i^0|
  // 3) sum {i in Vars} |P^i (x - x^0)| 
  // 4) sum {i in N}    |P^i (x - x^0)|
  //
  // where is x^0 is the optimal solution of a NLP problem. In the
  // last two, the l-1 norm is multiplied by the i-th row of a matrix
  // P obtained from the Hessian H of the Lagrangian of the problem. H
  // is a positive semidefinite matrix only on the null space of the
  // gradient g, i.e.
  //
  // x' H x >= 0 for all x: g'x = 0
  //
  // (where x' is x transposed), therefore we add a quadratic term to
  // make it positive semidefinite everywhere else. In order to do so,
  // the quadratic term (alpha is a scalar)
  //
  // alpha (g'x)^2
  // 
  // is zero in the null space of g, and it is strictly positive
  // everywhere else. Hence, we need an alpha so that 
  //
  // P = (H + alpha g g') 
  //
  // is PSD. In general, we might have a parameter beta in [0,1] such
  // that
  // 
  // P = beta I + (1-beta) (H + alpha g g')
  //
  // so that we can balance the hessian and the distance.

  CoinPackedMatrix P;

  if (!milp_) {

    milp_ = model_ -> solver () -> clone ();

    // copy bounding box and current solution to the problem (better
    // linearization cuts)
    problem_ -> domain () -> push (milp_ -> getNumCols (),
				   milp_ -> getColSolution (),
				   milp_ -> getColLower (),
				   milp_ -> getColUpper ());

    // no data is available so far, retrieve it from the MILP solver
    // used as the linearization

    // Construct an (empty) Hessian. It will be modified later, but
    // the changes should be relatively easy for the case when
    // betaMILP > 0 and there are no changes if betaMILP_ == 0

    P. setDimensions (problem_ -> nVars (), 
		      problem_ -> nVars ());

    // The MILP has to be changed the first time it is used.
    //
    // Suppose Ax >= b has m inequalities. In order to solve the
    // problem above, we need q new variables z_i and 2q inequalities
    //
    //   z_i >=   P^i (x - x^0)  or  P^i x - z_i <= P^i x^0 (*)
    //   z_i >= - P^i (x - x^0)                             (**)
    // 
    // (the latter being equivalent to
    //
    // - z_i <=   P^i (x - x^0)  or  P^i x + z_i >= P^i x^0 (***)
    //
    // so we'll use this instead as most coefficients don't change)
    // for each i, where q is the number of variables involved (either
    // q=|N|, the number of integer variables, or q=n, the number of
    // variables).
    //
    // We need to memorize the number of initial inequalities and of
    // variables, so that we know what (coefficients and rhs) to
    // change at every iteration.

    // Add q variables, each with coefficient 1 in the objective

    CoinPackedVector vec;

    for (int i=0; i<problem_ -> nVars (); i++)
      if (!compDistInt_ || milp_ -> isInteger (i))
	// (empty) coeff col vector, lb, ub, obj coeff
	milp_ -> addCol (vec, 0, COIN_DBL_MAX, 1.); 

    // Set to zero all other variables' obj coefficient. This means we
    // just do it for the single variable in the reformulated
    // problem's linear relaxation (all other variables do not appear
    // in the objective)

    milp_ -> setObjCoeff (problem_ -> Obj (0) -> Body () -> Index (), 0.);

  } else {

    // delete last rows and add them from scratch (common block below)

    int 
       nDeleted = compDistInt_ ? milp_ -> getNumIntegers () : 
                                 problem_ -> nVars (),
      *deleted  = new int [2*nDeleted],
       nCurRow  = milp_ -> getNumRows ();

    for (int i=2*nDeleted; i--;)
      deleted [i] = --nCurRow;

    milp_ -> deleteRows (nDeleted, deleted);

    printf ("FP:deleting last lines\n");
    milp_ -> writeLp ("fp-milp-del"); // !!

    delete [] deleted;
  }

  printf ("FP:add MILP ineq\n");
  milp_ -> writeLp ("fp-milp0");

  // Add 2q inequalities

  for (int i=0, j=problem_ -> nVars (); i<problem_ -> nVars (); i++)

    if (!compDistInt_ || milp_ -> isInteger (i)) {

      double val = 1.;
      CoinPackedVector vec (1, &i, val);

      if (betaMILP_ > 0.) {

	// reserved for non-UnitMatrix hessian (i.e. betaMILP_ > 0)

      }

      CoinPackedVector x0 (problem_ -> nVars (), milp_ -> getColSolution ());

      // right-hand side equals <P^i,x^0>
      double PiX0 = sparseDotProduct (vec, x0); 

      vec.insert     (j, -1.); milp_ -> addRow (vec, -COIN_DBL_MAX,         PiX0); // (*)
      vec.setElement (1, +1.); milp_ -> addRow (vec,          PiX0, COIN_DBL_MAX); // (***)

      ++j; // index of variable within problem (plus nVars_)
    }

  milp_ -> writeLp ("fp-milp");

  milp_ -> branchAndBound ();

  iSol = CoinCopyOfArray (milp_ -> getColSolution (), 
			  problem_ -> nVars ());

  return milp_ -> getObjValue ();
}
