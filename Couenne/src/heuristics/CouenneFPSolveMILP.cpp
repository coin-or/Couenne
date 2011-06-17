/* $Id$
 *
 * Name:    CouenneFPSolveMILP.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Solve the MILP within the Feasibility Pump 
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"

#include "CouenneConfig.h"
#include "CouenneFeasPump.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"

#include "CouenneFPpool.hpp"

#ifdef COIN_HAS_SCIP
/* general SCIP includes */
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#endif

using namespace Couenne;

#define NUMERICS_THRES 1e19

#ifdef COIN_HAS_SCIP
void CouenneFeasPump::checkInfinity(SCIP *scip, SCIP_Real val, double infinity){
   if( SCIPisInfinity(scip, val) && val < infinity)
      printf("Warning: %g will be considered to be Infinity by SCIP\n", val);
}
#endif


/// create clone of MILP and add variables for special objective
OsiSolverInterface *createCloneMILP (const CouenneFeasPump *fp, CbcModel *model, bool isMILP);


/// modify MILP or LP to implement distance by adding extra rows (extra cols were already added by createCloneMILP)
void addDistanceConstraints (const CouenneFeasPump *fp, OsiSolverInterface *lp, double *sol, bool isMILP);


/// find integer (possibly NLP-infeasible) point isol closest
/// (according to the l-1 norm of the Hessian) to the current
/// NLP-feasible (but fractional) solution nsol
CouNumber CouenneFeasPump::solveMILP (CouNumber *nSol0, CouNumber *&iSol, int niter, int* nsuciter, int depth) {

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
  // so that we can balance the Hessian and the distance.

  bool firstCall = (milp_ != NULL); // true if this is the first call to
  	                            // solveMILP; initialization will be
  	                            // necessary
  
  if (!milp_) {

    // create MILP

    milp_ = createCloneMILP (this, model_, true);

    // Post-processing LP: persistent if FP_DIST_POST, created on the
    // fly if FP_DIST_INT and numerics, not created if FP_DIST_ALL
    //
    // set up an LP as a copy of the original MILP. Don't do the
    // same for FP_DIST_INT as it is only necessary upon numerical
    // problems, which might not happen

    if (compDistInt_ == FP_DIST_POST)
      postlp_ = createCloneMILP (this, model_, false);
  }

  // Construct an (empty) Hessian. It will be modified later, but
  // the changes should be relatively easy for the case when
  // betaMILP > 0 and there are no changes if betaMILP_ == 0

  CoinPackedMatrix P;

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

  // ADD CODE HERE...

  // Add 2q inequalities

  int nInitRows = milp_ -> getNumRows ();

  CouNumber * nlpSolExp;

  if (nSol0) {

    nlpSolExp = new CouNumber [problem_ -> nVars ()];

    CoinCopyN (nSol0, problem_ -> nOrigVars (), nlpSolExp);
    problem_ -> getAuxs (nlpSolExp);

  } else 
    nlpSolExp = CoinCopyOfArray (milp_ -> getColSolution (), 
				 problem_ -> nVars ());

  // create constraints to define l_1 distance objective function
  addDistanceConstraints (this, milp_, nlpSolExp, true);

  delete [] nlpSolExp;

  int nFinalRows = milp_ -> getNumRows ();

  // The MILP is complete. We have several ways of solving it, or
  // better, to find feasible solutions to it. We have to interface
  // with each of them once at the very beginning, and later we loop
  // through them in order to find a feasible solution.

  if (firstCall)
    init_MILP ();

  if (false) {
    static int cntr = 0;
    char filename [30];
    sprintf (filename, "fp-milp%04d", cntr++);
    milp_ -> writeLp (filename);
  }

  double obj = findSolution (iSol, niter, nsuciter, depth);

  if ((nSol0 && iSol) &&
      (problem_ -> Jnlst () -> ProduceOutput (Ipopt::J_STRONGWARNING, J_NLPHEURISTIC))) {

    double dist = 0.;
    int nNonint = 0;

    for (int i = 0; i < problem_ -> nVars (); ++i) {

      if (problem_ -> Var (i) -> isInteger () &&
	  (fabs (iSol [i] - ceil (iSol [i] - .5)) > 1e-4))
	++nNonint;

      dist += 
	(iSol [i] - nSol0 [i]) * 
	(iSol [i] - nSol0 [i]);
    }

    printf ("FP: solution: distance %g, %d nonintegers\n", dist, nNonint);
  }

  //
  // POST PROCESSING 
  //
  // (if we got a solution from MILP, otherwise bail out)
  //

  if (iSol &&
      (compDistInt_ != FP_DIST_ALL)) {

    // check iSol for numerics (i.e. components whose fabs () is large,
    // or >= 1e20) or post-process to obtain a second solution by fixing
    // the integer coordinates and solving the resulting LP

      bool numerics = false;

      if (compDistInt_ == FP_DIST_INT) {

	for (std::vector <exprVar *>::iterator i = problem_ -> Variables (). begin (); 
	     i != problem_ -> Variables (). end (); ++i)

	  if ((  (*i) -> Multiplicity () > 0) &&
	      ! ((*i) -> isInteger    ())     &&
	      (fabs (iSol [(*i) -> Index ()]) > NUMERICS_THRES)) {

	    numerics = true;
	    break;
	  }
      }

      if (numerics || (compDistInt_ == FP_DIST_POST)) {

	// solve LP where integer variables have been fixed:
	//
	// a) check if postlp_ exists yet
	// 0) save integer bounds
	// 1) fix integer variables 
	// 2) add variables and inequalities
	// 3) solve LP
	// 4) if optimal, save solution
	// 5) restore IP bounds
	// 6) delete variables

	if (!postlp_)
	  postlp_ = createCloneMILP (this, model_, false);

	int nvars = postlp_ -> getNumCols ();

	// save integer bounds to restore them later
	double
	  *saveLB = CoinCopyOfArray (postlp_ -> getColLower (), nvars),
	  *saveUB = CoinCopyOfArray (postlp_ -> getColUpper (), nvars),
	  *newLB  = CoinCopyOfArray (postlp_ -> getColLower (), nvars),
	  *newUB  = CoinCopyOfArray (postlp_ -> getColUpper (), nvars);

	// fix integer variables

	for (int i = problem_ -> nVars (); i--;)
	  if (milp_ -> isInteger (i))
	    newLB [i] = newUB [i] = milp_ -> getColSolution () [i];

	postlp_ -> setColLower (newLB);
	postlp_ -> setColUpper (newUB);

	// add inequalities

	int nInitRowsLP  = postlp_ -> getNumRows ();
	addDistanceConstraints (this, postlp_, iSol, false);
	int nFinalRowsLP = postlp_ -> getNumRows ();

	// Solve the LP, obtain closest point with integer variables fixed

	postlp_ -> initialSolve ();

	// save as new solution

	if (postlp_ -> isProvenOptimal ()) 
	  CoinCopyN (postlp_ -> getColSolution (), problem_ -> nVars (), iSol);

	postlp_ -> setColLower (saveLB);
	postlp_ -> setColUpper (saveUB);

	// delete temp data

	delete [] saveLB;
	delete [] saveUB;
	delete [] newLB;
	delete [] newUB;

	// delete added rows

	int 
	  nDeleted = nFinalRowsLP - nInitRowsLP,
	  *deleted  = new int [nDeleted],
	  nCurRow  = nInitRowsLP;

	for (int i = nDeleted; i--;)
	  deleted [i] = nCurRow++;

	postlp_ -> deleteRows (nDeleted, deleted);

	delete [] deleted;
      }
  }

  // delete last rows and add them from scratch (common block below)

  int 
    nDeleted = nFinalRows - nInitRows,
   *deleted  = new int [nDeleted],
    nCurRow  = nInitRows;

  for (int i = nDeleted; i--;)
    deleted [i] = nCurRow++;

  milp_ -> deleteRows (nDeleted, deleted);

  delete [] deleted;

  return obj;// milp_ -> getObjValue ();
}


/// create clone of MILP and add variables for special objective
OsiSolverInterface *createCloneMILP (const CouenneFeasPump *fp, CbcModel *model, bool isMILP) {

  OsiSolverInterface *lp = model -> solver () -> clone ();

  // no data is available so far, retrieve it from the MILP solver
  // used as the linearization

  // Add q variables, each with coefficient 1 in the objective

  CoinPackedVector vec;

  for (int i = fp -> Problem () -> nVars (), j = 0; i--; ++j) {

    // column has to be added if:
    //
    // creating MIP AND (integer variable OR FP_DIST_ALL)
    // creating LP  AND fractional variable

    bool intVar = lp -> isInteger (j);

    if ((isMILP && (intVar || (fp -> compDistInt () == CouenneFeasPump::FP_DIST_ALL)))
	||
	(!isMILP && !intVar))
      // (empty) coeff col vector, lb, ub, obj coeff
      lp -> addCol (vec, 0., COIN_DBL_MAX, 1.); 
  }

  // Set to zero all other variables' obj coefficient. This means we
  // just do it for the single variable in the reformulated
  // problem's linear relaxation (all other variables do not appear
  // in the objective)

  lp -> setObjCoeff (fp -> Problem () -> Obj (0) -> Body () -> Index (), 0.);

  return lp;
}


/// modify MILP or LP to implement distance by adding extra rows (extra cols were already added by createCloneMILP)
void addDistanceConstraints (const CouenneFeasPump *fp, OsiSolverInterface *lp, double *sol, bool isMILP) {

  CoinPackedVector x0 (fp -> Problem () -> nVars (), sol);

  for (int i = 0, j = fp -> Problem () -> nVars (), k = j; k--; ++i) {

    // two rows have to be added if:
    //
    // amending MIP AND (integer variable OR FP_DIST_ALL)
    // amending ssLP  AND fractional variable

    bool intVar = lp -> isInteger (i);

    if ((isMILP && (intVar || (fp -> compDistInt () == CouenneFeasPump::FP_DIST_ALL)))
	||
	(!isMILP && !intVar)) {

      // create vector with single entry of 1 at i-th position 
      double val = 1.;
      CoinPackedVector vec (1, &i, val);

      if (fp -> betaMILP () > 0.) {

	// reserved for non-UnitMatrix hessian (i.e. betaMILP_ > 0)
      }

      // right-hand side equals <P^i,x^0>
      double PiX0 = sparseDotProduct (vec, x0); 

      vec.insert     (j, -1.); lp -> addRow (vec, -COIN_DBL_MAX,         PiX0); // (*)
      vec.setElement (1, +1.); lp -> addRow (vec,          PiX0, COIN_DBL_MAX); // (***)

      ++j; // index of variable within problem (plus nVars_)

    } else if (intVar) {

      // fix integer variable to its value in iSol      
    }
  }
}
