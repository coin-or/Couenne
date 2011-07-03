/* $Id$
 *
 * Name:    CouenneFPcreateMILP.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: create the MILP within the Feasibility Pump 
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CbcModel.hpp"

#include "IpLapack.hpp"

#include "CouenneSparseMatrix.hpp"
#include "CouenneTNLP.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"

#define COUENNE_EIG_RATIO .1 // how much smaller than the largest eigenvalue should the minimum be set at?

using namespace Couenne;

/// computes square root of a CouenneSparseMatrix
void ComputeSquareRoot (const CouenneFeasPump *fp, CouenneSparseMatrix *hessian, CoinPackedVector *P);

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
      // (empty) coeff col vector, lb = 0, ub = inf, obj coeff
      lp -> addCol (vec, 0., COIN_DBL_MAX, 1.); 
  }

  // Set to zero all other variables' obj coefficient. This means we
  // just do it for the single variable in the reformulated
  // problem's linear relaxation (all other variables do not appear
  // in the objective)

  int objInd = fp -> Problem () -> Obj (0) -> Body () -> Index ();

  if (objInd >= 0)
    lp -> setObjCoeff (objInd, 0.);

  return lp;
}


/// modify MILP or LP to implement distance by adding extra rows (extra cols were already added by createCloneMILP)
void addDistanceConstraints (const CouenneFeasPump *fp, OsiSolverInterface *lp, double *sol, bool isMILP) {

  // Construct an (empty) Hessian. It will be modified later, but
  // the changes should be relatively easy for the case when
  // multHessMILP_ > 0 and there are no changes if multHessMILP_ == 0

  int n = fp -> Problem () -> nVars ();

  CoinPackedVector *P = new CoinPackedVector [n];

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

  CoinPackedVector x0 (n, sol);

  // set objective coefficient if we are using a little Objective FP

  if (isMILP && (fp -> multObjFMILP () > 0.)) {

    int objInd = fp -> Problem () -> Obj (0) -> Body () -> Index ();

    if (objInd >= 0)
      lp -> setObjCoeff (objInd, fp -> multObjFMILP ());
  }

  if (isMILP && 
      (fp -> multHessMILP () > 0.) &&
      (fp -> nlp () -> optHessian ())) {

    // P is a convex combination, with weights multDistMILP_ and
    // multHessMILP_, of the distance and the Hessian respectively

    // obtain optHessian and compute its square root

    CouenneSparseMatrix *hessian = fp -> nlp () -> optHessian ();

    ComputeSquareRoot (fp, hessian, P);

  } else {

    // simply set P = I

    for (int i=0; i<n; i++)
      P[i].insert (i, 1.); 
  }

  // Add 2q inequalities

  for (int i = 0, j = n, k = j; k--; ++i) {

    // two rows have to be added if:
    //
    // amending MIP AND (integer variable OR FP_DIST_ALL)
    // amending ssLP  AND fractional variable

    bool intVar = lp -> isInteger (i);

    if ((isMILP && (intVar || (fp -> compDistInt () == CouenneFeasPump::FP_DIST_ALL)))
	||
	(!isMILP && !intVar)) {

      // create vector with single entry of 1 at i-th position 
      CoinPackedVector &vec = P [i];

      if (vec.getNumElements () == 0)
	continue;

      // right-hand side equals <P^i,x^0>
      double PiX0 = sparseDotProduct (vec, x0); 

      // j is the index of the j-th extra variable z_j, used for z_j >=  P (x - x0)  ===> z_j - Px >= - Px_0 ==> -z_j + Px <= Px_0
      // Second inequality is                                    z_j >= -P (x - x0)                          ==>  z_j + Px >= Px_0
      vec.insert     (j,                         -1.); lp -> addRow (vec, -COIN_DBL_MAX,         PiX0); // (*)
      vec.setElement (vec.getNumElements () - 1, +1.); lp -> addRow (vec,          PiX0, COIN_DBL_MAX); // (***)

      ++j; // index of variable within problem (plus nVars_)

    } else if (intVar) { // implies (!isMILP)

      // fix integer variable to its value in iSol      

#define INT_LP_BRACKET 0

      lp -> setColLower (i, sol [i] - INT_LP_BRACKET);
      lp -> setColUpper (i, sol [i] + INT_LP_BRACKET);
    }
  }

  delete [] P;
}


#define GRADIENT_WEIGHT 1

void ComputeSquareRoot (const CouenneFeasPump *fp, 
			CouenneSparseMatrix *hessian, 
			CoinPackedVector *P) {
  int 
    objInd = fp -> Problem () -> Obj (0) -> Body () -> Index (),
    n      = fp -> Problem () -> nVars ();

  //assert (objInd >= 0);

  double *val = hessian -> val ();
  int    *row = hessian -> row ();
  int    *col = hessian -> col ();
  int     num = hessian -> num ();

  //printf ("compute square root:\n");

  // Remove objective's row and column (equivalent to taking the
  // Lagrangian's Hessian, setting f(x) = x_z = c, and recomputing the
  // hessian). 

  double maxElem = 0.; // used in adding diagonal element of x_z

  for (int i=num; i--; ++row, ++col, ++val) {

    //printf ("elem: %d, %d --> %g\n", *row, *col, *val);

    if ((*row == objInd) || 
	(*col == objInd))

      *val = 0;

    else if (fabs (*val) > maxElem)
      maxElem = fabs (*val);
  }

  val -= num;
  row -= num;
  col -= num;

  // fill an input to Lapack/Blas procedures using hessian

  double *A = (double *) malloc (n*n * sizeof (double));

  CoinZeroN (A, n*n);

  // Add Hessian part -- only lower triangular part
  for (int i=0; i<num; ++i, ++row, ++col, ++val)
    if (*col <= *row)
      A [*col * n + *row] = fp -> multHessMILP () * *val;

  // Add distance part
  for (int i=0; i<n; ++i)
    A [i * (n+1)] += fp -> multDistMILP ();

  // Add gradient-parallel term to the hessian, (x_z - x_z_0)^2. This
  // amounts to setting the diagonal element to GRADIENT_WEIGHT. Don't
  // do it directly on hessian

  if (objInd >= 0)
    A [objInd * (n+1)] = maxElem * GRADIENT_WEIGHT * n;

  // call Lapack/Blas routines
  double *eigenval = (double *) malloc (n   * sizeof (double));
  int status;

  // compute eigenvalues and eigenvectors
  Ipopt::IpLapackDsyev (true, n, A, n, eigenval, status);

  if      (status < 0) printf ("Couenne: warning, argument %d illegal\n",                     -status);
  else if (status > 0) printf ("Couenne: warning, dsyev did not converge (error code: %d)\n",  status);

  // define a new matrix B = E' * D, where E' is E transposed,
  //
  // E = eigenvector matrix
  // D = diagonal with square roots of eigenvalues
  //
  // Eventually, the square root is given by E' D E

  double *B = (double *) malloc (n*n * sizeof(double));

  double *eigenvec = A; // as overwritten by dsyev;

  // 
  // eigenvec is column major, hence post-multiplying it by D equals
  // multiplying each column i by the i-th eigenvalue
  //

  // if all eigenvalues are nonpositive, set them all to one

  double
    MinEigVal = eigenval [0],
    MaxEigVal = eigenval [n-1];

  for (int j=1; j<n; j++)
    assert (eigenval [j-1] <= eigenval [j]);

  if (MaxEigVal <= 0.)

    // in this case it makes sense to invert each eigenvalue
    // (i.e. take its inverse) and change its sign, as the steepest
    // descent should correspond to the thinnest direction

    for (int j=0; j<n; j++)
      eigenval [j] = 1. / (.1 - eigenval [j]);

  else {

    // set all not-too-positive ones to a fraction of the maximum
    // ("un-thins" the level curves defined by the HL)

    MinEigVal = MaxEigVal * COUENNE_EIG_RATIO;

    if (eigenval [0] <= MinEigVal) 
      for (int j=0; eigenval [j] <= MinEigVal; j++)
	eigenval [j] = MinEigVal;
  }

  // Now obtain sqrt (A)

  for (int j=0; j<n; ++j) {

    register double sqrtEig = sqrt (eigenval [j]);

    for (int i=n; i--;)
      *B++ = sqrtEig * eigenvec [i*n+j];
  }

  B -= n*n;

  // TODO: set up B as    row-major sparse matrix and 
  //              E as column-major sparse matrix
  //
  // Otherwise this multiplication is O(n^3)

  // Now compute B * E

  for     (int i=0; i<n; ++i)
    for   (int j=0; j<n; ++j) { 

      // multiply i-th row of B by j-th column of E

      double elem = 0.;

      for (int k=0; k<n; ++k)
	elem += B [i + k * n] * eigenvec [j * n + k];

      if (fabs (elem) > COUENNE_EPS)
	P [i]. insert (j, elem);
    }

  if (fp -> Problem () -> Jnlst () -> ProduceOutput (Ipopt::J_STRONGWARNING, J_NLPHEURISTIC)) {

    printf ("P:\n");


    printf ("P^{1/2}:\n");

  }

  free (eigenval);
  free (A);
  free (B);
}
