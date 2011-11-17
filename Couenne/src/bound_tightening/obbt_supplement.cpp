/* $Id$
 *
 * Name:    obbt_supplement.cpp
 * Author:  Pietro Belotti
 * Purpose: Post-OBBT computations that use dual information to infer
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiSolverInterface.hpp"

// Use dual information lambda to obtain, from solution to this
// problem, a dual bound to OBBT subproblem (min|max) for every
// other variable.

int obbt_supplement (const OsiSolverInterface *csi, /// interface to use as a solver
		     int index,                     /// variable being looked at
		     int sense) {                   /// 1: minimize, -1: maximize

  if (!(csi -> isProvenOptimal ()))
    return 0;

  int
    n = csi -> getNumCols (),
    m = csi -> getNumRows ();

  // get dual

  const double *lambda = csi -> getRowPrice ();

  double alpha_i = (sense==1 ? 1. : -1.);

  double *beta = new double [m];

  // The problem just solved was either (depending on sense):
  //
  // LB) min {  x_i: Ax>=b, l<=x<=u}
  // UB) min { -x_i: Ax>=b, l<=x<=u}
  //
  // lambda contains the dual variables at the optimal solution,
  // i.e., the Lagrangian multipliers of the problems
  //
  // L_lb (lambda) = min { x_i + lambda^T (b-Ax): l<=x<=u}
  // U_lb (lambda) = min {-x_i + lambda^T (b-Ax): l<=x<=u}
  //
  // Suppose M={1..m} and N={1..n} index the set of rows and columns,
  // respectively. Rewrite both problems above by setting alpha_i as 1
  // for the LB problem and -1 for the UB problem.
  //
  // L (i, lambda) = min {alpha_i x_i + sum {h in M} lambda_h (b_h - sum {k in N} a_hk x_k):        l <= x <= u}
  //               = lambda^T b + min {alpha_i x_i - sum {h in M} sum {k in N}  lambda_h a_hk x_k:  l <= x <= u}
  //               = lambda^T b + min {alpha_i x_i - sum {k in N} (sum {h in M} lambda_h a_hk) x_k: l <= x <= u}
  //               = lambda^T b + min {alpha_i x_i - sum {k in N} beta_k x_k:                       l <= x <= u}
  //               = lambda^T b + min { sum {k in N} gamma_i_k x_k:                                 l <= x <= u}
  //
  //               = lambda^T b + sum {k in N: gamma_i_k > 0} gamma_i_k l_k + 
  //                              sum {k in N: gamma_i_k < 0} gamma_i_k u_k,
  //
  // where
  //
  // beta_k = sum {h in M} lambda_h a_hk and 
  //
  // and
  //
  // gamma_i_k = - beta_i + alpha_i   if i==k
  //             - beta_k             otherwise.
  //
  // Then any dual solution to the i-th LB or UB problems above can be
  // used to get a DUAL solution to all other j-th LB/UB problems, for
  // j!=i. Just compute L (j,lambda) for all j, for alpha_i in {-1,1},
  // and compare with previous bound.

  for (int j=0; j<n; j++) {

  }
}
