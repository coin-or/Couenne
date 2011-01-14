/* $Id$
 *
 * Name:    CouenneExprHess.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Hessian of the Lagrangian, definition
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneExprHess_HPP
#define CouenneExprHess_HPP

namespace Couenne {

  class expression;
  class CouenneProblem;

  /// expression matrices
  class ExprHess {

  private:

    int           nnz_;  ///< number of (symbolic) nonzeroes
    int          *iRow_; ///< col indices
    int          *jCol_; ///< row starts

    /// there are m+1 (m constraints + 1 obj).
    ///
    /// Implementing a FP requires adding one for gg', the gradient
    /// again being set up at the beginning (at least its expression
    /// members are known).
    ///
    /// This can simply be hacked by the FP itself. Same for the
    /// changed hessian, simply replace the CouenneProblem's objective
    /// with sum (objective, norm)
    ///
    /// Actually, we could do the gg' trick by replacing the objective
    /// with sum (objective, norm, gg')

    int          *numL_; ///< size of each lamI_

    int         **lamI_; ///< vector of indices in the lambda vector
		 	 ///< whose constraint has nonzero entry in
		 	 ///< this position of the hessian

    expression ***expr_; ///< list of lists of pointers to expression

  public:

    ExprHess  ();
    ExprHess  (CouenneProblem *);
    ~ExprHess ();

    int   nnz  () {return nnz_;}
    int  *iRow () {return iRow_;}
    int  *jCol () {return jCol_;}
    int  *numL () {return numL_;}
    int **lamI () {return lamI_;}

    expression ***expr () {return expr_;}
  };
}

#endif
