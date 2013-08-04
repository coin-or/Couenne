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

  /// expression matrices. Used to evaluate the Hessian of the
  /// Lagrangian function at an optimal solution of the NLP

  class ExprHess {

  private:

    int   nnz_;  ///< number of (symbolic) nonzeroes
    int  *iRow_; ///< row indices (read this way by eval_h)
    int  *jCol_; ///< col indices

    /// There are m+1 (m constraints + 1 obj) components:
    ///
    /// \f$ \nabla^2 \mathcal L (x,\lambda) = \nabla^2 f(x) + \lambda^\top \nabla^2 g(x) \f$
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

    int  *numL_; ///< size of each lamI_

    int **lamI_; ///< vector of indices in the lambda vector
                 ///< whose constraint has nonzero entry in
                 ///< this position of the hessian

    expression ***expr_; ///< list of lists of pointers to expression

  public:

    ExprHess  ();
    ExprHess  (CouenneProblem *);

    ExprHess  (const ExprHess &);
    ExprHess  &operator=(const ExprHess &);
    ExprHess  *clone ();

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
