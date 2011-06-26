/* $Id$
 *
 * Name:    CouenneExprJac.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Jacobian expression
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneExprJac_HPP
#define CouenneExprJac_HPP

namespace Couenne {

  class expression;
  class CouenneProblem;

  /// Jacobian of the problem (computed through Couenne expression
  /// classes).

  class ExprJac {

  private:

    int          nnz_;   ///< number of (symbolic) nonzeroes
    int         *iRow_;  ///< row indices (read this way by eval_jac_g)
    int         *jCol_;  ///< col indices

    expression **expr_;  ///< nonzero expression elements (there are nnz_ of them)

    int          nRows_; ///< number of actual constraints

  public:

    ExprJac  ();
    ExprJac  (CouenneProblem *);
    ~ExprJac ();

    ExprJac  (const ExprJac &);
    ExprJac *clone ();
    ExprJac &operator= (const ExprJac &);

    int  nnz  () const {return nnz_;}
    int *iRow () const {return iRow_;}
    int *jCol () const {return jCol_;}

    expression **expr () const {return expr_;}

    int nRows () const {return nRows_;}
  };
}

#endif
