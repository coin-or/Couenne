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

  /// expression Jacobian
  class ExprJac {

  private:

    int          nnz_;  ///< number of (symbolic) nonzeroes
    int         *iRow_; ///< col indices
    int         *jCol_; ///< row starts
    expression **expr_; ///< expression

    int          nRows_; ///< number of actual constraints

  public:

    ExprJac  ();
    ExprJac  (CouenneProblem *);
    ~ExprJac ();

    int  nnz  () const {return nnz_;}
    int *iRow () const {return iRow_;}
    int *jCol () const {return jCol_;}

    expression **expr () const {return expr_;}

    int nRows () const {return nRows_;}
  };
}

#endif
