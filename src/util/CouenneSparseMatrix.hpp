/* $Id$
 *
 * Name:    CouenneSparseMatrix.hpp
 * Authors: Pietro Belotti, Clemson University
 * Purpose: Definition of a sparse Matrix for use in distance
 *          measurements in Feasibility Pump
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNESPARSEMATRIX_HPP
#define COUENNESPARSEMATRIX_HPP

namespace Couenne {

  /// Class for sparse Matrixs (used in modifying distances in FP)
  class CouenneSparseMatrix {

  public:

    /// Constructor 
    CouenneSparseMatrix ();

    /// Copy constructor 
    CouenneSparseMatrix (const CouenneSparseMatrix &);

    /// Assignment
    CouenneSparseMatrix &operator= (const CouenneSparseMatrix &rhs);

    /// Clone
    CouenneSparseMatrix *clone ();

    /// Destructor
    virtual ~CouenneSparseMatrix ();

    /// Get methods
    int     &num () {return num_;} ///< number of elements
    double *&val () {return val_;} ///< values
    int    *&col () {return col_;} ///< column indices
    int    *&row () {return row_;} ///< row indices

  private:

    /// Stores the values of the Matrix of the Lagrangian at optimum for later use
    int     num_; ///< number of elements
    double *val_; ///< values
    int    *col_; ///< column indices
    int    *row_; ///< row indices
  };
}

#endif
