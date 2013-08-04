/* $Id$
 *
 * Name:    CouenneSparseMatrix.cpp
 * Authors: Pietro Belotti, Clemson University
 * Purpose: Implementation of a sparse Matrix for use in distance
 *          measurements in Feasibility Pump
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneSparseMatrix.hpp"
#include "CoinHelperFunctions.hpp"

using namespace Couenne;

/// Empty constructor
CouenneSparseMatrix::CouenneSparseMatrix ():

  num_  (0),
  val_  (NULL),
  col_  (NULL),
  row_  (NULL) {}


/// Empty constructor
CouenneSparseMatrix::~CouenneSparseMatrix () {

  if (val_) {

    free (val_);
    free (col_);
    free (row_);
  }
}

/// Copy constructor 
CouenneSparseMatrix::CouenneSparseMatrix (const CouenneSparseMatrix &rhs)
{operator= (rhs);}

/// Assignment 
CouenneSparseMatrix &CouenneSparseMatrix::operator= (const CouenneSparseMatrix &rhs) {

  num_ = rhs.num_;

  val_ = rhs.val_ && num_ ? CoinCopyOfArray (rhs.val_, num_) : NULL;
  col_ = rhs.col_ && num_ ? CoinCopyOfArray (rhs.col_, num_) : NULL;
  row_ = rhs.row_ && num_ ? CoinCopyOfArray (rhs.row_, num_) : NULL;

  return *this;
}

/// Clone
CouenneSparseMatrix *CouenneSparseMatrix::clone () 
{return new CouenneSparseMatrix (*this);}
