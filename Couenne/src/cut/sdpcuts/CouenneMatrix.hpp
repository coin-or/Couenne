/* $Id$
 *
 * Name:    CouenneMatrix.hpp
 * Author:  Pietro Belotti
 * Purpose: define the class of expression matrices
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneMatrix_hpp
#define CouenneMatrix_hpp

#include <set>
#include <vector>

namespace Couenne {

  class expression;
  class CouenneExprMatrix;

  // Base class for elements of our sparse structures ////////////////////////////

  class CouenneScalar {

  protected:

    int         index_; ///< index of element in vector
    expression *elem_;  ///< element

  public:

    CouenneScalar (int index, expression *elem):
      index_ (index),
      elem_  (elem) {}

    ~CouenneScalar ();

    CouenneScalar (const CouenneScalar &rhs):
      index_ (rhs.index_),
      elem_  (rhs.elem_) {}

    CouenneScalar &operator= (const CouenneScalar &rhs) {
      index_ = rhs.index_;
      elem_  = rhs.elem_;
      return *this;
    }

    CouenneScalar *clone () {return new CouenneScalar (*this);}

    inline int         getIndex () const {return index_;}
    inline expression *getElem  () const {return elem_;}

    bool operator< (const CouenneScalar &rhs) const {return (index_ < rhs.index_);}

    void print () const;
  };


  // Sparse vector of expressions /////////////////////////////////////////////////

  class CouenneSparseVector {

  public:

    struct compare_scalars {
      inline bool operator() (register CouenneScalar * const &a, 
			      register CouenneScalar * const &b)
      {return a -> getIndex () < b -> getIndex ();}
    };

  protected:

    std::set <CouenneScalar *, compare_scalars> elem_;

  public:

    CouenneSparseVector () {}
   ~CouenneSparseVector ();

    CouenneSparseVector            (const CouenneSparseVector &rhs) {elem_ = rhs.elem_;}
    CouenneSparseVector &operator= (const CouenneSparseVector &rhs) {elem_ = rhs.elem_; return *this;}
    CouenneSparseVector *clone () {return new CouenneSparseVector (*this);}

    void add_element (int index, expression *elem);

    void print () const;

    const std::set <CouenneScalar *, compare_scalars> &getElements () {return elem_;}

    double               operator *     (const CouenneSparseVector &factor)           const; ///< vector * vector (dot product)
    CouenneSparseVector &operator *     (const CouenneExprMatrix &post)             const; ///< vector * matrix

    double               multiply_thres (const CouenneSparseVector &v2, double thres) const; ///< stops multiplication if above threshold
  };


  // Sparse matrix of expressions ///////////////////////////////////////////////////

  class CouenneExprMatrix {

  public:

    struct compare_pair_ind {
      inline bool operator() (register const std::pair <int, CouenneSparseVector *> &a, 
			      register const std::pair <int, CouenneSparseVector *> &b)
      {return a. first < b. first;}
    };

  protected:

    std::set <std::pair <int, CouenneSparseVector *>, compare_pair_ind> row_; ///< row major
    std::set <std::pair <int, CouenneSparseVector *>, compare_pair_ind> col_; ///< col major

    std::vector <expression *> varIndices_; ///< if used in sdp cuts, contains indices of x_i used in X_ij = x_i * x_j

  public:

    CouenneExprMatrix () {}
   ~CouenneExprMatrix ();

    CouenneExprMatrix (const CouenneExprMatrix &rhs);

    CouenneExprMatrix *clone () {return new CouenneExprMatrix (*this);}

    const std::set <std::pair <int, CouenneSparseVector *>, compare_pair_ind> &getRows () const {return row_;}
    const std::set <std::pair <int, CouenneSparseVector *>, compare_pair_ind> &getCols () const {return col_;}

    std::vector <expression *> &varIndices () {return varIndices_;}

    void add_element (int row, int column, expression *elem);
    void print () const;
    long unsigned int size ();

    CouenneSparseVector &operator * (const CouenneSparseVector &factor) const; ///< matrix * vector 
    CouenneExprMatrix   &operator * (const CouenneExprMatrix   &post)   const; ///< matrix * matrix
  };
}

#endif
