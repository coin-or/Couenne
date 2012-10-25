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

namespace Couenne {

  // Base class for elements of our sparse structures ////////////////////////////

  class expression;

  class CouenneScalar {

  protected:

    int         index_; ///< index of element in vector
    expression *elem_;  ///< element

  public:

    CouenneScalar (int index, expression *elem):
      index_ (index),
      elem_  (elem) {}

    ~CouenneScalar () {}

    CouenneScalar (const CouenneScalar &rhs):
      index_ (rhs.index_),
      elem_  (rhs.elem_) {}

    CouenneScalar &operator= (const CouenneScalar &rhs) {
      index_ = rhs.index_;
      elem_  = rhs.elem_;
      return *this;
    }

    inline int         getIndex () const {return index_;}
    inline expression *getElem  () const {return elem_;}

    bool operator< (const CouenneScalar &rhs) const {return (index_ < rhs.index_);}
  };


  class CouenneSparseMatrix;

  // Sparse vector of expressions /////////////////////////////////////////////////

  class CouenneSparseVector {

  protected:

    std::set <CouenneScalar *> elem_;

  public:

    CouenneSparseVector () {}
   ~CouenneSparseVector () {}

    CouenneSparseVector (const CouenneSparseVector &rhs) {elem_ = rhs.elem_;}
    CouenneSparseVector *clone () {return new CouenneSparseVector (*this);}

    void add_element (int index, expression *elem);

    //void print ();

    const std::set <CouenneScalar *> &getElements () {return elem_;}

    double               operator *     (const CouenneSparseVector &factor)           const; ///< vector * vector (dot product)
    CouenneSparseVector &operator *     (const CouenneSparseMatrix &post)             const; ///< vector * matrix

    double               multiply_thres (const CouenneSparseVector &v2, double thres) const; ///< stops multiplication if above threshold
  };


  // Sparse matrix of expressions ///////////////////////////////////////////////////

  class CouenneSparseMatrix {

  protected:

    std::set <std::pair <int, CouenneSparseVector *> > row_;
    std::set <std::pair <int, CouenneSparseVector *> > col_;

  public:

    CouenneSparseMatrix () {}
   ~CouenneSparseMatrix () {}

    CouenneSparseMatrix (const CouenneSparseMatrix &rhs) {row_ = rhs.row_; col_ = rhs.col_;}
    CouenneSparseMatrix *clone () {return new CouenneSparseMatrix (*this);}

    const std::set <std::pair <int, CouenneSparseVector *> > &getRows    () const {return row_;}
    const std::set <std::pair <int, CouenneSparseVector *> > &getColumns () const {return col_;}

    void add_element (int row, int column, expression *elem);

    //void print ();

    CouenneSparseVector &operator * (const CouenneSparseVector &factor) const; ///< matrix * vector 

    //CouenneSparseMatrix &operator * (CouenneSparseMatrix &post);   ///< matrix * matrix
  };
}

#endif
