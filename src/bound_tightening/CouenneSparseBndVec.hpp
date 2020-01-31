/* $Id$
 *
 * Name:    CouenneSparseBndVec.hpp
 * Author:  Pietro Belotti
 * Purpose: Keep track of tightened bounds with a sparse vector
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNESPARSEBNDVEC_HPP
#define COUENNESPARSEBNDVEC_HPP

namespace Couenne {

  template <class T> class CouenneSparseBndVec {

    /// Implements a fast sparse+dense vector data structure, with a
    /// size of n and a number k of nonzero elements. Usually, k<<n as
    /// happens in FBBT where k is the number of tightened variables
    /// and n is the number of variables. The main purpose is that of
    /// having a vector with
    ///
    /// 1) no O(n) initialization;
    ///
    /// 2) easy scan of the list of nonzero elements, i.e., O(k)
    /// rather than O(n).
    ///
    /// Implemented based on the (simple but beautiful) idea found at
    ///
    /// http://research.swtch.com/2008/03/using-uninitialized-memory-for-fun-and.html
    ///
    /// which in turn refers to a paper by Briggs and Torczon:
    ///
    /// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.30.7319
    ///
    /// NOTE: This will make valgrind complain on every assignment to
    /// non-previously-assigned entries. Get over it.

  private:

    unsigned int  size_; ///< maximum size
    unsigned int  n_;    ///< current number of elements
    unsigned int *dInd_; ///< indices vector, dense (garbage exists past n_)
    unsigned int *sInd_; ///< indices vector, sparse (lots of garbage in between entries)
    T            *data_; ///< data vector, sparse (lots of garbage in between entries)

  public:

    /// Constructor
    CouenneSparseBndVec (unsigned int size):

      size_ (size),
      n_    (0) {

      dInd_ = new unsigned int [size_];
      sInd_ = new unsigned int [size_];
      data_ = new T   [size_];
    }

    /// Copy constructor
    CouenneSparseBndVec (CouenneSparseBndVec &src):

      size_ (src.size_),
      n_    (src.n_) {

      for (register unsigned int i=0; i<n_; i++) {

	register unsigned int ind = (dInd_ [i] = src.dInd_ [i]);
	data_ [ind] = src.data_ [ind];	
	sInd_ [ind] = i; /// assert: src.sInd [ind] == i
      }
    }

    /// Destructor
    ~CouenneSparseBndVec () {
      delete [] sInd_;
      delete [] dInd_;
      delete [] data_;
    }

    /// Reset (eeeeasy!)
    void reset () 
    {n_ = 0;}

    /// Access -- the only chance for garbage to be returned (and for
    /// valgrind to complain) is when object[ind] is READ without
    /// making sure it has been written. This should not happen to the
    /// end user as read operations are only performed on the dense
    /// structure, after this object has been populated.
    T &operator[] (register unsigned int index) {

      register unsigned int &sind = sInd_ [index];

      if ((sind >= n_) || 
	  (dInd_ [sind] != index))
	dInd_ [sind = n_++] = index;	// this entry is new and has to be initialized	

      return data_ [sind];
    }

    /// Return data in DENSE format -- use with care
    T *data ()
    {return data_;}

    /// Return indices in DENSE format -- for use with data()
    unsigned int *indices ()
    {return dInd_;}

    /// Return current size
    unsigned int nElements ()
    {return n_;}

    /// Resize
    void resize (unsigned int newsize) 
    {size_ = newsize;}
  };
}

#endif

// #include <stdio.h>
// #include <stdlib.h>

// int main () {

//   Couenne::CouenneSparseBndVec <int> v (100);

//   v[84] = 10;
//   v[0]  = 60;
//   v[99] = 63;
//   v[72] = 16;
//   v[84] = 70;
//   v[25] = 33;
//   v[21] = 15;
//   v[21] = 12;
//   v[21] = 12;
//   v[8]  = 22;
//   v[4]  = 66;

//   srand48(1243235);
//   for (int i=0; i<1e9; i++)
//     v [(int)(99.999 * drand48())] = (int)(10000 * drand48());

//   for (int i=0; i< v.nElements(); i++)
//     printf ("v [%d] = %d\n", v.indices () [i], v.data () [i]);
// }
