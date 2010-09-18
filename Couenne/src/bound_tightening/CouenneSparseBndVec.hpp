/* $Id$
 *
 * Name:    CouenneSparseBndVec.hpp
 * Author:  Pietro Belotti
 * Purpose: Keep track of tightened bounds with a sparse vector
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNESPARSEBNDVEC_HPP
#define COUENNESPARSEBNDVEC_HPP

namespace Couenne {

  template <class T> class CouenneSparseBndVec {

  private:

    /// maximum size
    int size_;

    /// current number of elements
    int n_; 

    /// indices vector, dense
    int *dInd_;

    /// indices vector, sparse (lots of garbage in between entries)
    int *sInd_;

    /// data vector, sparse (lots of garbage in between entries)
    T *data_;

  public:

    /// Constructor
    CouenneSparseBndVec (int size):

      size_ (size),
      n_    (0) {

      dInd_ = new int [size_];
      sInd_ = new int [size_];
      data_ = new T   [size_];
    }

    /// Copy constructor
    CouenneSparseBndVec (CouenneSparseBndVec &src):

      size_ (src.size_),
      n_    (src.n_) {

      for (int i=0; i<n_; i++) {

	int ind = (dInd_ [i] = src.dInd_ [i]);
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

    /// reset (eeeeasy!)
    void reset () 
    {n_ = 0;}

    /// access -- the only chance for garbage to be returned (and for
    /// valgrind to complain) is when object[ind] is READ without
    /// making sure it has been written
    T &operator[] (int index) {

      int sind = sInd_ [index];

      if ((sind < 0)   ||
	  (sind >= n_) || 
	  (dInd_ [sind] != index)) {
	
	dInd_ [n_] = index;
	sInd_ [index] = n_++;
      }

      return data_ [index];
    }

    /// return data -- use with care
    T *data ()
    {return data_;}

    /// return indices -- for use with data ()
    int *indices ()
    {return dInd_;}

    /// return current size
    int nElements ()
    {return n_;}
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
//     printf ("v [%d] = %d\n", v.indices () [i], v.data () [v.indices () [i]]);
// }
