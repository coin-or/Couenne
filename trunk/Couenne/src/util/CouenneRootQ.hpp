/* $Id$
 *
 * Name:    rootQ.hpp
 * Author:  Pietro Belotti
 * Purpose: find roots of polynomial Q^k(x) (see Liberti and Pantelides, 2003)
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneRootQ_hpp
#define CouenneRootQ_hpp

#include "CouenneTypes.hpp"
#include <map>
#include <iostream>

namespace Couenne
{

/// Find roots of polynomial $Q^k(x) = \f$\sum_{i=1}^{2k} i x^{i-1}\f$. Used
/// in convexification of powers with odd exponent
CouNumber rootQ (int k);


/// class that stores result of previous calls to rootQ into a map for
/// faster access

class Qroot {

 protected:

  /// Maps an integer k with the root of \f$Q^k(x)\f$.

  static std::map <int, CouNumber> Qmap;

 public:

  /// Empty constructor -- we only need the method to work on the static
  /// structure.
  Qroot  () {}

  /// Empty destructor.
  ~Qroot () {}

  /// Retrieve root of Q with order = k. If no such computation has
  /// been performed yet, do it here

  inline CouNumber operator () (int k) {

    std::map <int, CouNumber>:: iterator pos;
    CouNumber root;

    if( k % 2 == 1 )
    {
       k /= 2; // becomes true index

       if ((pos = Qmap.find (k)) == Qmap.end()) {

          std::pair <int, CouNumber> newpair;

          newpair.first  = k;
          newpair.second = (root = rootQ (k));

          Qmap.insert (newpair);
       }
       else
          root = pos -> second;
    }
    else
    {
       switch(k) {
          case 2: return -(sqrt(2.0) - 1.0);
          case 4: return -0.56042566045031785945;
          case 6: return -0.64146546982884663257;
          case 8: return -0.69428385661425826738;
          case 10: return -0.73192937842370733350;
          default:
             std::cerr << "Need to implement root finding for even k" << std::endl;
             throw -1;
       }
    }

    return root;
  }
};

}

#endif
