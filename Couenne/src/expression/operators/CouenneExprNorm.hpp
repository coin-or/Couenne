/* $Id$
 *
 * Name:    exprNorm.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the exprNorm class implementing $l_k$ norms
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef EXPRNORM_HPP
#define EXPRNORM_HPP

#include "CouenneExprOp.hpp"

namespace Couenne {

  /// Class for \f$ p \f$-norms, \f$ || f(x)||_p = \left(\sum_{i=1}^n f_i(x)^p\right)^{\frac{1}{p}} \f$

  class exprNorm: public exprOp {

  };

}

#endif
