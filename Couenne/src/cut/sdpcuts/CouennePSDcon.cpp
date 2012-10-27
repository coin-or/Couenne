/* $Id$
 *
 * Name:    CouennePSDcon.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation, positive semidefinite constraints
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneMatrix.hpp"
#include "CouennePSDcon.hpp"

using namespace Couenne;

/// Destructor
CouennePSDcon::~CouennePSDcon ()
{if (X_) delete X_;}

/// Copy constructor
CouennePSDcon::CouennePSDcon (const CouennePSDcon &c, Domain *d)
{X_ = c.X_ -> clone ();}

/// decompose body of constraint through auxiliary variables
exprAux *CouennePSDcon::standardize (CouenneProblem *p)
{return NULL;}

/// print constraint
void CouennePSDcon::print (std::ostream &os) {}
