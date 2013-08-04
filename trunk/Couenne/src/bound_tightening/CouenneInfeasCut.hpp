/* $Id$
 *
 * Name:    CouenneInfeasCut.hpp
 * Author:  Pietro Belotti
 * Purpose: An infeasible cut to tell the node solver this node is infeasible
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEINFEASCUT_HPP
#define COUENNEINFEASCUT_HPP

#include "OsiCuts.hpp"

/// Add a fictitious cut 1<= x_0 <= -1 as a signal to the node solver
/// that this node is deemed infeasible by this cut generator (most
/// likely a bound tightener).

void WipeMakeInfeas (OsiCuts &cs);


/// Check whether the previous cut generators have added an infeasible
/// cut.

bool isWiped (OsiCuts &cs);

#endif
