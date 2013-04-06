/* $Id$
 *
 * Name:    CrossConvGenCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: Convexification cuts on redundant relationships between auxiliaries
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCrossConv.hpp"

using namespace Couenne;

/// the main CglCutGenerator
void CouenneCrossConv::generateCuts (const OsiSolverInterface &, 
				     OsiCuts &, 
				     const CglTreeInfo)
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
   const
#endif
{

}
