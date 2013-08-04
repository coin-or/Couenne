/** $Id$
 *
 * @file   cons_rowcuts.h
 * @brief  constraint handler for rowcuts constraints
 *         enables separation of convexification cuts during SCIP solution procedure
 * @author Pietro Belotti
 * @author Timo Berthold
 * @license This file is licensed under the Eclipse Public License (EPL)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ROWCUTS_H__
#define __SCIP_CONS_ROWCUTS_H__

#include "CouenneCutGenerator.hpp"
#include "OsiSolverInterface.hpp"

#ifdef COIN_HAS_SCIP

#include "scip/scip.h"

using namespace Couenne;

/** creates the handler for rowcuts constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrRowcuts(
   SCIP*                 scip,               /**< SCIP data structure */
   CouenneCutGenerator*  cutgenerator,       /**< CouenneCutGenerator for linearization cuts */
   OsiSolverInterface*   milp                /**< Couenne's MILP relaxation of Couenne's MINLP */
   );

#endif
#endif
