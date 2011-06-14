/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_rowcuts.h,v 1.27.2.1 2011/01/02 11:19:45 bzfheinz Exp $"

/**@file   cons_rowcuts.h
 * @brief  constraint handler for rowcuts constraints
 * @author Pietro Belotti
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ROWCUTS_H__
#define __SCIP_CONS_ROWCUTS_H__

#include "CouenneCutGenerator.hpp"
#include "OsiSolverInterface.hpp"

#include "scip/scip.h"
namespace Couenne
{
/** creates the handler for rowcuts constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrRowcuts(
   SCIP*                 scip,               /**< SCIP data structure */
   CouenneCutGenerator*  cutgenerator,       /**< CouenneCutGenerator for linearization cuts */
   OsiSolverInterface*   milp                /**< Couenne's MILP relaxation of Couenne's MINLP */
   );
}
#endif
