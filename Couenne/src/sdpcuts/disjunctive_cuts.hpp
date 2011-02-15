/* $Id$
 *
 * Name:    disjunctive_cuts.hpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef DISJUNCTIVE_CUTS_HPP
#define DISJUNCTIVE_CUTS_HPP


#include <CglCutGenerator.hpp>
#include <OsiSolverInterface.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>


#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))


void disjunctiveCutGen(const OsiSolverInterface &si, OsiCuts &cs, const double *sol, int n, Tracer *tracer);

#endif
