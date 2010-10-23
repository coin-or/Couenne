#ifndef RLT_CUTS_HPP
#define RLT_CUTS_HPP

#include <CglCutGenerator.hpp>
#include <tracer.hpp>


void rltCutsGen(const double *sol, int n, OsiCuts &cs, double *lb, double *ub, int m, Tracer *tracer);

#endif
