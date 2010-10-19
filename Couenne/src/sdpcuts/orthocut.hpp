#ifndef ORTHOCUT_HPP
#define ORTHOCUT_HPP

#include <CglCutGenerator.hpp>
#include <OsiSolverInterface.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>



#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))





void orthoCutGen(const double *sol, int n, OsiCuts &cs, double *z, double *w, int m, Tracer *tracer);






#endif


