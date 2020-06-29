# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Andreas Waechter            IBM    2006-12-06
#           Pietro Belotti

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Couenne],[devel],[https://github.com/coin-or/Couenne/issues/new],[coin-or-couenne],[https://github.com/coin-or/Couenne])

AC_COPYRIGHT([
Copyright 2006 Carnegie Mellon University, International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/expression/expression.cpp)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize libtool
AC_COIN_PROG_LIBTOOL

#############################################################################
#                            COIN-OR components                             #
#############################################################################

AC_COIN_CHK_PKG(Bonmin,[CouenneLib])
if test $coin_has_bonmin != yes ; then
  AC_MSG_ERROR([Required package Bonmin not available.])
fi

AC_COIN_CHK_PKG(ASL,[CouenneLib],[bonminamplinterface])

# more workaround for having ASL deep inside Couenne
if test $coin_has_asl = yes ; then
  AC_DEFINE(COUENNEINTERFACE_FROM_ASL, 1, [whether CouenneInterface is derived from Bonmins AmplInterface])
fi

# Osi interfaces to LP solvers
AC_COIN_CHK_PKG(OsiCpx,[CouenneLib],[osi-cplex])  # Osi/Cplex
AC_COIN_CHK_PKG(OsiGrb,[CouenneLib],[osi-gurobi]) # Osi/Gurobi
AC_COIN_CHK_PKG(OsiMsk,[CouenneLib],[osi-mosek])  # Osi/Mosek
AC_COIN_CHK_PKG(OsiSpx,[CouenneLib],[osi-soplex]) # Osi/SoPlex
AC_COIN_CHK_PKG(OsiXpr,[CouenneLib],[osi-xpress]) # Osi/Xpress

#############################################################################
#                       Other third party software                          #
#############################################################################

# Nauty library (for symmetry detection)
AC_COIN_CHK_PKG([Nauty],[CouenneLib],[nauty])

AC_LANG_PUSH(C)

# SCIP (used in Feasibility Pump)
AC_COIN_CHK_LIB([SCIP],[CouenneLib],[-lscipsolver],[],[],[SCIPcreate])

# CPLEX (so we can include cplex.h)
AC_COIN_CHK_LIB([Cplex],[CouenneLib],[-lcplex -lpthread -lm -ldl],[],[],[CPXgetstat])

AC_LANG_POP(C)

# for dsyevx in sdpcuts
AC_COIN_CHK_LAPACK(CouenneLib)

#############################################################################
#                               System header                               #
#############################################################################

AC_LANG_PUSH(C++)
# for readnl
AC_CHECK_HEADERS([cstdint stdint.h],[break],[])
AC_LANG_POP(C++)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Cbc Ipopt Bonmin)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

AC_COIN_FINALIZE_FLAGS([CouenneLib])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 test/Makefile
                 src/main/Makefile
                 src/bound_tightening/Makefile
                 src/bound_tightening/twoImpliedBT/Makefile
                 src/branch/Makefile
                 src/convex/Makefile
                 src/interfaces/Makefile
                 src/heuristics/Makefile
                 src/disjunctive/Makefile
                 src/cut/crossconv/Makefile
                 src/cut/sdpcuts/Makefile
                 src/expression/Makefile
                 src/problem/Makefile
                 src/standardize/Makefile
                 src/util/Makefile
                 src/Makefile
                 test/grandTest.sh
                 couenne.pc 
                 doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_couenne.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
