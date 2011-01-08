# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Andreas Waechter            IBM    2006-12-06
#           Pietro Belotti

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Couenne],[trunk],[couenne@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 Carnegie Mellon University, International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/expression/expression.cpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# Initialization specifically for project main directories
AC_COIN_PROJECTDIR_INIT(Couenne)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Couenne)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# for readnl.cpp
AC_COIN_CHECK_CXX_CHEADER(stdint)

#############################################################################
#                            COIN-OR components                             #
#############################################################################

AC_COIN_CHECK_PACKAGE(Bonmin, [bonmin], [CouenneLib])
if test $coin_has_bonmin != yes ; then
  AC_MSG_ERROR([Required package Bonmin not available.])
fi

#############################################################################
#                       Other third party software                          #
#############################################################################

# Bonmin AMPL interface
AC_COIN_CHECK_PACKAGE(ASL, [bonminamplinterface])

# Nauty library (for symmetry detection)
AC_LANG_PUSH(C)
AC_COIN_CHECK_USER_LIBRARY([Nauty],[NTY],[nauty.h],[nauty],,[CouenneLib])
AC_LANG_POP(C)

# CPLEX: check for OsiCplex, and check for cplex include dir, so we can include cplex.h
AC_COIN_CHECK_PACKAGE(Cpx, [osi-cplex], [CouenneLib])
AC_ARG_WITH([cplex-incdir],
  [AS_HELP_STRING([cplex-incdir], [specify the header file directory for library cplex])],
  [CPXINCDIR=`cd $withval; pwd`])

# other LP solvers
AC_COIN_CHECK_PACKAGE(Grb, [osi-gurobi], [CouenneLib])  # Osi/Gurobi
AC_COIN_CHECK_PACKAGE(Msk, [osi-mosek],  [CouenneLib])  # Osi/Mosek
AC_COIN_CHECK_PACKAGE(Spx, [osi-soplex], [CouenneLib])  # Osi/Soplex
AC_COIN_CHECK_PACKAGE(Xpr, [osi-xpress], [CouenneLib])  # Osi/XpressMP

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 test/Makefile
                 src/main/Makefile
                 src/bound_tightening/Makefile
                 src/branch/Makefile
                 src/convex/Makefile
                 src/interfaces/Makefile
                 src/heuristics/Makefile
                 src/disjunctive/Makefile
                 src/sdpcuts/Makefile
                 src/ellipcuts/Makefile
                 src/expression/Makefile
                 src/problem/Makefile
                 src/standardize/Makefile
                 src/readnl/Makefile
                 src/util/Makefile
                 src/crossconv/Makefile
                 src/two_implied_bt/Makefile
                 src/Makefile
		 test/grandTest.sh
		 couenne.pc couenne-uninstalled.pc])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_couenne.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
