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
AC_COIN_PROJECTDIR_INIT

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
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(ThirdParty)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(Bonmin)
AC_COIN_HAS_PROJECT(Ipopt)
AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Couenne)

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_HAS_ASL

# Nauty library (for symmetry detection)
AC_LANG_PUSH(C)
AC_COIN_HAS_USER_LIBRARY([Nauty],[NTY],[nauty.h],[nauty])
AC_LANG_POP(C)

AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])         # Cplex library
AC_COIN_HAS_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBgetattrinfo]) # Gurobi library
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_openmosek])      # Mosek
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])      # OSL
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h])                 # Soplex library
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])  # XpressMP

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

#AC_COIN_VPATH_LINK(examples/hello.mps)
#AC_COIN_VPATH_LINK(examples/input.130)

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
                 src/Makefile
		 test/grandTest.sh
		 couenne_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_couenne.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
