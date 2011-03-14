# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneProblem.la

# List all source files for this library, including headers
libCouenneProblem_la_SOURCES = \
	problem.cpp \
	CouenneProblem.cpp \
	problemIO.cpp \
	writeAMPL.cpp \
	writeGAMS.cpp \
	constraint.cpp \
	fillQuadIndices.cpp \
	fillDependence.cpp \
	checkNLP.cpp \
	getIntegerCandidate.cpp \
	testIntFix.cpp \
	CouenneSOS.cpp \
	CouenneProblemConstructors.cpp \
	CouenneRestoreUnused.cpp \
	reformulate.cpp \
	depGraph/depGraph.cpp \
	depGraph/checkCycles.cpp

# Note: the files below are #included in CouenneSolverInterface.hpp
# (it's the only way to make templates work -- see also discussions on
# the "export" keyword that, as of 2009, hasn't been implemented on
# most compilers)
#
# 	CouenneSolverInterface.cpp \
# 	CouenneLPtightenBounds.cpp \
# 	CouenneLPtightenBoundsCLP.cpp \
# 	CouenneLPtightenBoundsCLP-light.cpp 

if COIN_HAS_ASL
libCouenneProblem_la_SOURCES += \
	../readnl/readnl.cpp \
	../readnl/nl2e.cpp \
	../readnl/invmap.cpp \
	../readnl/readASLfg.cpp
endif

libCouenneProblem_la_DEPENDENCIES = $(libCouenneProblem_la_LIBADD)

# This is for libtool (on Windows)
libCouenneProblem_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/standardize` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem/depGraph` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem`
endif

if COIN_HAS_ASL
AM_CPPFLAGS += $(ASL_CPPFLAGS)
endif
