# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

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
	writeLP.cpp \
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
	CouenneGlobalCutOff.cpp \
	CouenneSymmetry.cpp \
	CouenneRecordBestSol.cpp \
	depGraph/depGraph.cpp \
	depGraph/checkCycles.cpp 

# Note: the files below are #included in CouenneSolverInterface.hpp
# (it's the only way to make templates work -- see also discussions on
# C++'s "export" keyword that, as of 2009, hasn't been implemented on
# most compilers).
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

# This is for libtool
libCouenneProblem_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../standardize` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../problem/depGraph` \
	-I`$(CYGPATH_W) $(srcdir)/../cut/sdpcuts` \
	$(COUENNELIB_CFLAGS) $(ASL_CFLAGS)

if COIN_HAS_NTY
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)/`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
