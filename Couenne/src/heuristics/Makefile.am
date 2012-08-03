# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneHeuristics.la

# List all source files for this library, including headers
libCouenneHeuristics_la_SOURCES = \
	BonInitHeuristic.cpp \
	BonNlpHeuristic.cpp \
	CouenneFeasPumpConstructors.cpp \
	CouenneFeasPump.cpp \
	CouenneFPcreateMILP.cpp \
	CouenneFPSolveMILP.cpp \
	CouenneFPSolveNLP.cpp \
	CouenneFPFindSolution.cpp \
	CouenneFPpool.cpp \
	CouenneFPscipSolve.cpp \
	CouenneIterativeRounding.cpp 

if COIN_HAS_SCIP
libCouenneHeuristics_la_SOURCES += cons_rowcuts.cpp
endif

# This is for libtool
libCouenneHeuristics_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/partial` \
	-I`$(CYGPATH_W) $(srcdir)/../interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../util` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	$(COUENNELIB_CFLAGS)

if COIN_HAS_CPX
AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
