# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneHeuristics.la

# List all source files for this library
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

# Here list all include flags, relative to this "srcdir" directory
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../expression/partial \
	-I$(srcdir)/../interfaces \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../util \
	-I$(srcdir)/../problem \
	$(COUENNELIB_CFLAGS)
