# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneProblem.la

# List all source files for this library
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

if COIN_HAS_BONMINAMPL
libCouenneProblem_la_SOURCES += \
	../readnl/readnl.cpp \
	../readnl/nl2e.cpp \
	../readnl/invmap.cpp \
	../readnl/readASLfg.cpp \
	../readnl/CouenneAmplInterface.cpp
endif

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../standardize \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../problem/depGraph \
	-I$(srcdir)/../cut/sdpcuts \
	-I$(srcdir)/../main \
	-I$(srcdir)/../interfaces \
	$(COUENNELIB_CFLAGS)

