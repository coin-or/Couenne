# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneDisjunctive.la

# List all source files for this library
libCouenneDisjunctive_la_SOURCES = \
        CouenneDisjCuts.cpp \
        generateDisjCuts.cpp \
        disjCut.cpp \
        getDisjunctions.cpp \
        separateWithDisjunction.cpp \
        singleDisjunctions.cpp \
        dummyVecMatr.cpp \
        OsiLinear2MatrVec.cpp \
        OsiCuts2MatrVec.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../heuristics \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../problem \
	$(COUENNELIB_CFLAGS)
