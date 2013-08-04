# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneDisjunctive.la

# List all source files for this library, including headers
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

# This is for libtool
libCouenneDisjunctive_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../heuristics` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	$(COUENNELIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
