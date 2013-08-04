# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneBoundTightening.la

# List all source files for this library, including headers
libCouenneBoundTightening_la_SOURCES = \
	boundTightening.cpp \
	aggressiveBT.cpp \
	fake_tightening.cpp \
	obbt.cpp \
	obbt_iter.cpp \
	tightenBounds.cpp \
	impliedBounds.cpp \
	FixPointConstructors.cpp \
	FixPointGenCuts.cpp \
	CouenneInfeasCut.cpp \
	obbt_supplement.cpp \
	BTPerfIndicator.cpp \
	BTPerfIndicatorConstr.cpp \
	operators/impliedBounds-mul.cpp \
	operators/impliedBounds-sum.cpp \
	operators/impliedBounds-exprSum.cpp \
	operators/impliedBounds-exprDiv.cpp \
	operators/impliedBounds-exprMul.cpp \
	operators/impliedBounds-exprTrilinear.cpp \
	operators/impliedBounds-exprQuad.cpp \
	operators/impliedBounds-exprPow.cpp \
	../main/BonCouenneInfo.cpp \
	CouenneAggrProbing.cpp 

# This is for libtool
libCouenneBoundTightening_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../heuristics` \
	-I`$(CYGPATH_W) $(srcdir)/../main` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/partial` \
	-I`$(CYGPATH_W) $(srcdir)/../interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	$(COUENNELIB_CFLAGS) 

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
