# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneBranch.la

# List all source files for this library, including headers
libCouenneBranch_la_SOURCES = \
	CouenneThreeWayBranchObj.cpp \
	CouenneBranchingObject.cpp \
	CouenneObject.cpp \
	CouenneVarObject.cpp \
	CouenneChooseVariable.cpp \
	CouenneChooseStrong.cpp \
	doStrongBranching.cpp \
	projections.cpp \
	infeasibility.cpp \
	infeasibilityVT.cpp \
	feasibleRegion.cpp \
	CouenneComplObject.cpp \
	CouenneComplBranchingObject.cpp \
	CouenneOrbitObj.cpp \
	CouenneOrbitBranchingObj.cpp \
	BranchCore.cpp \
	CouenneSOSObject.cpp \
	StrongBranchingSetupList.cpp

if COIN_HAS_NTY
libCouenneBranch_la_SOURCES += Nauty.cpp
endif

# 	operators/branchExprAbs.cpp \
# 	operators/branchExprExp.cpp \
# 	operators/branchExprDiv.cpp \
# 	operators/branchExprInv.cpp \
# 	operators/branchExprLog.cpp \
# 	operators/branchExprMul.cpp \
# 	operators/branchExprPow.cpp \
# 	operators/branchExprQuad.cpp \
# 	operators/branchExprSinCos.cpp \
# 	operators/minMaxDelta.cpp \
# 	operators/computeMulBrDist.cpp

# This is for libtool
libCouenneBranch_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../util` \
	-I`$(CYGPATH_W) $(srcdir)/../standardize` \
	$(COUENNELIB_CFLAGS)

if COIN_HAS_NTY
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)/`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
