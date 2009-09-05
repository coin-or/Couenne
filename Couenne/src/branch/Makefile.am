# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

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
	CouenneSOSObject.cpp
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

libCouenneBranch_la_DEPENDENCIES = $(libCouenneBranch_la_LIBADD)

# This is for libtool (on Windows)
libCouenneBranch_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/util` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/standardize` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch`
endif
