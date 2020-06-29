# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneBranch.la

# List all source files for this library
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

if COIN_HAS_NAUTY
libCouenneBranch_la_SOURCES += CouenneNauty.cpp
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

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../problem \
	-I$(srcdir)/../util \
	-I$(srcdir)/../standardize \
	$(COUENNELIB_CFLAGS)
