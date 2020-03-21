# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneBoundTightening.la

# List all source files for this library
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

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../problem \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../heuristics \
	-I$(srcdir)/../main \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/partial \
	-I$(srcdir)/../interfaces \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../expression/operators \
	$(COUENNELIB_CFLAGS)
