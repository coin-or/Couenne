# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneConvex.la

# List all source files for this library
libCouenneConvex_la_SOURCES = \
	addEnvelope.cpp \
	createCuts.cpp \
	CouenneCutGenerator.cpp \
	generateCuts.cpp \
	genColCuts.cpp \
	genRowCuts.cpp \
	updateBranchInfo.cpp \
	isOptimumCut.cpp

# 	operators/conv-exprAbs.cpp \
# 	operators/conv-exprDiv.cpp \
# 	operators/conv-exprMul.cpp \
# 	operators/conv-exprMul-genCuts.cpp \
# 	operators/unifiedProdCuts.cpp \
# 	operators/conv-exprOpp.cpp \
# 	operators/conv-exprPow.cpp \
# 	operators/conv-exprPow-getBounds.cpp \
# 	operators/conv-exprPow-envelope.cpp \
# 	operators/powNewton.cpp \
# 	operators/conv-exprSub.cpp \
# 	operators/conv-exprSum.cpp \
# 	operators/conv-exprInv.cpp \
# 	operators/conv-exprSinCos.cpp \
# 	operators/conv-exprExp.cpp \
# 	operators/conv-exprLog.cpp \
# 	operators/conv-exprGroup.cpp \
# 	operators/conv-exprQuad.cpp \
# 	operators/trigNewton.cpp \
# 	operators/alphaConvexify.cpp \
# 	operators/quadCuts.cpp 

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../expression/operators/bounds \
	-I$(srcdir)/../util \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../problem \
	$(COUENNELIB_CFLAGS)
