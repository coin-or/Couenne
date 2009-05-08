# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneConvex.la

# List all source files for this library, including headers
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

libCouenneConvex_la_DEPENDENCIES = $(libCouenneConvex_la_LIBADD)

# This is for libtool (on Windows)
libCouenneConvex_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators/bounds` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/util` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem`
endif
