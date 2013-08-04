# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

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

# This is for libtool
libCouenneConvex_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators/bounds` \
	-I`$(CYGPATH_W) $(srcdir)/../util` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	$(COUENNELIB_CFLAGS)

if COIN_HAS_NTY
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)/`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
