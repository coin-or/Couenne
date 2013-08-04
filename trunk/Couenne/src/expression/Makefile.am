# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneExpression.la

# List all source files for this library, including headers
libCouenneExpression_la_SOURCES = \
	domain.cpp \
	expression.cpp \
	exprOp.cpp \
	exprUnary.cpp \
	exprVar.cpp \
	exprAux.cpp \
	CouExpr.cpp \
	operators/exprAbs.cpp \
	operators/exprDiv.cpp \
	operators/exprExp.cpp \
	operators/exprInv.cpp \
	operators/exprLog.cpp \
	operators/exprMul.cpp \
	operators/exprTrilinear.cpp \
	operators/exprOpp.cpp \
	operators/exprPow.cpp \
	operators/exprSin.cpp \
	operators/exprCos.cpp \
	operators/exprSub.cpp \
	operators/exprSum.cpp \
	operators/exprMinMax.cpp \
	operators/exprGroup.cpp \
	operators/exprQuad.cpp \
	operators/compQuadFinBounds.cpp \
	operators/bounds/exprBQuad.cpp \
	partial/CouenneExprJac.cpp \
	partial/CouenneExprHess.cpp \
	../convex/operators/conv-exprAbs.cpp \
	../convex/operators/conv-exprDiv.cpp \
	../convex/operators/conv-exprMul.cpp \
	../convex/operators/conv-exprMul-genCuts.cpp \
	../convex/operators/conv-exprMul-reformulate.cpp \
	../convex/operators/conv-exprTrilinear.cpp \
	../convex/operators/conv-exprTrilinear-gencuts.cpp \
	../convex/operators/unifiedProdCuts.cpp \
	../convex/operators/exprMul-upperHull.cpp \
	../convex/operators/conv-exprOpp.cpp \
	../convex/operators/conv-exprPow.cpp \
	../convex/operators/conv-exprPow-getBounds.cpp \
	../convex/operators/conv-exprPow-envelope.cpp \
	../convex/operators/powNewton.cpp \
	../convex/operators/conv-exprSub.cpp \
	../convex/operators/conv-exprSum.cpp \
	../convex/operators/conv-exprInv.cpp \
	../convex/operators/conv-exprSinCos.cpp \
	../convex/operators/conv-exprExp.cpp \
	../convex/operators/conv-exprLog.cpp \
	../convex/operators/conv-exprGroup.cpp \
	../convex/operators/conv-exprQuad.cpp \
	../convex/operators/trigNewton.cpp \
	../convex/operators/alphaConvexify.cpp \
	../convex/operators/quadCuts.cpp \
	../branch/operators/branchExprAbs.cpp \
	../branch/operators/branchExprExp.cpp \
	../branch/operators/branchExprDiv.cpp \
	../branch/operators/branchExprInv.cpp \
	../branch/operators/branchExprLog.cpp \
	../branch/operators/branchExprMul.cpp \
	../branch/operators/branchExprTrilinear.cpp \
	../branch/operators/branchExprPow.cpp \
	../branch/operators/branchExprQuad.cpp \
	../branch/operators/branchExprSinCos.cpp \
	../branch/operators/minMaxDelta.cpp \
	../branch/operators/computeMulBrDist.cpp

# This is for libtool
libCouenneExpression_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../convex/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators/bounds` \
	-I`$(CYGPATH_W) $(srcdir)/../standardize` \
	-I`$(CYGPATH_W) $(srcdir)/../util` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../problem/depGraph` \
	$(COUENNELIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
