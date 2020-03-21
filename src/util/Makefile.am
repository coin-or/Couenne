# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneUtil.la

# List all source files for this library, including headers
libCouenneUtil_la_SOURCES = \
	drawCuts.cpp \
	rootQ.cpp \
	CouenneSparseMatrix.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../problem \
	-I$(srcdir)/../util \
	$(COUENNELIB_CFLAGS)
