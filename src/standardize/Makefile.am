# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneStandardize.la

# List all source files for this library
libCouenneStandardize_la_SOURCES = \
	decomposeTerm.cpp \
	linStandardize.cpp \
	sumStandardize.cpp \
	analyzeSparsity.cpp \
	constrStandardize.cpp \
	elementBreak.cpp \
	splitAux.cpp \
	auxiliarize.cpp \
	standardize.cpp \
	flattenMul.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../problem \
	-I$(srcdir)/../problem/depGraph \
	-I$(srcdir)/../cut/sdpcuts \
	$(COUENNELIB_CFLAGS) -DCOUENNE_BUILD
