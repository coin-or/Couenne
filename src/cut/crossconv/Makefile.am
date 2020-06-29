# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneCrossConv.la

# List all source files for this library
libCouenneCrossConv_la_SOURCES = \
	CrossConvConstructors.cpp \
	CrossConvGenCuts.cpp \
	CrossConvSetup.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../.. \
	-I$(srcdir)/../../convex \
	-I$(srcdir)/../../heuristics \
	-I$(srcdir)/../../expression \
	-I$(srcdir)/../../expression/operators \
	-I$(srcdir)/../../main \
	-I$(srcdir)/../../branch \
	-I$(srcdir)/../../problem \
	$(COUENNELIB_CFLAGS)
