# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneTwoImplied.la

# List all source files for this library
libCouenneTwoImplied_la_SOURCES = \
	TwoImpliedConstructors.cpp \
	TwoImpliedIneqCombine.cpp \
	TwoImpliedGenCuts.cpp 

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../.. \
	-I$(srcdir)/../../expression \
	-I$(srcdir)/../../branch \
	-I$(srcdir)/../../problem \
	-I$(srcdir)/../../bound_tightening \
	$(COUENNELIB_CFLAGS)
