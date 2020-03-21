# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

noinst_LTLIBRARIES = libCouenneInterfaces.la

# List all source files for this library
libCouenneInterfaces_la_SOURCES = \
	BonCouenneInterface.cpp \
	CouenneTNLP.cpp \
	CouenneMINLPInterface.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../expression/partial \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../util \
	-I$(srcdir)/../problem \
	$(COUENNELIB_CFLAGS)
