# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libCouenneReadnl.la

# List all source files for this library
libCouenneReadnl_la_SOURCES = CouenneAmplInterface.cpp
#libCouenneReadnl_la_SOURCES = \
#	readnl.cpp \
#	nl2e.cpp \
#	invmap.cpp \
#	readASLfg.cpp

# List all additionally required libraries
libCouenneReadnl_la_LIBADD = ../libCouenne.la

# Add libtool flags
LDFLAGS += $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../main \
	-I$(srcdir)/../interfaces \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../problem \
	$(COUENNELIB_CFLAGS)
