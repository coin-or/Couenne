# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneCrossConv.la

# List all source files for this library, including headers
libCouenneCrossConv_la_SOURCES = \
	CrossConvConstructors.cpp \
	CrossConvGenCuts.cpp \
	CrossConvSetup.cpp

# This is for libtool
libCouenneCrossConv_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../../heuristics` \
	-I`$(CYGPATH_W) $(srcdir)/../../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../../main` \
	-I`$(CYGPATH_W) $(srcdir)/../../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../../problem` \
	$(COUENNELIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src
