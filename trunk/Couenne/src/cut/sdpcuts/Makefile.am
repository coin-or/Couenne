# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneSdpCuts.la

# List all source files for this library, including headers

libCouenneSdpCuts_la_SOURCES = \
	CouenneSdpCuts.cpp \
	CouennePSDcon.cpp \
	CouenneMatrix.cpp \
	CutGen.cpp \
	CutGenSparse.cpp \
	dsyevx_wrapper.cpp

# This is for libtool
libCouenneSdpCuts_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../..` \
	-I`$(CYGPATH_W) $(srcdir)/../../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../../problem` \
	$(COUENNELIB_CFLAGS)

if COIN_HAS_CPX
AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src
