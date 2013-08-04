# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneStandardize.la

# List all source files for this library, including headers
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

# This is for libtool
libCouenneStandardize_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../problem/depGraph` \
	-I`$(CYGPATH_W) $(srcdir)/../cut/sdpcuts` \
	$(COUENNELIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
