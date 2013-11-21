# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libCouenneReadnl.la

# List all source files for this library, including headers
libCouenneReadnl_la_SOURCES = CouenneAmplInterface.cpp
#libCouenneReadnl_la_SOURCES = \
#	readnl.cpp \
#	nl2e.cpp \
#	invmap.cpp \
#	readASLfg.cpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libCouenneReadnl_la_LIBADD = $(COUENNELIB_LIBS) ../libCouenne.la
endif

# This is for libtool
libCouenneReadnl_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../main` \
	-I`$(CYGPATH_W) $(srcdir)/../interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	$(COUENNELIB_CFLAGS) $(ASL_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
