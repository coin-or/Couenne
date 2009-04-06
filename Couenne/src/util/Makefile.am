# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneUtil.la

# List all source files for this library, including headers
libCouenneUtil_la_SOURCES = \
	drawCuts.cpp \
	rootQ.cpp

libCouenneUtil_la_DEPENDENCIES = $(libCouenneUtil_la_LIBADD)

# This is for libtool (on Windows)
libCouenneUtil_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/util` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch`
endif
