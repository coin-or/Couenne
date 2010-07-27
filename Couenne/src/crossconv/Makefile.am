# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneCrossConv.la

# List all source files for this library, including headers
libCouenneCrossConv_la_SOURCES = \
	CrossConvConstructors.cpp \
	CrossConvGenCuts.cpp \
	CrossConvSetup.cpp

libCouenneCrossConv_la_DEPENDENCIES = $(libCouenneCrossConv_la_LIBADD)

# This is for libtool (on Windows)
libCouenneCrossConv_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src`  \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src`

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/heuristics` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/main` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` 
endif

if COIN_HAS_NTY
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)`
endif
