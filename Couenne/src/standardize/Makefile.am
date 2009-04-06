# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

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

# This is for libtool (on Windows)
libCouenneStandardize_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem/depGraph` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` 
endif
