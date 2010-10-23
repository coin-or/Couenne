# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneSdpCuts.la

# List all source files for this library, including headers
libCouenneSdpCuts_la_SOURCES =

# libCouenneSdpCuts_la_SOURCES = \
# 	CutGen.cpp \
# 	dsyevx_wrapper.cpp \
# 	Heuristics.cpp \
# 	linquad_cuts.cpp \
# 	misc_util.cpp \
# 	orthocut.cpp \
# 	populate.cpp \
# 	quadratic_cuts_check.cpp \
# 	report2.cpp \
# 	report.cpp \
# 	sdpcuts.cpp \
# 	sdpsol.cpp \
# 	sedumi.cpp \
# 	tracerreport.cpp

#	disjunctive_cuts.cpp 
#	tracer.cpp 
#	rlt_cuts.cpp 

libCouenneSdpCuts_la_DEPENDENCIES = $(libCouenneSdpCuts_la_LIBADD)

# This is for libtool (on Windows)
libCouenneSdpCuts_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` 

if COIN_HAS_COUENNE
  AM_CPPFLAGS += 
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem`
endif

if COIN_HAS_NTY
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)`
endif
