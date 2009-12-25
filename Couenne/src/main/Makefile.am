# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 105 2006-11-03 19:31:19Z pbonami $

# Author: Pietro Belotti, on a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

# what is built at this level

if COIN_HAS_ASL
  bin_PROGRAMS = couenne
endif

lib_LTLIBRARIES = libBonCouenne.la

libBonCouenne_la_SOURCES = \
	BonCouenneSetup.cpp \
	BonCouenneInterface.cpp \
	BonNlpHeuristic.cpp \
	BonInitHeuristic.cpp \
	CouenneOSInterface.cpp

COINLIBS = \
	$(CBCOBJDIR)/src/libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
        $(BONMINOBJDIR)/src/CbcBonmin/libbonmin.la \
        $(BONMINOBJDIR)/src/Algorithms/Branching/libbonbranching.la \
	$(BONMINOBJDIR)/src/Interfaces/libbonmininterfaces.la \
        $(IPOPTOBJDIR)/src/Interfaces/libipopt.la \
        $(IPOPTOBJDIR)/src/LinAlg/liblinalg.la

if COIN_HAS_CPX
     COINLIBS += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif

COINLIBS += \
     $(OSIOBJDIR)/src/libOsi.la \
     $(CLPOBJDIR)/src/libClp.la \
     $(COINUTILSOBJDIR)/src/libCoinUtils.la

IPOPTLIBS = \
        $(IPOPTOBJDIR)/src/Interfaces/libipopt.la \
        $(IPOPTOBJDIR)/src/LinAlg/liblinalg.la

couenne_LDADD         = libBonCouenne.la
couenne_DEPENDENCIES  = libBonCouenne.la
if COIN_HAS_ASL
couenne_LDADD        += ../readnl/libCouenneReadnl.la
couenne_DEPENDENCIES += ../readnl/libCouenneReadnl.la
endif
couenne_LDADD        += ../libCouenne.la
couenne_DEPENDENCIES += ../libCouenne.la

couenne_SOURCES = \
	BonCouenne.cpp

couenne_LDADD += \
	$(LIBCBCCOPY) \
	$(COINLIBS) \
	$(ADDLIBS) \
	$(ASLLIB) \
	`cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`


# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
couenne_LDFLAGS = -rpath $(libdir)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/examples` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`\
	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLandP` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglRedSplit` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess`

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

AM_CPPFLAGS += $(ASL_CPPFLAGS)

#if COIN_HAS_FILTERSQP
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces/Filter` 
#endif

AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/disjunctive` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/main` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc
