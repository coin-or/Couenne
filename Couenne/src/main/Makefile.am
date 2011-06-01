# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author: Pietro Belotti, on a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

# what is built at this level

lib_LTLIBRARIES = libBonCouenne.la

libBonCouenne_la_SOURCES = \
	BonCouenneSetup.cpp \
	SmartAslDestr.cpp

if COIN_HAS_ASL
  bin_PROGRAMS = couenne 
endif

couenne_SOURCES = BonCouenne.cpp
couenne_LDADD   = libBonCouenne.la ../libCouenne.la ../readnl/libCouenneReadnl.la $(COUENNELIB_LIBS) $(ASL_LIBS) 

if COIN_HAS_SCP
couenne_LDADD += \
	$(SCPLIBDIR)/libscip-2.0.1.linux.x86_64.gnu.opt-shared.so \
	$(SCPLIBDIR)/liblpispx-2.0.1.linux.x86_64.gnu.opt-shared.so \
	$(SCPLIBDIR)/libnlpi-2.0.1.linux.x86_64.gnu.opt-shared.so \
	$(SCPLIBDIR)/libobjscip-2.0.1.linux.x86_64.gnu.opt-shared.so
endif

if COIN_HAS_SPX
couenne_LDADD += \
	$(SPXLIBDIR)/libsoplex-1.5.0.linux.x86_64.gnu.opt-shared.so
endif

couenne_DEPENDENCIES  = libBonCouenne.la ../libCouenne.la ../readnl/libCouenneReadnl.la $(COUENNELIB_DEPENDENCIES) $(ASL_DEPENDENCIES)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
couenne_LDFLAGS = -rpath $(libdir)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = $(COUENNELIB_CFLAGS) $(ASL_CFLAGS) \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../heuristics` \
	-I`$(CYGPATH_W) $(srcdir)/../branch`\
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../crossconv` \
	-I`$(CYGPATH_W) $(srcdir)/../two_implied_bt` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../disjunctive` \
	-I`$(CYGPATH_W) $(srcdir)/../readnl`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc
