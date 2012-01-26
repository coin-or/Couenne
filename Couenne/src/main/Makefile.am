# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author: Pietro Belotti, on a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

# what is built at this level

lib_LTLIBRARIES = libBonCouenne.la

libBonCouenne_la_SOURCES = \
	CouenneBab.cpp \
	BonCouenneSetup.cpp \
	SmartAslDestr.cpp

if COIN_HAS_ASL
  bin_PROGRAMS = couenne 
endif

couenne_SOURCES = BonCouenne.cpp
couenne_LDADD   = libBonCouenne.la ../libCouenne.la ../readnl/libCouenneReadnl.la $(COUENNELIB_LIBS) $(ASL_LIBS) 
couenne_DEPENDENCIES  = libBonCouenne.la ../libCouenne.la ../readnl/libCouenneReadnl.la $(COUENNELIB_DEPENDENCIES) $(ASL_DEPENDENCIES)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../convex` \
	-I`$(CYGPATH_W) $(srcdir)/../expression` \
	-I`$(CYGPATH_W) $(srcdir)/../expression/operators` \
	-I`$(CYGPATH_W) $(srcdir)/../interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../heuristics` \
	-I`$(CYGPATH_W) $(srcdir)/../branch` \
	-I`$(CYGPATH_W) $(srcdir)/../problem` \
	-I`$(CYGPATH_W) $(srcdir)/../cut/crossconv` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening` \
	-I`$(CYGPATH_W) $(srcdir)/../bound_tightening/twoImpliedBT` \
	-I`$(CYGPATH_W) $(srcdir)/../disjunctive` \
	-I`$(CYGPATH_W) $(srcdir)/../readnl` \
	$(COUENNELIB_CFLAGS) $(ASL_CFLAGS) 

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..
