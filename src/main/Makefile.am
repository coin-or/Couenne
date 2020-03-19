# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author: Pietro Belotti, on a template by Pierre Bonami

if COIN_HAS_ASL
  bin_PROGRAMS = couenne 
endif

couenne_SOURCES = BonCouenne.cpp
couenne_LDADD   = ../libCouenne.la

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../convex \
	-I$(srcdir)/../expression \
	-I$(srcdir)/../expression/operators \
	-I$(srcdir)/../expression/partial \
	-I$(srcdir)/../interfaces \
	-I$(srcdir)/../heuristics \
	-I$(srcdir)/../branch \
	-I$(srcdir)/../problem \
	-I$(srcdir)/../bound_tightening \
	-I$(srcdir)/../bound_tightening/twoImpliedBT \
	-I$(srcdir)/../disjunctive \
	-I$(srcdir)/../cut/sdpcuts \
	-I$(srcdir)/../cut/crossconv \
	-I$(srcdir)/../readnl \
	$(COUENNELIB_CFLAGS) -DCOUENNE_BUILD

AM_LDFLAGS = $(LT_LDFLAGS)
