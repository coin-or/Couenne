# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Pietro Belotti

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libCouenneSdpCuts.la

# List all source files for this library
libCouenneSdpCuts_la_SOURCES = \
	CouenneSdpCuts.cpp \
	CouennePSDcon.cpp \
	CouenneMatrix.cpp \
	CutGen.cpp \
	CutGenSparse.cpp \
	dsyevx_wrapper.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../.. \
	-I$(srcdir)/../../expression \
	-I$(srcdir)/../../problem \
	$(COUENNELIB_CFLAGS) -DCOUENNE_BUILD
