# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src/util \
	src/expression \
	src/standardize \
	src/problem \
	src/bound_tightening \
	src/bound_tightening/twoImpliedBT \
	src/convex \
	src/branch \
	src/disjunctive \
	src/interfaces \
	src/heuristics \
	src/cut/sdpcuts \
	src/cut/crossconv \
	src \
	src/main

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test && $(MAKE) test

unitTest: test

# Doxygen documentation

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )

clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                  Installation of the pkg-config file                 #
########################################################################

pkgconfiglib_DATA = couenne.pc
