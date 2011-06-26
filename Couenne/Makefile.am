# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src/util \
	  src/expression

if COIN_HAS_ASL
SUBDIRS += src/readnl
endif

SUBDIRS += src/standardize \
	  src/problem \
	  src/bound_tightening \
	  src/bound_tightening/twoImpliedBT \
	  src/convex \
	  src/branch \
	  src/disjunctive \
	  src/interfaces \
	  src/heuristics \
	  src \
	  src/main 

#	  src/cut/crossconv
#	  src/cut/sdpcuts
#	  src/cut/ellipcuts

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = test/grandTest.sh.in

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = couenne.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Couenne

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs couenne > $(addlibsdir)/couenne_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "/libpath:`$(CYGPATH_W) @abs_lib_dir@` libBonCouenne.lib libCouenne.lib @COUENNELIB_LIBS_INSTALLED@" > $(addlibsdir)/couenne_addlibs_cpp.txt
else
	echo -L@abs_lib_dir@ -lBonCouenne -lCouenne @COUENNELIB_LIBS_INSTALLED@ > $(addlibsdir)/couenne_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/couenne_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
