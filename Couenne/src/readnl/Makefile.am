# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# Author: Pietro Belotti, from a template by Pierre Bonami

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libCouenneReadnl.la

# List all source files for this library, including headers
if COIN_HAS_ASL
libCouenneReadnl_la_SOURCES =
#libCouenneReadnl_la_SOURCES = \
#	readnl.cpp \
#	nl2e.cpp \
#	invmap.cpp \
#	readASLfg.cpp

libCouenneReadnl_la_LIBADD = \
	../../../Bonmin/src/Interfaces/Ampl/BonAmplInterface.lo \
	../../../Bonmin/src/Interfaces/Ampl/BonAmplTMINLP.lo \
	../../../Bonmin/src/Interfaces/Ampl/sos_kludge.lo \
	../../../Bonmin/src/Algorithms/Ampl/BonAmplSetup.lo \
        $(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la
endif

libCouenneReadnl_la_DEPENDENCIES = $(libCouenneReadnl_la_LIBADD)

# This is for libtool (on Windows)
libCouenneReadnl_la_LDFLAGS = $(LT_LDFLAGS)


# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` 

if COIN_HAS_ASL
AM_CPPFLAGS += $(ASL_CPPFLAGS)
endif

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` 
endif

if COIN_HAS_ASL
install-exec-local:
	$(install_sh_DATA) ../../../ThirdParty/ASL/amplsolv* $(DESTDIR)$(libdir)
endif
