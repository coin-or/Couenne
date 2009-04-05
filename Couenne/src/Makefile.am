# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 405 2006-07-10 18:52:52Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCouenne                             #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCouenne.la

libCouenne_la_SOURCES = 

libCouenne_la_LIBADD = \
	./problem/libCouenneProblem.la \
	./readnl/libCouenneReadnl.la \
	./expression/libCouenneExpression.la \
	./branch/libCouenneBranch.la \
	./convex/libCouenneConvex.la \
	./util/libCouenneUtil.la \
	./standardize/libCouenneStandardize.la \
	./bound_tightening/libCouenneBoundTightening.la \
	./disjunctive/libCouenneDisjunctive.la \
	../../Bonmin/src/Interfaces/Ampl/BonAmplInterface.lo \
	../../Bonmin/src/Interfaces/Ampl/BonAmplTMINLP.lo \
	../../Bonmin/src/Interfaces/Ampl/sos_kludge.lo \
	../../Bonmin/src/Algorithms/Ampl/BonAmplSetup.lo \
        $(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la

libCouenne_la_DEPENDENCIES = $(libCouenne_la_LIBADD)

# This is for libtool (on Windows)
libCouenne_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS =

# If we have the Ampl solver library, we need to add additional things
if COIN_HAS_ASL
AM_CPPFLAGS += $(ASL_CPPFLAGS)
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_couenne.h \
	util/rootQ.hpp \
	util/funtriplets.hpp \
	main/BonCouenneInfo.hpp \
	main/BonCouenneSetup.hpp \
	main/BonCouenneInterface.hpp \
	main/BonNlpHeuristic.hpp \
	main/BonInitHeuristic.hpp \
	branch/CouenneObject.hpp \
	branch/CouenneVarObject.hpp \
	branch/CouenneVTObject.hpp \
	branch/CouenneThreeWayBranchObj.hpp \
	branch/CouenneBranchingObject.hpp \
	branch/CouenneChooseVariable.hpp \
	branch/CouenneChooseStrong.hpp \
	branch/CouenneComplObject.hpp \
	branch/CouenneComplBranchingObject.hpp \
	branch/CouenneSOSObject.hpp \
	branch/projections.hpp \
	standardize/lqelems.hpp \
	expression/domain.hpp \
	expression/operators/exprAbs.hpp \
	expression/operators/exprExp.hpp \
	expression/operators/exprInv.hpp \
	expression/operators/exprMul.hpp \
	expression/operators/exprOpp.hpp \
	expression/operators/exprLog.hpp \
	expression/operators/exprDiv.hpp \
	expression/operators/exprMin.hpp \
	expression/operators/exprSum.hpp \
	expression/operators/exprMax.hpp \
	expression/operators/exprPow.hpp \
	expression/operators/exprSin.hpp \
	expression/operators/exprCos.hpp \
	expression/operators/exprSub.hpp \
	expression/operators/exprGroup.hpp \
	expression/operators/exprQuad.hpp \
	expression/operators/exprIf.hpp \
	expression/operators/exprNorm.hpp \
	expression/operators/exprPWLinear.hpp \
	expression/operators/bounds/exprBMul.hpp \
	expression/operators/bounds/exprBDiv.hpp \
	expression/operators/bounds/exprBSin.hpp \
	expression/operators/bounds/exprBCos.hpp \
	expression/exprConst.hpp \
	expression/exprIVar.hpp \
	expression/exprVar.hpp \
	expression/exprAux.hpp \
	expression/exprUnary.hpp \
	expression/expression.hpp \
	expression/exprClone.hpp \
	expression/exprBound.hpp \
	expression/exprCopy.hpp \
	expression/exprStore.hpp \
	expression/exprOp.hpp \
	expression/CouennePrecisions.hpp \
	problem/depGraph/depGraph.hpp \
	problem/CouenneProblemElem.hpp \
	problem/CouenneProblem.hpp \
	problem/CouenneSolverInterface.hpp \
	problem/CouenneJournalist.hpp \
	expression/CouenneTypes.hpp \
	disjunctive/CouenneDisjCuts.hpp \
	convex/CouenneCutGenerator.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CouenneConfig.h

if COIN_HAS_ASL
install-asl:
	$(install_sh_DATA) ../../ThirdParty/ASL/amplsolv* $(DESTDIR)$(libdir)
else
install-asl:

endif

install-exec-local: install-asl
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
