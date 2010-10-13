# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 405 2006-07-10 18:52:52Z andreasw $

# Authors:  Andreas Waechter           IBM    2006-04-13
#           Pietro Belotti

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCouenne                             #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCouenne.la

libCouenne_la_SOURCES = 

libCouenne_la_LIBADD = \
	./interfaces/libCouenneInterfaces.la \
	./heuristics/libCouenneHeuristics.la \
	./problem/libCouenneProblem.la \
	./standardize/libCouenneStandardize.la \
	./expression/libCouenneExpression.la \
	./branch/libCouenneBranch.la \
	./convex/libCouenneConvex.la \
	./bound_tightening/libCouenneBoundTightening.la \
	./util/libCouenneUtil.la \
	./disjunctive/libCouenneDisjunctive.la \
	./sdpcuts/libCouenneSdpCuts.la \
	./ellipcuts/libCouenneEllipCuts.la \
	./crossconv/libCouenneCrossConv.la \
	./two_implied_bt/libCouenneTwoImplied.la

#	./readnl/libCouenneReadnl.la 

libCouenne_la_DEPENDENCIES = $(libCouenne_la_LIBADD)

# This is for libtool (on Windows)
libCouenne_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
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
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`\
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

if COIN_HAS_NTY
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(NTYINCDIR)`
endif

#if COIN_HAS_FILTERSQP
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces/Filter` 
#endif

AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/interfaces` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/heuristics` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/disjunctive` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/sdpcuts` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/main` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/sdpcuts` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/ellipcuts` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/crossconv` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/two_implied_bt`

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
	util/CouenneRootQ.hpp \
	util/CouenneFunTriplets.hpp \
	main/BonCouenneInfo.hpp \
	main/BonCouenneSetup.hpp \
	interfaces/BonCouenneInterface.hpp \
	interfaces/CouenneTNLP.hpp \
	interfaces/CouenneMINLPInterface.hpp \
	heuristics/BonNlpHeuristic.hpp \
	heuristics/BonInitHeuristic.hpp \
	heuristics/CouenneFeasPump.hpp \
	branch/CouenneObject.hpp \
	branch/CouenneVarObject.hpp \
	branch/CouenneVTObject.hpp \
	branch/CouenneThreeWayBranchObj.hpp \
	branch/CouenneBranchingObject.hpp \
	branch/CouenneChooseVariable.hpp \
	branch/CouenneChooseStrong.hpp \
	branch/CouenneComplObject.hpp \
	branch/CouenneComplBranchingObject.hpp \
	branch/CouenneOrbitObj.hpp \
	branch/CouenneOrbitBranchingObj.hpp \
	branch/CouenneSOSObject.hpp \
	branch/CouenneProjections.hpp \
	standardize/CouenneLQelems.hpp \
	bound_tightening/CouenneFixPoint.hpp \
	expression/CouenneDomain.hpp \
	expression/operators/CouenneExprAbs.hpp \
	expression/operators/CouenneExprExp.hpp \
	expression/operators/CouenneExprInv.hpp \
	expression/operators/CouenneExprMul.hpp \
	expression/operators/CouenneExprOpp.hpp \
	expression/operators/CouenneExprLog.hpp \
	expression/operators/CouenneExprDiv.hpp \
	expression/operators/CouenneExprMin.hpp \
	expression/operators/CouenneExprSum.hpp \
	expression/operators/CouenneExprMax.hpp \
	expression/operators/CouenneExprPow.hpp \
	expression/operators/CouenneExprSin.hpp \
	expression/operators/CouenneExprCos.hpp \
	expression/operators/CouenneExprSub.hpp \
	expression/operators/CouenneExprGroup.hpp \
	expression/operators/CouenneExprQuad.hpp \
	expression/operators/CouenneExprIf.hpp \
	expression/operators/CouenneExprNorm.hpp \
	expression/operators/CouenneExprPWLinear.hpp \
	expression/operators/bounds/CouenneExprBMul.hpp \
	expression/operators/bounds/CouenneExprBDiv.hpp \
	expression/operators/bounds/CouenneExprBSin.hpp \
	expression/operators/bounds/CouenneExprBCos.hpp \
	expression/CouenneExprConst.hpp \
	expression/CouenneExprIVar.hpp \
	expression/CouenneExprVar.hpp \
	expression/CouenneExprAux.hpp \
	expression/CouenneExprUnary.hpp \
	expression/CouenneExpression.hpp \
	expression/CouenneExprClone.hpp \
	expression/CouenneExprBound.hpp \
	expression/CouenneExprCopy.hpp \
	expression/CouenneExprStore.hpp \
	expression/CouenneExprOp.hpp \
	expression/CouExpr.hpp \
	expression/CouennePrecisions.hpp \
	expression/partial/CouenneExprHess.hpp \
	expression/partial/CouenneExprJac.hpp \
	problem/depGraph/CouenneDepGraph.hpp \
	problem/CouenneProblemElem.hpp \
	problem/CouenneProblem.hpp \
	problem/CouenneSolverInterface.hpp \
	problem/CouenneJournalist.hpp \
	problem/CouenneGlobalCutOff.hpp \
	expression/CouenneTypes.hpp \
	disjunctive/CouenneDisjCuts.hpp \
	convex/CouenneCutGenerator.hpp \
	sdpcuts/CouenneSdpCuts.hpp \
	ellipcuts/CouenneEllipCuts.hpp \
	crossconv/CouenneCrossConv.hpp \
	two_implied_bt/CouenneTwoImplied.hpp

if COIN_HAS_NTY
includecoin_HEADERS += \
	branch/Nauty.h
endif

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CouenneConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
