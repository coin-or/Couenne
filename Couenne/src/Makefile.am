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
	./standardize/libCouenneStandardize.la \
	./expression/libCouenneExpression.la \
	./branch/libCouenneBranch.la \
	./convex/libCouenneConvex.la \
	./bound_tightening/libCouenneBoundTightening.la \
	./util/libCouenneUtil.la \
	./disjunctive/libCouenneDisjunctive.la

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
	-I`$(CYGPATH_W) $(CBCSRCDIR)/examples` \
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
	expression/CouExpr.hpp \
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

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
