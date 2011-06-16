# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
## $Id$
#
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
	./bound_tightening/twoImpliedBT/libCouenneTwoImplied.la \
	./util/libCouenneUtil.la \
	./disjunctive/libCouenneDisjunctive.la 

#	./cut/crossconv/libCouenneCrossConv.la
#	./cut/sdpcuts/libCouenneSdpCuts.la
#	./ellipcuts/libCouenneEllipCuts.la 
#	./readnl/libCouenneReadnl.la 

libCouenne_la_DEPENDENCIES = $(libCouenne_la_LIBADD)

# This is for libtool
libCouenne_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.

AM_CPPFLAGS =

# AM_CPPFLAGS = $(COUENNELIB_CFLAGS) \
# 	-I`$(CYGPATH_W) $(srcdir)/convex` \
# 	-I`$(CYGPATH_W) $(srcdir)/interfaces` \
# 	-I`$(CYGPATH_W) $(srcdir)/heuristics` \
# 	-I`$(CYGPATH_W) $(srcdir)/expression` \
# 	-I`$(CYGPATH_W) $(srcdir)/expression/operators` \
# 	-I`$(CYGPATH_W) $(srcdir)/disjunctive` \
# 	-I`$(CYGPATH_W) $(srcdir)/main` \
# 	-I`$(CYGPATH_W) $(srcdir)/problem` \
# 	-I`$(CYGPATH_W) $(srcdir)/branch` \
# 	-I`$(CYGPATH_W) $(srcdir)/cut/sdpcuts` \
# 	-I`$(CYGPATH_W) $(srcdir)/cut/crossconv` \
# 	-I`$(CYGPATH_W) $(srcdir)/two_implied_bt`

#	-I`$(CYGPATH_W) $(srcdir)/ellipcuts` 

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
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
	bound_tightening/CouenneAggrProbing.hpp \
	bound_tightening/twoImpliedBT/CouenneTwoImplied.hpp \
	expression/CouenneDomain.hpp \
	expression/operators/CouenneExprAbs.hpp \
	expression/operators/CouenneExprExp.hpp \
	expression/operators/CouenneExprInv.hpp \
	expression/operators/CouenneExprMul.hpp \
	expression/operators/CouenneExprTrilinear.hpp \
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
	problem/CouenneRecordBestSol.hpp \
	expression/CouenneTypes.hpp \
	disjunctive/CouenneDisjCuts.hpp \
	convex/CouenneCutGenerator.hpp

# cut/crossconv/CouenneCrossConv.hpp
# cut/sdpcuts/CouenneSdpCuts.hpp 
# cut/sdpcuts/dsyevx_wrapper.hpp 
# cut/sdpcuts/linquad_cuts.hpp 
# cut/sdpcuts/populate.hpp 
# cut/sdpcuts/quadratic_cuts_check.hpp 
# cut/sdpcuts/rlt_cuts.hpp

# ellipcuts/CouenneEllipCuts.hpp 

if COIN_HAS_NTY
includecoin_HEADERS += branch/Nauty.h
endif

install-exec-local:
	$(install_sh_DATA) config_couenne.h $(DESTDIR)$(includecoindir)/CouenneConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/CouenneConfig.h
