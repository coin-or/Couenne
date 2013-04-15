/* $Id$
 *
 * Name:    generateDisjCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: separation method for disjunctive cuts
 *
 * (C) Carnegie-Mellon University, 2008-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>
#include "CoinTime.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneDisjCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneInfeasCut.hpp"

using namespace Ipopt;
using namespace Couenne;

/// generate disjunctive cuts
void CouenneDisjCuts::generateCuts (const OsiSolverInterface &si, 
				    OsiCuts &cs, 
				    const CglTreeInfo info)
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
  const
#endif
{

  // Check if cs contains only one cut and if it is of the form 1 <=
  // x0 <= -1. That means a previous cut generator has determined that
  // this node is infeasible and we shouldn't take the pain of running
  // this CGL.

  if (isWiped (cs))
    return;
  
  if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS)) 
    printf ("--- generateDisjCuts: level = %d, pass = %d, intree = %d [%d]\n",
	    info.level, info.pass, info.inTree, depthStopSeparate_);

  if ((depthStopSeparate_ >= 0) &&        // if -1 no limit on depth
      (info.level > depthStopSeparate_))  // check if too deep for adding these cuts
    return;

  int nInitCuts = nrootcuts_;

  if ((info.level <= 0) && !(info.inTree)) {

    jnlst_ -> Printf (J_ERROR, J_COUENNE, 
		      "Disjunctive cuts (root node): ");
    fflush (stdout);
  }

  double time = CoinCpuTime ();

  // use clone of solver interface
  OsiSolverInterface *csi = si.clone ();

  int
    initRowCuts = cs.sizeRowCuts (),
    initColCuts = cs.sizeColCuts ();

  // get disjunctions /////////////////////////////
  //
  // start:
  //
  // consider problem Ax <= a_0, x in [l,u]
  //
  // // A and a_0 may be, depending on depth of BB tree and size of the
  // // problem, from the linearization at the root node (in this case
  // // a clone() is needed at the first iteration and branching rules
  // // should be applied explicitly) or the current LP as passed by
  // // si.
  //
  // Get set of disjunctions (in the form of OsiBranchingObjects) in
  // number limited by depth, problem size, and sorted according to
  // branching method (HotInfo if strong branching?)
  //
  // preprocessing /////////////////////////////////
  //
  // for each disjunction (x_i <= or >= x_i^d) {
  //
  //   1) apply left  disj., bound tightening returns x in [l_1,u_1]
  //   2) apply right disj., bound tightening returns x in [l_2,u_2]
  //
  //   3) if both infeasible, bail out
  //   4) if one infeasible only, save column cut, apply tightening to
  //      si and goto start
  // }
  //
  // CGLP //////////////////////////////////////////
  //
  // for each disjunction (x_i <= or >= x_i^d) {
  //
  //   5) get cuts Bx <= b_0 for left  disj.
  //   6) get cuts Cx <= c_0 for right disj.
  //
  //   7) if both infeasible, bail out
  //   8) if one infeasible only, save column cut, apply tightening to
  //      si and continue
  //
  //   9) generate CGLP:
  //         consider subset of Ax <= a_0: 
  //           a) active only, or
  //           b) {root linearization cuts} union {currently active}, or
  //           c) all cuts (!)
  //
  //   10) solve CGLP
  //
  //   11) add corresponding cut, possibly to CGLP as well?
  // }

  int maxDisj = (initDisjNumber_ >= 0) ? 
    CoinMin ((int) (csi -> numberObjects () * initDisjPercentage_), initDisjNumber_) :
    (int) (csi -> numberObjects () * initDisjPercentage_);

  // number of disjunctions to consider (branching objects)
  numDisjunctions_ = (depthLevelling_ < 0 || info.level < depthLevelling_) ? 
    (int) (maxDisj) :
    (int) (maxDisj / (2 + info.level - depthLevelling_));

  if (numDisjunctions_ < 1) numDisjunctions_ = 1;

  const int 
    exc_infeasible = 1, 
    exc_normal     = 2, 
    max_iterations = 1;

  couenneCG_ -> Problem () -> domain () -> push (&si, &cs, false);

  std::vector <std::pair <OsiCuts *, OsiCuts *> > disjunctions;

  bool infeasNode = false;

  try {

    // get disjunctions (rows and cols) ////////////////////////////////////////////////////

    bool start_over;
    int iterations = 0;

    do { // repeat as long as preprocessing or separation of convCuts
	 // shrink the whole problem

      ++iterations;

      start_over = false;

      // preprocess, get column cuts (disjoint bounding boxes)
      int result = getDisjunctions (disjunctions, *csi, cs, info);

      if      (result == COUENNE_INFEASIBLE) throw exc_infeasible; // fathom node
      else if (result == COUENNE_TIGHTENED && 
	       iterations < max_iterations) start_over = true;     // tightened node

      if (disjunctions.empty ())
	throw exc_normal;

      // generate convexification cuts for each disjunction
      for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	   disjI != disjunctions.end (); ++disjI) {

	if (CoinCpuTime () > couenneCG_ -> Problem () -> getMaxCpuTime ()) {
	  start_over = false;
	  break;
	}

	// separate on single disjunction

	// left
	result = separateWithDisjunction (disjI -> first, *csi, cs, info);
	if      (result == COUENNE_INFEASIBLE) throw exc_infeasible;           // fathom node
	else if (result == COUENNE_TIGHTENED && iterations < max_iterations) { // tightened node
	  start_over = true;
	  break;
	}

	// right
	result = separateWithDisjunction (disjI -> second, *csi, cs, info);
	if      (result == COUENNE_INFEASIBLE) throw exc_infeasible;           // fathom node
	else if (result == COUENNE_TIGHTENED && iterations < max_iterations) { // tightened node
	  start_over = true;
	  break;
	}
      }

      if (start_over) {

	for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	     disjI != disjunctions.end (); ++disjI) {
	  delete disjI -> first;
	  delete disjI -> second;
	}

	disjunctions.erase (disjunctions.begin (), disjunctions.end ());
      }

      if (!start_over && jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))

	// generate convexification cuts for each disjunction
	for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	     disjI != disjunctions.end (); ++disjI) {

	  printf ("=========================== CUTS for the LEFT part\n");
	  for (int i=0; i<disjI->first->sizeColCuts (); i++) disjI->first->colCutPtr(i)->print();
	  for (int i=0; i<disjI->first->sizeRowCuts (); i++) disjI->first->rowCutPtr(i)->print();
	  printf ("=========================== CUTS for the RIGHT part\n");
	  for (int i=0; i<disjI->second->sizeColCuts (); i++) disjI->second->colCutPtr(i)->print();
	  for (int i=0; i<disjI->second->sizeRowCuts (); i++) disjI->second->rowCutPtr(i)->print();
	  printf ("===========================\n");
        }

    } while (start_over && (iterations < max_iterations));

    // si contains the tightest bounding box. Use it to update
    // CouenneProblem's bounds AND add to cs

    // already done above

    // CGLP //////////////////////////////////////////////////////////////////////////

    // maybe one last FBBT before big CGLP?

    // generate all cuts
    if (generateDisjCuts (disjunctions, *csi, cs, info) == COUENNE_INFEASIBLE) // node can be fathomed
      throw exc_infeasible;
  }

  catch (int exception) {

    if (exception == exc_infeasible) { // add infeasible column cut 1 <= x_0 <= -1

      jnlst_ -> Printf (J_DETAILED, J_DISJCUTS, "--- Disjunctive Cut separator: infeasible node\n");
      WipeMakeInfeas (cs);
      infeasNode = true;
    }
  }

  // cleanup
  for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
       disjI != disjunctions.end (); ++disjI) {

    delete disjI -> first;
    delete disjI -> second;
  }

  couenneCG_ -> Problem () -> domain () -> pop ();

  if (!infeasNode) {

    // tighten bounds of si based on those tightened of csi

    CoinPackedVector 
      tighterLower, 
      tighterUpper;

    const double 
      *oldLo = si. getColLower (),   *newLo = csi -> getColLower (), 
      *oldUp = si. getColUpper (),   *newUp = csi -> getColUpper ();

    int ncols = si.getNumCols ();

    bool tightened = false;

    for (int i=0; i<ncols; i++, newLo++, newUp++) {

      if (*newLo > *oldLo++ + COUENNE_EPS) {tighterLower.insert (i, *newLo); tightened = true;}
      if (*newUp < *oldUp++ - COUENNE_EPS) {tighterUpper.insert (i, *newUp); tightened = true;}
    }

    if (tightened) {
      OsiColCut tighter;
      tighter.setLbs (tighterLower);
      tighter.setUbs (tighterUpper);
      if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS)) {
	printf ("tightened bounds in disjunctive cuts:");
	tighter.print ();
      }
      cs.insert (tighter);
    }

    int deltaNcuts = 
      cs.sizeRowCuts () - initRowCuts + 
      cs.sizeColCuts () - initColCuts;

    if (info.level <= 0 && !(info.inTree)) 
      nrootcuts_ += deltaNcuts;
    ntotalcuts_ += deltaNcuts;

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS)) {

      if (cs.sizeRowCuts()>initRowCuts) printf ("added %d row cuts\n", cs.sizeRowCuts () - initRowCuts);
      if (cs.sizeColCuts()>initColCuts) printf ("added %d col cuts\n", cs.sizeColCuts () - initColCuts);
    }

    if ((info.level <= 0) && !(info.inTree))
      jnlst_ -> Printf (J_ERROR, J_COUENNE, 
			"%d cuts\n", CoinMax (0, nrootcuts_ - nInitCuts));
    else if (deltaNcuts)
      jnlst_ -> Printf (J_WARNING, J_COUENNE, 
			"In-BB disjunctive cuts: %d row cuts, %d col cuts\n",
			cs.sizeRowCuts () - initRowCuts,
			cs.sizeColCuts () - initColCuts);
  }

  // else {
  //     jnlst_ -> Printf (J_STRONGWARNING, J_COUENNE, 
  // 			"In-BB disjunctive cuts: %d row cuts, %d col cuts\n",
  // 			cs.sizeRowCuts () - initRowCuts,
  // 			cs.sizeColCuts () - initColCuts);
  //   if ((info.level <= 0) && !(info.inTree))
  //     jnlst_ -> Printf (J_ERROR, J_COUENNE, 
  // 			"infeasible node\n");
  // }

  delete csi;

  septime_ += (CoinCpuTime () - time);
}
