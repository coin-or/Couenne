/* $Id$ */
/*
 * Name:    checkCycles.cpp
 * Author:  Pietro Belotti
 * Purpose: check for cycles in dependence graph
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDepGraph.hpp"

//#define DEBUG

/// check for cycles in dependence graph

using namespace Couenne;

bool DepGraph::checkCycles () {

  for (std::set <DepNode *, compNode>::iterator 
	 i = vertices_.begin ();
       i  != vertices_.end   (); ++i) {

    int xi = (*i) -> Index ();

    std::set <DepNode *, compNode> *gen2 = (*i) -> DepList ();

    for (std::set <DepNode *, compNode>::iterator j = gen2 -> begin (); 
	 j != gen2 -> end (); ++j) {

      std::set <DepNode *, compNode> already_visited;

      if ((*j) -> depends (xi, true, &already_visited)) {

#ifdef DEBUG
	printf ("(%d -> %d) ", xi, (*j) -> Index ());
	fflush (stdout);
#endif

	return true;
      }
      else already_visited.insert (*j);
    }
  }

  return false;
}
