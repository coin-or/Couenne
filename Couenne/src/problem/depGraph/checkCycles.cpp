/* $Id$
 *
 * Name:    checkCycles.cpp
 * Author:  Pietro Belotti
 * Purpose: check for cycles in dependence graph
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

//#include <stdio.h>
#include "CouenneDepGraph.hpp"

//#define DEBUG

using namespace Couenne;

static bool visit (std::set <DepNode *, compNode>::iterator &v);

/// check for cycles in dependence graph

bool DepGraph::checkCycles () {

  for (std::set <DepNode *, compNode>::iterator 
	 i = vertices_.begin ();
       i  != vertices_.end   (); ++i)
    (*i) -> color () = DepNode::DEP_WHITE; 

  // simple DFS that checks for cycles

  for (std::set <DepNode *, compNode>::iterator 
	 i = vertices_.begin ();
       i  != vertices_.end   (); ++i)

    if (((*i) -> color () == DepNode::DEP_WHITE) &&
	(visit (i)))
      return true;

  return false;
}

// subroutine that visits a single node
bool visit (std::set <DepNode *, compNode>::iterator &v) {

  //printf ("%d is gray\n", (*v) -> Index ());
  (*v) -> color () = DepNode::DEP_GRAY; 
  std::set <DepNode *, compNode> *list = (*v) -> DepList ();

  // gen2 contains the adjacency list of this node

  for (std::set <DepNode *, compNode>::iterator 
	 j = list -> begin ();
       j  != list -> end   (); ++j)

    if ((*j) -> color () == DepNode::DEP_GRAY) {
      //printf ("%d is gray\n", (*j) -> Index ());
      return true;
    }
    else {

      if ((*j) -> color () == DepNode::DEP_WHITE)
	//printf ("visiting %d\n", (*j) -> Index ());

      if (((*j) -> color () == DepNode::DEP_WHITE) && (visit (j))) {
	//printf ("%d is true\n", (*j) -> Index ());
	return true;
      }
    }

  (*v) -> color () = DepNode::DEP_BLACK;
  //printf ("%d is black\n", (*v) -> Index ());
  return false;
}
