/*
 *
 * Name:    depGraph.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for manipulating dependencies between variables
 *
 * (C) Carnegie-Mellon University, 2007-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <cstdlib>
#include <cstdio>

#include "CouenneDepGraph.hpp"

using namespace Couenne;

//#define DEBUG

// Methods of the class DepNode ///////////////////////////////////////////////////////////

/// does this variable depend on variable with index xi? (used in checkCycles)
bool DepNode::depends (int xi, bool recursive,
		       std::set <DepNode *, compNode> *already_visited) const {

  // check if any node of the forward star has index xi
  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin ();
       i != depList_ -> end (); ++i) {

#ifdef DEBUG
    printf ("checking dependence of %d on %d\n", (*i) -> Index (), xi);
    //    fflush (stdout);
#endif

    if (!already_visited ||
	(already_visited -> find (*i) == already_visited -> end ())) {

      if (((*i) -> Index () == xi) || // check current node
	(recursive &&
	 ((*i) -> depends (xi, recursive, already_visited)))) // check deplist recursively
      {
#ifdef DEBUG
	printf ("%d <- ", (*i) -> Index ());
	fflush (stdout);
#endif
	return true;
      } else {
	if (already_visited) {
	  already_visited -> insert (*i);
	  /*printf ("checked [%d]: ", (*i) -> Index ());
	  for (std::set <DepNode *, compNode>::iterator j = already_visited -> begin ();
	       j != already_visited -> end (); ++j)
	    printf ("%d ", (*j) -> Index ());
	    printf ("\n");*/
	}
      }
    }
  }

  return false;
}


/// assign numbering to all nodes of graph
void DepNode::createOrder (DepGraph *g) {

  if (order_ != -1) return;

  if (order_ == -2) {

    printf ("detected cycle in creating order, exiting\n");
    exit (-1);
  }

  order_ = -2;

  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
       i != depList_ -> end (); ++i)
    if ((*i) -> Order () == -1)
      (*i) -> createOrder (g);

  if (order_ == -2)
    order_ = g -> Counter () ++;
}


/// debugging procedure
void DepNode::print (int indent, bool descend) const {

  printf ("%d ", index_);
  if (order_ >= 0) printf ("[%d]", order_);
  fflush (stdout);

  if (depList_ -> size () > 0) {
    printf ("("); fflush (stdout);

    for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
	 i != depList_ -> end (); ++i)
      if (descend)
	(*i) -> print (indent + 1, descend);
      else printf ("%d ", (*i) -> Index ());

    printf (") "); fflush (stdout);
  }
}

/// replace the index of a variable with another in the entire
/// graph. Used when redundant constraints w := x are discovered
void DepNode::replaceIndex (DepNode *oldVarNode, DepNode *newVarNode) {

  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
       i != depList_ -> end (); ++i)

    if ((*i) -> Index () == oldVarNode -> Index ()) {

      depList_ -> erase (i);

      if (depList_ -> find (newVarNode) == depList_ -> end ())
	depList_ -> insert (newVarNode);

      break;
    }
}

// Methods of the class DepGraph ////////////////////////////////////////////////////////////


/// insert new variable if new
void DepGraph::insert (exprVar *var) {

  DepNode *el = new DepNode (var -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  if (i == vertices_ . end ())
    vertices_.insert (el);
  else delete el;
}


/// insert new auxiliary if new
void DepGraph::insert (exprAux *aux) {

  if (!aux) return;

  DepNode *el = new DepNode (aux -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  if (i == vertices_ . end ()) {
    vertices_.insert (el);
    aux -> Image () -> fillDepSet (el -> DepList (), this);
  } else {
    aux -> Image () -> fillDepSet ((*i) -> DepList (), this);
    delete el;
  }
}


/// erase element from graph
void DepGraph::erase (exprVar *var) {

  DepNode *el = new DepNode (var -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  if (i != vertices_ . end ())
    vertices_.erase (i);
  delete el;
}

/// does w depend on x?
bool DepGraph::depends (int wi, int xi, bool recursive) {

  DepNode *el = new DepNode (wi);
  std::set <DepNode *, compNode>::iterator i = vertices_. find (el);
  delete el;

  if (i != vertices_. end ()) {              // if such element is in the set
    std::set <DepNode *, compNode> already_visited;
    return (*i) -> depends (xi, recursive, &already_visited); // then search it
  }
  else return false;
}


/// assign numbering to all nodes of graph
void DepGraph::createOrder () {

  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); ++i)
    (*i) -> createOrder (this);
}


/// debugging procedure
void DepGraph::print (bool descend) {

  printf ("Dependence graph: \n");
  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); ++i) {
    (*i) -> print (0, descend);
    printf ("\n");
  }
}


/// search for node in vertex set
DepNode *DepGraph::lookup (int index) {

  DepNode *el = new DepNode (index), *ret;
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  ret = (i == vertices_.end ()) ? NULL : (*i);

  delete el;
  return ret;
}

/// replace, throughout the whole graph, the index of a variable with
/// another in the entire graph. Used when redundant constraints w :=
/// x are discovered
void DepGraph::replaceIndex (int oldVar, int newVar) {

  DepNode
    *copyOld = new DepNode (oldVar),
    *copyNew = new DepNode (newVar);

  std::set <DepNode *, compNode>::iterator
    oldNode = vertices_.find (copyOld),
    newNode = vertices_.find (copyNew);

  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); ++i)
    (*i) -> replaceIndex (*oldNode, *newNode);

  delete copyOld;
  delete copyNew;
}
