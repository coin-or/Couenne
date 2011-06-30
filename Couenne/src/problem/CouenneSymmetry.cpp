/* $Id$
 *
 * Name:    CouenneSymmetry.cpp
 * Author:  Jim Ostrowski
 * Purpose: methods for exploiting symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <cassert>
#include <vector>
#include <algorithm>
#include <ostream>
#include <iterator>
#include <stdio.h>

#include "CouenneExprVar.hpp"
#include "CouenneExprGroup.hpp"

#include "CouenneProblem.hpp"

using namespace Couenne;

#ifdef COIN_HAS_NTY

#include "Nauty.h"

void Node::node(int i, double c , double l, double u, int cod, int s){
  index = i;
  coeff = c;
  lb = l;
  ub = u;
  color = -1;
  code = cod;
  sign = s;
}

inline bool CouenneProblem::compare (register Node &a, register Node &b) const {
  if(a.get_code() == b.get_code() )
    if(a.get_coeff() == b.get_coeff() )
      if(a.get_sign() == b.get_sign() )
	if( fabs ( a.get_lb() - b.get_lb() ) <= COUENNE_EPS )
	  if( fabs ( a.get_ub() - b.get_ub() ) <= COUENNE_EPS )
	    return 1; 
  return 0;
}

/*
bool CouenneProblem::node_sort (Node a, Node  b){
  bool is_less = 0;

  if(a.get_code() < b.get_code() )
    is_less = 1;
  else {
    if(a.get_code() == b.get_code() )
      if(a.get_coeff() < b.get_coeff() )
	is_less = 1;
      else{
	if(a.get_coeff() ==  b.get_coeff() )
	  if(a.get_lb() < b.get_lb())
	    is_less = 1;
	  else{
	    if(a.get_lb() == b.get_lb())
	      if(a.get_ub() < b.get_ub())
		is_less = 1;
	      else{
		if(a.get_index() < b.get_index())
		  is_less = 1;
	      }
	  }
      }
  }
  return is_less;
}
bool CouenneProblem::index_sort (Node a, Node b){
  return (a.get_index() < b.get_index() );
}
*/

void CouenneProblem::sym_setup (){

  //  // Find Coefficients

  /// initialize nauty

  int num_affine = 0;

  for (std::vector <exprVar *>:: iterator i = Variables (). begin ();
       i != Variables (). end (); ++i) {

    if ((*i) -> Type () == AUX) {
      if ((*i) -> Image () -> code () == COU_EXPRDIV) {
	      num_affine ++;
      }
      if ((*i) -> Image () -> code () != COU_EXPRGROUP) {
	if ((*i) -> Image () -> Type () == N_ARY) {
	  for (int a=0; a < (*i) -> Image () -> nArgs(); a++) {
	    expression *arg = (*i) -> Image () -> ArgList () [a];

	    if (arg -> Type () == CONST) {
	      num_affine ++;

	    }
	  }
	}
      }
      if ((*i) -> Image () -> code () == COU_EXPRGROUP) {

	exprGroup *e = dynamic_cast <exprGroup *> ((*i) -> Image ());

	// add a node for e -> getC0 ();
	if (e -> getc0 () != 0 ){
	  num_affine ++;
	}

	// for each term add nodes for their non-one coefficients and their variable

	for (exprGroup::lincoeff::iterator el = e ->lcoeff().begin (); el != e -> lcoeff ().end (); ++el) {
	  if ( el -> second !=1){
	    num_affine ++;
	  }
	}
      }
    }
  }

  // Create global Nauty object

  int nc = num_affine + nVars ();
  // printf (" There are   %d  coefficient vertices in the graph \n", num_affine);
  //printf (" Graph has    %d  vertices \n", nc);

  nauty_info = new Nauty(nc);
  // create graph

  int coef_count= nVars ();
  for (std::vector <exprVar *>:: iterator i =  Variables (). begin ();
       i != Variables (). end (); ++i) {

    //    printf ("I have code %d \n",  (*i) ->  Image() -> code() );

    if ((*i) -> Type () == AUX) {
      //printf ("aux is %d with code %d \n", (*i) -> Index (), (*i) -> Image () -> code() );
      // this is an auxiliary variable

      Node vertex;
      vertex.node( (*i) -> Index () , 0.0 , (*i) -> lb () , (*i) -> ub () ,  (*i) -> Image () -> code(), (*i)-> sign() );
      //printf(" sign of aux %d \n", (*i) -> sign () );
      node_info.push_back( vertex);

      // add node in nauty graph for its index, (*i) -> Index ()

      if ((*i) -> Image () -> Type () == N_ARY) {

	if ((*i) -> Image () -> code () == COU_EXPRDIV) {
	  expression *arg = (*i) -> Image () -> ArgList () [0];
	  nauty_info->addElement((*i) -> Index (),  arg -> Index ());
	  expression *arg2 = (*i) -> Image () -> ArgList () [1];
	  nauty_info->addElement((*i) -> Index (),  coef_count);
	  nauty_info->addElement( coef_count,  arg2 -> Index ());
	  Node coef_vertex;
	  coef_vertex.node( coef_count, -1, -1 ,-1, -2 , 0);
	  node_info.push_back(coef_vertex);
	  coef_count ++;
	}
	
	if ((*i) -> Image () -> code () != COU_EXPRGROUP) {

	  for (int a=0; a < (*i) -> Image () -> nArgs(); a++) {
	    expression *arg = (*i) -> Image () -> ArgList () [a];

	    if (arg -> Type () != CONST) {
	      //printf (" add edge  %d , %d\n", (*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement((*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement( arg -> Index (), (*i) -> Index ());
	    }

	    else{

	      assert (arg -> Type () == CONST);

	      // printf (" add new vertex to graph, coef # %d, value %g \n", coef_count, arg -> Value() );
	      // printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (),  coef_count);
	      nauty_info->addElement((*i) -> Index (),  coef_count);
	      nauty_info->addElement( coef_count, (*i) -> Index ());


	      Node coef_vertex;
	      coef_vertex.node( coef_count, arg -> Value(), arg -> Value() , arg -> Value(), -2 , 0);
	      node_info.push_back(coef_vertex);
	      coef_count ++;
	    }

	  }
	}


	if ((*i) -> Image () -> code () == COU_EXPRGROUP) {

	  // dynamic_cast it to an exprGroup
	  exprGroup *e = dynamic_cast <exprGroup *> ((*i) -> Image ());

	  // add a node for e -> getC0 ();
	  if (e -> getc0 () != 0 ){
	    Node coef_vertex;
	    coef_vertex.node( coef_count, e -> getc0(), e -> getc0() , e -> getc0(), -2, 0 );
	    node_info.push_back(coef_vertex);

	    //printf ("Add coef vertex to graph (coef value   %f) \n", e -> getc0 () );
	    //printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), coef_count);
	    nauty_info->addElement((*i) -> Index (),  coef_count);
	    nauty_info->addElement( coef_count, (*i) -> Index ());


	    coef_count ++;
	  }

	  // for each term add nodes for their non-one coefficients and their variable

	  for (exprGroup::lincoeff::iterator el = e ->lcoeff().begin (); el != e -> lcoeff ().end (); ++el) {

	    if ( el -> second ==1){
	      //printf (" add edge index %d ,  index %d\n", (*i) -> Index (), el -> first -> Index()    );
	      nauty_info->addElement((*i) -> Index (),  el -> first -> Index());
	      nauty_info->addElement( el -> first -> Index (), (*i) -> Index ());
	    }
	    else{
	      //printf (" add new vertex to graph, coef # %d with coef %f \n", coef_count, el -> second);
	      Node coef_vertex;
	      coef_vertex.node( coef_count, el -> second, el -> second, el -> second, -2, 0 );
	      node_info.push_back(coef_vertex);

	      //printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), coef_count);
	      nauty_info->addElement((*i) -> Index (),  coef_count);
	      nauty_info->addElement( coef_count, (*i) -> Index ());

	      // printf (" add edge coef index %d ,  2nd index %d\n", coef_count,  el -> first -> Index()  );
	      nauty_info->addElement(coef_count,  el -> first -> Index());
	      nauty_info->addElement( el -> first -> Index (), coef_count);
	      coef_count ++;
	    }
	    // coefficient = el -> second

	    // variable index is el -> first -> Index ()
	  }

	}
	
      }
      else if ((*i) -> Image () -> Type () == UNARY) {
	//	printf ("variable is unary  %d\n", (*i) -> Index ());
	expression *arg = (*i) -> Image () -> Argument () ;
	nauty_info->addElement( arg-> Index(), (*i) -> Index() );
	//printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), arg-> Index()); 
      }
      else if ((*i) -> Image () -> Type () == AUX) {
	//printf ("variable is AUX  %d\n", (*i) -> Index ());
	nauty_info->addElement((*i) -> Index (), (*i) -> Image() -> Index());
	nauty_info->addElement( (*i) -> Image() -> Index(), (*i) -> Index() );
	//printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), (*i) -> Image() -> Index()); 
      }
      else if ((*i) -> Image () -> Type () == VAR) {
	//printf ("variable is VAR  %d, image %d \n", (*i) -> Index (), (*i) -> Image() -> Index());
	nauty_info->addElement((*i) -> Index (), (*i) -> Image() -> Index());
	nauty_info->addElement( (*i) -> Image() -> Index(), (*i) -> Index() );

	//printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), (*i) -> Image() -> Index()); 
      }
    }
    else {
      // printf ("variable is %d\n", (*i) -> Index ());
      Node var_vertex;

      // Bounds of +- infinity make the compare function likely to return a false negative. Rather than add inf as a boud, I use lb-1 (or ub +1 
      if( (*i) -> ub() >= COUENNE_INFINITY && (*i) -> lb() <= - COUENNE_INFINITY){
	var_vertex.node( (*i) -> Index () , 0 , 1 , 0 ,  -1, -1 );
	node_info.push_back(var_vertex);
	//printf( "var info index %d, lb %f, ub %f \n",(*i) -> Index () , 1 , 0 ) ; 
      }
      else  if( (*i) -> ub() >= COUENNE_INFINITY ){
	var_vertex.node( (*i) -> Index () , 0 , (*i) -> lb () , (*i) -> lb() -1 ,  -1, -1 );
	node_info.push_back(var_vertex);
	//printf( "var info index %d, lb %f, ub %f \n",(*i) -> Index () , (*i) -> lb () , (*i) -> lb () -1 ) ; 
      }
      else  if( (*i) -> lb() <= - COUENNE_INFINITY){
	var_vertex.node( (*i) -> Index () , 0 , (*i) -> ub () +1 , (*i) -> ub () ,  -1, -1 );
	node_info.push_back(var_vertex);
	//printf( "var info index %d, lb %f, ub %f \n",(*i) -> Index () , (*i) -> ub () +1 , (*i) -> ub () ) ; 
      }
      else{
	var_vertex.node( (*i) -> Index () , 0 , (*i) -> lb () , (*i) -> ub () ,  -1, -1 );
	//printf( "var info index %d, lb %f, ub %f \n",(*i) -> Index () , (*i) -> lb () , (*i) -> ub () ) ; 
	// var_vertex.get_index() , var_vertex.get_coeff() , var_vertex.get_lb() , var_vertex.get_ub() ,  var_vertex.get_code() );
	node_info.push_back(var_vertex);
	// this is an original variable
      }
    }
  }
  
}


void CouenneProblem::Compute_Symmetry() const{

  ChangeBounds (Lb (), Ub (), nVars ());

  // jnlst_ -> Printf(Ipopt::J_VECTOR, J_BRANCHING,"== Computing Symmetry\n");
  // for (int i = 0; i < nVars (); i++)
  //   if (Var (i) -> Multiplicity () > 0)
  //     jnlst_->Printf(Ipopt::J_VECTOR, J_BRANCHING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
  // 		     X  (i), Lb (i), Ub (i));
  // jnlst_->Printf(Ipopt::J_VECTOR, J_BRANCHING,"=============================\n");

  std::sort(node_info. begin (), node_info. end (), node_sort);

  for (std::vector <Node>:: iterator i = node_info. begin ();  i != node_info. end (); ++i) 
    (*i).color_vertex(-1);
  
  int color = 1;
  for (std::vector <Node>:: iterator i = node_info. begin (); i != node_info. end (); ++i) {
    if( (*i).get_color() == -1){
      (*i).color_vertex(color);
      //printf ("Graph vertex %d is given color %d\n", (*i).get_index(), color);
      nauty_info -> color_node((*i).get_index(), color);
      for (std::vector <Node>:: iterator j = i+1; j != node_info. end (); ++j)
	if( compare( (*i) , (*j) ) ==1){
	  (*j).color_vertex(color);
	  nauty_info -> color_node((*j).get_index(),color);
	  //printf ("Graph vertex %d is given color %d, the same as vertex %d\n", (*j).get_index(), color, (*i).get_index());
	}
      //       else
      // j = node_info. end();
      color++;
    }
  }

  //Print_Orbits ();

  nauty_info -> computeAuto();
}

  
void CouenneProblem::Print_Orbits () const {

  //printf ("num gens = %d, num orbits = %d \n", nauty_info -> getNumGenerators(), nauty_info -> getNumOrbits() );

  std::vector<std::vector<int> > *new_orbits = nauty_info->getOrbits();

  //printf("There were %d orbits and %d generators\n",
  //nauty_info->getNumOrbits(),
  //nauty_info->getNumGenerators());

  int nNonTrivialOrbits = 0;

  for (unsigned int i = 0; i < new_orbits -> size(); i++) {
    if ((*new_orbits)[i].size() > 1) 
      nNonTrivialOrbits++;
    else continue;

    // int orbsize = (*new_orbits)[i].size();
    // printf( "Orbit %d [size: %d] [", i, orbsize);
    // copy ((*new_orbits)[i].begin(), (*new_orbits)[i].end(),
    // 	  std::ostream_iterator<int>(std::cout, " "));
    // printf("] \n");
  }

  printf ("%d non-trivial orbits\n", nNonTrivialOrbits);

#if 0
  if (nNonTrivialOrbits)
    for (int i=0; i< nVars (); i++) {

      std::vector< int > *branch_orbit = Find_Orbit (i);

      if (branch_orbit -> size () > 1) {
  	printf ("x%04d: ", i);

  	for (std::vector<int>::iterator it = branch_orbit -> begin (); it != branch_orbit -> end (); ++it) 
  	  printf ("%d ", *it);
  	printf ("\n");
      }
    }
#endif
  delete new_orbits;
}

std::vector<int> *CouenneProblem::Find_Orbit(int index) const{

  std::vector<int> *orbit = new std::vector <int>;
  int which_orbit = -1;
  std::vector<std::vector<int> > *new_orbits = nauty_info->getOrbits();

  for (unsigned int i = 0; i < new_orbits -> size(); i++) {
    for (unsigned int j = 0; j < (*new_orbits)[i].size(); j++) {
      //   for (std::vector <int>:: iterator j = new_orbits[i].begin(); new_orbits[i].end(); ++j){
      if( (*new_orbits)[i][j] ==  index)
	which_orbit = i;
    }
  }
  
  //  for (std::vector <int>:: iterator j = new_orbits[which_orbit].begin(); new_orbits[which_orbit].end(), ++j)
  for (unsigned int j = 0; j < (*new_orbits)[which_orbit].size(); j++) 
    orbit -> push_back ((*new_orbits)[which_orbit][j]);

  delete new_orbits;

  return orbit;
}


void CouenneProblem::ChangeBounds (const double * new_lb, const double * new_ub, int num_cols) const {
  assert (num_cols == nVars ()); // replaced Variables () . size () as Variables () is not a const method
  std::sort(node_info. begin (), node_info. end (), index_sort);

  for (int  i = 0; i < num_cols; i++) {
    //   printf("Var %d  lower bound: %f   upper bound %f \n", i, new_lb[i], new_ub[i]);
    
    assert (node_info[i].get_index () == i);
    node_info[i ].bounds ( new_lb[i] , new_ub[i] );
    // printf("Var %d  INPUT lower bound: %f   upper bound %f \n", i, node_info[i].get_lb(), node_info[i].get_ub());
  }
}
#endif

void CouenneProblem::setupSymmetry () {

#ifdef COIN_HAS_NTY
  sym_setup ();
  Compute_Symmetry ();
  Print_Orbits ();
#else
  if (orbitalBranching_) {
    printf ("\
Couenne: Warning, you have set orbital_branching but Nauty is not available.\n\
Reconfigure with appropriate options --with-nauty-lib=/path/to/libnauty.* and --with-nauty-incdir=/path/to/nauty/include/files/\n");
  }
#endif
}
