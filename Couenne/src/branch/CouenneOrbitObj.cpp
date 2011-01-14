/* $Id$
 *
 * Name:    CouenneOrbitObj.cpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 * Purpose: Base object for variables (to be used in branching)
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

/*
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

#include "CouenneProblem.hpp"
#include "CouenneOrbitObj.hpp"
#include "CouenneBranchingObject.hpp"

#include "CouenneExprGroup.hpp"
using namespace Couenne;

void Node::node(int i, double c , double l, double u, int cod){
  index = i;
  coeff = c;
  lb = l;
  ub = u;
  color = -1;
  code = cod;
}
void Node::color_vertex(int k){
  color = k;
  }
 
bool compare (Node a, Node b){
  if(a.get_code() == b.get_code() )
    if(a.get_coeff() == b.get_coeff() )
      if(a.get_lb() == b.get_lb())
	if(a.get_ub() == b.get_ub())
	    return 1;
  return 0;   
  }

bool node_sort (Node a, Node b){
  bool is_less = 0;
  
  if(a.get_code() < b.get_code() )
    is_less = 1;
  else if(a.get_code() == b.get_code() )
      if(a.get_coeff() < b.get_coeff() )
	is_less = 1;
      else if(a.get_coeff() ==  b.get_coeff() )
	  if(a.get_lb() < b.get_lb()) 
	    is_less = 1;
	  else if(a.get_lb() == b.get_lb())
	      if(a.get_ub() < b.get_ub())
		is_less = 1;
	      else if(a.get_index() < b.get_index())
		  is_less = 1;
  
  return is_less;   
}
bool index_sort (Node a, Node b){
  return (a.get_index() < b.get_index() );
}





/// Empty constructor
CouenneOrbitObj::CouenneOrbitObj ():

  CouenneObject () {}

CouenneOrbitObj::CouenneOrbitObj (CouenneCutGenerator *cutgen,
				  CouenneProblem *p, 
				  exprVar *ref, 
				  Bonmin::BabSetupBase *base, 
				  JnlstPtr jnlst):
  CouenneObject (cutgen, p, ref, base, jnlst){

//  // Find Coefficients

/// initialize nauty

  int num_affine = 0;

  for (std::vector <exprVar *>:: iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end (); ++i) {
    
    if ((*i) -> Type () == AUX) {
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

 int nc = num_affine + p-> nVars ();
 printf (" There are   %d  coefficient vertices in the graph \n", num_affine);
 printf (" Graph has    %d  vertices \n", nc);


 nauty_info = new Nauty(nc);

  // create graph

 int coef_count= p-> nVars ();
 for (std::vector <exprVar *>:: iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end (); ++i) {

    //    printf ("I have code %d \n",  (*i) ->  Image() -> code() );


    
    if ((*i) -> Type () == AUX) {
      printf ("aux is %d with code %d \n", (*i) -> Index (), (*i) -> Image () -> code() );
      // this is an auxiliary variable
      
      Node vertex;
      vertex.node( (*i) -> Index () , 0.0 , (*i) -> lb () , (*i) -> ub () ,  (*i) -> Image () -> code() );
      node_info.push_back( vertex);
      		
      // add node in nauty graph for its index, (*i) -> Index ()

      if ((*i) -> Image () -> Type () == N_ARY) {
	
	if ((*i) -> Image () -> code () != COU_EXPRGROUP) {
	  
	  for (int a=0; a < (*i) -> Image () -> nArgs(); a++) {
	    expression *arg = (*i) -> Image () -> ArgList () [a];
	    
	    if (arg -> Type () != CONST) {
	      printf (" add edge  %d , %d\n", (*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement((*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement( arg -> Index (), (*i) -> Index ());
	    }

	    else{
	      
	      assert (arg -> Type () == CONST);	

	      printf (" add new vertex to graph, coef # %d, value %g \n", coef_count, arg -> Value() );
	      printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (),  coef_count);
	      nauty_info->addElement((*i) -> Index (),  coef_count);
	      nauty_info->addElement( coef_count, (*i) -> Index ());
	      

	      Node coef_vertex;
	      coef_vertex.node( coef_count, arg -> Value(), arg -> Value() , arg -> Value(), -2 );
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
	    coef_vertex.node( coef_count, e -> getc0(), e -> getc0() , e -> getc0(), -2 );
	    node_info.push_back(coef_vertex);

	    printf ("Add coef vertex to graph (coef value   %f) \n", e -> getc0 () );
	    printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), coef_count);
	    nauty_info->addElement((*i) -> Index (),  coef_count);
	    nauty_info->addElement( coef_count, (*i) -> Index ());


	    coef_count ++;
	  }
	  
	  // for each term add nodes for their non-one coefficients and their variable

	  for (exprGroup::lincoeff::iterator el = e ->lcoeff().begin (); el != e -> lcoeff ().end (); ++el) {

	    if ( el -> second ==1){ 
	      printf (" add edge index %d ,  index %d\n", (*i) -> Index (), el -> first -> Index()    );
	      nauty_info->addElement((*i) -> Index (),  el -> first -> Index());
	      nauty_info->addElement( el -> first -> Index (), (*i) -> Index ());
	    }
	    else{
	      printf (" add new vertex to graph, coef # %d with coef %f \n", coef_count, el -> second);
	      Node coef_vertex;
	      coef_vertex.node( coef_count, el -> second, el -> second, el -> second, -2 );
	      node_info.push_back(coef_vertex);

	      printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (), coef_count);
	      nauty_info->addElement((*i) -> Index (),  coef_count);
	      nauty_info->addElement( coef_count, (*i) -> Index ());
	      
	      printf (" add edge coef index %d ,  2nd index %d\n", coef_count,  el -> first -> Index()  );
	      nauty_info->addElement(coef_count,  el -> first -> Index());
	      nauty_info->addElement( el -> first -> Index (), coef_count);
	      coef_count ++;
	    }
	    // coefficient = el -> second
	    
	    // variable index is el -> first -> Index ()
	  }
	  
	}

      } else if ((*i) -> Image () -> Type () == UNARY) {

      }

    } else {
      //  printf ("variable is %d\n", (*i) -> Index ());
      Node var_vertex;
      var_vertex.node( (*i) -> Index () , 0 , (*i) -> lb () , (*i) -> ub () ,  -1 );
      //      printf( "var info index %d, coef %f, lb %f, ub %f, code %d \n", var_vertex.get_index() , var_vertex.get_coeff() , var_vertex.get_lb() , var_vertex.get_ub() ,  var_vertex.get_code() );
      node_info.push_back(var_vertex);
       // this is an original variable
      
    }
  }

}


/// Constructor with lesser information, used for infeasibility only
CouenneOrbitObj::CouenneOrbitObj (exprVar *ref, 
				  Bonmin::BabSetupBase *base, 
				  JnlstPtr jnlst):

  CouenneObject (ref, base, jnlst) {}


/// Copy constructor
CouenneOrbitObj::CouenneOrbitObj (const CouenneOrbitObj &src):
  CouenneObject       (src) {}


/// apply the branching rule
OsiBranchingObject *CouenneOrbitObj::createBranch (OsiSolverInterface *si,
						   const OsiBranchingInformation *info,
						   int way) const {

  return NULL;
}





// set point at current LP solution
double CouenneOrbitObj::feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const {
  return 0;
}



/// non-linear infeasibility -- not called by independent's CouenneVarObject
double CouenneOrbitObj::infeasibility (const OsiBranchingInformation *info, int &way) const {

  return 0;
}


/// non-linear infeasibility -- no need for the domain's push
/// instruction as this is called from
/// CouenneVarObject::infeasibility()
double CouenneOrbitObj::checkInfeasibility (const OsiBranchingInformation *info) const {

  return 0;

}

void CouenneOrbitObj::Compute_Symmetry(){
  sort(node_info. begin (), node_info. end (), node_sort);
  
  int color = 1;
  for (std::vector <Node>:: iterator i = node_info. begin (); 
       i != node_info. end (); ++i) {
    if( (*i).get_color() == -1){
      (*i).color_vertex(color);
      printf ("Graph vertex %d is given color %d\n", (*i).get_index(), color);
      nauty_info -> color_node((*i).get_index(), color);
      for (std::vector <Node>:: iterator j = i+1; j <= node_info. end (); ++j) 
	if( compare( (*i) , (*j) ) ==1){
	  (*j).color_vertex(color);
	  nauty_info -> color_node((*j).get_index(),color);
	  printf ("Graph vertex %d is given color %d, the same as vertex %d\n", (*j).get_index(), color, (*i).get_index());
	}
      //       else
      // j = node_info. end();
      color++;
      
    }
  }
  
  nauty_info -> computeAuto();
}

void CouenneOrbitObj::Print_Orbits(){

  printf("num gens = %d, num orbits = %d \n", nauty_info -> getNumGenerators(), nauty_info -> getNumOrbits() );
  
  std::vector<std::vector<int> > new_orbits = nauty_info->getOrbits();
  
  printf("There were %d orbits and %d generators\n",
	 nauty_info->getNumOrbits(),
	 nauty_info->getNumGenerators());
  
  for (unsigned int i = 0; i < new_orbits.size(); i++) {
    printf( "Orbit %d [", i);
    copy(new_orbits[i].begin(), new_orbits[i].end(),
	 std::ostream_iterator<int>(std::cout, " "));
    printf("] \n");
  }
  
  
  
}


  
void CouenneOrbitObj::ChangeBounds (CouenneProblem  * p ){
  sort(node_info. begin (), node_info. end (), index_sort);
 
  for (std::vector <exprVar *>:: iterator i =  p->Variables (). begin (); 
       i !=  p->Variables (). end (); ++i) {
    node_info[(*i)->Index() ].bounds ( (*i) -> lb () , (*i) -> ub () );

  }

  
}
*/
