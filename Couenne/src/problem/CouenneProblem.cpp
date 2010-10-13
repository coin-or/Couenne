/* $Id$
 *
 * Name:    CouenneProblem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CbcBranchActual.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprIVar.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprOpp.hpp"

#include "CouenneProblem.hpp"
#include "CouenneGlobalCutOff.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneLQelems.hpp"

#include "CouenneObject.hpp"

using namespace Couenne;

/// methods to add objective function. 

void CouenneProblem::addObjective (expression *newobj, const std::string &sense) {
  objectives_ . push_back
    (new CouenneObjective ((sense == "min") ? 
			   newobj : new exprOpp (new exprClone (newobj))));
}


/// methods to add nonlinear constraints:

/// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs) {

  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

/// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (COUENNE_INFINITY)));
}

/// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (-COUENNE_INFINITY), rhs));
}

/// Add (non linear) objective function
void CouenneProblem::setObjective (int indObj, expression * newObj, const std::string &sense) {
  objectives_ [indObj] = (new CouenneObjective ((sense == "min") ? 
						newObj : new exprOpp (new exprClone (newObj))));
}


/// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb, expression *ub) {
  if (!lb) lb = new exprConst (0.);
  if (!ub) ub = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}



/// add variable to the problem -- check whether it is integer (isDiscrete)

expression *CouenneProblem::addVariable (bool isDiscrete, Domain *d) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size (), d)) :
    (new exprVar  (variables_ . size (), d));

  variables_ . push_back (var);

  if (isDiscrete) 
    nIntVars_++;

  nOrigVars_++;

  return var;
}


/// add auxiliary variable and associate it with pointer to expression
/// given as argument
exprAux *CouenneProblem::addAuxiliary (expression *symbolic) {

  // check if image is already in the expression database auxSet_
  std::set <exprAux *, compExpr>::iterator i;

  int var_ind = variables_ . size ();
  domain_. current () -> resize (var_ind + 1);

  symbolic -> getBounds (domain_. lb (var_ind), 
			 domain_. ub (var_ind));

  // create new aux associated with that expression
  exprAux *w = new exprAux (symbolic,
			    var_ind,
			    1 + symbolic -> rank (), 
			    exprAux::Unset, 
			    &domain_);
  //symbolic -> isInteger () ? exprAux::Integer : exprAux::Continuous);

  //  w -> linkDomain (&domain_);

  // seek expression in the set
  if ((i = auxSet_ -> find (w)) == auxSet_ -> end ()) {

    // no such expression found in the set, create entry therein
    variables_ . push_back (w);
    auxSet_ -> insert (w); // insert into repetition checking structure
    graph_  -> insert (w); // insert into acyclic structure

  } else {  // otherwise, just return the entry's pointer

    w->Image(NULL); // otherwise "delete w" will also delete user given expression "symbolic"
    delete w;
    w = *i;
    (*i) -> increaseMult ();
  }

  return w;
}


/// translates pair (indices, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexL, 
				    CouNumber *coeff,
				    std::vector <std::pair <exprVar *, CouNumber> > &lcoeff) {
  // TODO: sort

  for (int i=0; indexL [i] >= 0; i++)
    lcoeff.push_back (std::pair <exprVar *, CouNumber> (Var (indexL [i]), coeff [i]));
}


/// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexI,
				    int *indexJ,
				    CouNumber *coeff,
				    std::vector <quadElem> &qcoeff) {
  // TODO: sort

  for (int i=0; indexI [i] >= 0; i++)
    qcoeff.push_back (quadElem (Var (indexI [i]), Var (indexJ [i]), coeff [i]));
}


/// fill in the integerRank_ array
void CouenneProblem::fillIntegerRank () const {

  if (integerRank_)
    return;

  int nvars = nVars ();

  integerRank_ = new int [nvars];

  // 0: fractional
  // 1: integer
  // k: integer,    depending on at least one integer with associated value k-1, or
  // k: fractional, depending on at least one integer with associated value k

  for (int ii = 0; ii < nvars; ii++) {

    int index = numbering_ [ii];

    if (Var (index) -> Multiplicity () <= 0) {
      integerRank_ [index] = 0;
      continue;
    }

    bool isInt = Var (index) -> isDefinedInteger ();

    integerRank_ [index] = (isInt) ? 1 : 0;

    if (Var (index) -> Type () == AUX) {

      std::set <int> deplist;

      if (Var (index) -> Image () -> DepList (deplist, STOP_AT_AUX) != 0) // depends on something
	for (std::set <int>::iterator i = deplist.begin (); i != deplist.end (); ++i) {

	  int token = integerRank_ [*i];
	  if (isInt) token++;

	  if (token > integerRank_ [index]) // there's a free integer below us
	    integerRank_ [index] = token;
	}
    }
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "Free (original) integers\n");
  for (int i=0; i<nOrigVars_; i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, integerRank_ [i]);

  // fill in numberInRank_
  for (int i=0; i<nOrigVars_; i++)
    if ((variables_ [i] -> isDefinedInteger ()) &&
	(variables_ [i] -> Multiplicity () > 0)) {

    int rank = integerRank_ [i];

    if (numberInRank_.size () <= (unsigned int) rank)
      for (int j=numberInRank_.size (); j <= rank; j++)
	numberInRank_ .push_back (0);

    numberInRank_ [rank] ++;
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "numInteger [neglect non-originals]\n");
  for (unsigned int i=0; i<numberInRank_.size(); i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, numberInRank_ [i]);
}


/// Called from simulateBranch when object is not CouenneObject and
/// therefore needs explicit FBBT
bool Couenne::BranchingFBBT (CouenneProblem *problem,
			     OsiObject *Object,
			     OsiSolverInterface *solver) {

  // do not perform this if object is not a variable object

  bool feasible = true;

  if (problem -> doFBBT ()) {

    int 
      indVar = Object  -> columnNumber (),
      nvars  = problem -> nVars ();

    // do not perform this if object is not a variable object

    if (indVar >= 0) {

      t_chg_bounds *chg_bds = new t_chg_bounds [nvars];
      chg_bds [indVar].setUpper (t_chg_bounds::CHANGED);
      problem -> installCutOff ();

      if ((feasible = problem -> btCore (chg_bds))) {

	const double
	  *lb = solver -> getColLower (),
	  *ub = solver -> getColUpper ();
	  
	for (int i=0; i<nvars; i++) {
	  if (problem -> Lb (i) > lb [i]) solver -> setColLower (i, problem -> Lb (i));
	  if (problem -> Ub (i) < ub [i]) solver -> setColUpper (i, problem -> Ub (i));
	}
      }

      delete [] chg_bds;
    }
  }

  return feasible;
}


//Symmetry Stuff --------------------------------------------------------------------

#ifdef COIN_HAS_NTY
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

bool CouenneProblem::compare ( Node a, Node b) const{
  if(a.get_code() == b.get_code() )
    if(a.get_coeff() == b.get_coeff() )
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
      // printf ("aux is %d with code %d \n", (*i) -> Index (), (*i) -> Image () -> code() );
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
	      //printf (" add edge  %d , %d\n", (*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement((*i) -> Index (),  arg -> Index ());
	      nauty_info->addElement( arg -> Index (), (*i) -> Index ());
	    }

	    else{

	      assert (arg -> Type () == CONST);

	      //printf (" add new vertex to graph, coef # %d, value %g \n", coef_count, arg -> Value() );
	      //printf (" add edge aux index %d ,  coef index %d\n", (*i) -> Index (),  coef_count);
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

	    // printf ("Add coef vertex to graph (coef value   %f) \n", e -> getc0 () );
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
	      coef_vertex.node( coef_count, el -> second, el -> second, el -> second, -2 );
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


void CouenneProblem::Compute_Symmetry() const{
  //printf("Computing Symmetry\n");
  std::sort(node_info. begin (), node_info. end (), node_sort);

  for (std::vector <Node>:: iterator i = node_info. begin ();
       i != node_info. end (); ++i) 
    (*i).color_vertex(-1);
  
  int color = 1;
  for (std::vector <Node>:: iterator i = node_info. begin ();
       i != node_info. end (); ++i) {
    if( (*i).get_color() == -1){
      (*i).color_vertex(color);
      // printf ("Graph vertex %d is given color %d\n", (*i).get_index(), color);
      nauty_info -> color_node((*i).get_index(), color);
      for (std::vector <Node>:: iterator j = i+1; j <= node_info. end (); ++j)
	if( compare( (*i) , (*j) ) ==1){
	  (*j).color_vertex(color);
	  nauty_info -> color_node((*j).get_index(),color);
	  //	  printf ("Graph vertex %d is given color %d, the same as vertex %d\n", (*j).get_index(), color, (*i).get_index());
	}
      //       else
      // j = node_info. end();
      color++;

    }
  }

  nauty_info -> computeAuto();
}
  
void CouenneProblem::Print_Orbits(){

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
std::vector<int>  CouenneProblem::Find_Orbit(int index){

  std::vector<int> orbit;
  int which_orbit = -1;
  std::vector<std::vector<int> > new_orbits = nauty_info->getOrbits();

  for (unsigned int i = 0; i < new_orbits.size(); i++) {
    for (unsigned int j = 0; j < new_orbits[i].size(); j++) {
      //   for (std::vector <int>:: iterator j = new_orbits[i].begin(); new_orbits[i].end(); ++j){
      if( new_orbits[i][j] ==  index)
	which_orbit = i;
    }
  }
  
  //  for (std::vector <int>:: iterator j = new_orbits[which_orbit].begin(); new_orbits[which_orbit].end(), ++j)
  for (unsigned int j = 0; j < new_orbits[which_orbit].size(); j++) 
    orbit.push_back(new_orbits[which_orbit][j]);
    
  return orbit;
}

    


void CouenneProblem::ChangeBounds (const double * new_lb, const double * new_ub, int num_cols) const {
  assert (num_cols == Variables().size());
  std::sort(node_info. begin (), node_info. end (), index_sort);

  for (int  i = 0; i < num_cols; i++) {
    //   printf("Var %d  lower bound: %f   upper bound %f \n", i, new_lb[i], new_ub[i]);
    
    assert (node_info[i].index == i);
    node_info[i ].bounds ( new_lb[i] , new_ub[i] );
    // printf("Var %d  INPUT lower bound: %f   upper bound %f \n", i, node_info[i].get_lb(), node_info[i].get_ub());
  }
  

}

#endif



/// Get cutoff
CouNumber CouenneProblem::getCutOff () const
{return pcutoff_ -> getCutOff ();}


/// Get cutoff solution
CouNumber *CouenneProblem::getCutOffSol () const
{return pcutoff_ -> getCutOffSol ();}

