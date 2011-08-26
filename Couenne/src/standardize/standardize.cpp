/* $Id$
 *
 * Name:    standardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize all expressions in a problem
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"

#include "BonBabSetupBase.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprIVar.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneDepGraph.hpp"

using namespace Ipopt;
using namespace Couenne;

void replace (CouenneProblem *p, int wind, int xind);

/// standardize (nonlinear) common expressions, objectives, and constraints

bool CouenneProblem::standardize () {

  if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
    printf ("Reformulation. current point: %d vars -------------------\n", 
	    domain_.current () -> Dimension ());
    for (int i=0; i<domain_.current () -> Dimension (); i++)
      printf ("%3d %20g [%20g %20g]\n", i, domain_.x (i), domain_.lb (i), domain_.ub (i));
  }

  bool retval = true;

  // create dependence graph to assign an order to the evaluation (and
  // bound propagation, and -- in reverse direction -- implied bounds)
  graph_ = new DepGraph;

  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); ++i)
    graph_ -> insert (*i);

  // allocate space in auxiliaries_ from commonexprs_

  //int initVar = variables_ . size () - commonexprs_ . size ();

  //std::set <int> DefVarDiffSet;

  // // DEFINED VARIABLES -----------------------------------------------------------------------

  // // standardize initial aux variables (aka defined variables, aka
  // // common expression)

  // if (commonexprs_.size ()) jnlst_ -> Printf (J_ALL, J_REFORMULATE,
  // 					      "%d common exprs, initVar = %d = %d - %d\n", 
  // 					      commonexprs_ . size (), 
  // 					      initVar, 
  // 					      variables_   . size (), 
  // 					      commonexprs_ . size ());

  // for (std::vector <expression *>::iterator i = commonexprs_ . begin ();
  //      i != commonexprs_ . end (); ++i) {

  //   if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
  //     printf ("\n=====> [1] stdz common expr [%d] := ", initVar); fflush (stdout);
  //     (*i) -> print (); printf ("\n");
  //   }

  //   exprAux    *naux = (*i) -> standardize (this, false);
  //   expression *img  = (*i);

  //   if (naux)
  //     img = naux -> Image ();

  //   if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
  //     printf ("\n=====> [2] "); fflush (stdout);
  //     if (*i)   (*i) -> print (); else printf ("(null)"); printf (" [*i] ");   fflush (stdout);
  //     if (naux) naux -> print (); else printf ("(null)"); printf (" [naux] "); fflush (stdout);
  //     if (img)  img  -> print (); else printf ("(null)"); printf (" [img]\n");
  //   }

  //   // trick to obtain same index as common expression: create exprAux
  //   // with index initVar and replace occurrences with address of
  //   // newly created exprAux through auxiliarize()

  //   exprVar *newvar;
  //   //img -> isInteger () ? exprAux::Integer : exprAux::Continuous);

  //   //auxiliarize (newvar); // puts newvar at right position in variables_

  //   if (((*i) -> Type () == VAR) ||
  // 	((*i) -> Type () == AUX)) {

  //     newvar = new exprAux (img, initVar, 1 + img -> rank (), exprAux::Unset, &domain_);
  //     replace (this, initVar, img -> Index ());
  //     //auxiliarize (variables_ [initVar], 
  //     //variables_ [img -> Index ()]);

  //     //delete variables_ [initVar];

  //     variables_ [initVar] = newvar;
  //     variables_ [initVar] -> zeroMult ();

  //   } else {

  //     if (img -> dependsOn (&initVar, 1, TAG_AND_RECURSIVE)) {

  // 	//printf ("depends! "); img -> print (); 

  // 	expression *diff = new exprSub (new exprClone (variables_ [initVar]), img);

  // 	//printf ("; diff = "); diff -> print ();

  // 	exprAux *diffAux = diff -> standardize (this, false);

  // 	//printf ("; aux: "); if (diffAux) diffAux -> print (); 

  // 	//if (diffAux)
  // 	exprAux *newAux = addAuxiliary (diff);

  // 	//printf ("; real aux: "); if (newAux) newAux -> print (); putchar ('\n');
	
  // 	//Lb (newAux -> Index ()) = 
  // 	//Ub (newAux -> Index ()) = 0.;

  // 	DefVarDiffSet. insert (newAux -> Index ());

  //     } else {

  // 	newvar = new exprAux (img, initVar, 1 + img -> rank (), exprAux::Unset, &domain_);
  // 	//replace (this, initVar, newvar -> Index ());

  // 	auxiliarize (newvar);

  // 	//delete variables_ [initVar];
  // 	variables_ [initVar] = newvar;

  // 	graph_ -> insert (newvar);
  // 	//if (naux) 
  // 	graph_ -> erase (naux);
  //     }
  //   }

  //   //variables_ . erase (variables_ . end () - 1);

  //   if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
  //     if (naux) {
  // 	printf ("done: "); fflush (stdout);
  // 	naux -> print ();
  // 	printf (" := "); fflush (stdout);
  // 	naux -> Image () -> print (); printf ("\n..."); fflush (stdout);
  //     } else if (*i) {
  // 	(*i) -> print ();
  // 	//printf (" := "); fflush (stdout);
  // 	//aux -> Image () -> print (); 
  // 	printf ("\n");
  //     } else printf ("[n]aux NULL!\n");
  //   }

  //   initVar++;
  // }

  // OBJECTIVES ------------------------------------------------------------------------------

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) {

    if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
      printf ("Objective [code: %d]", (*i) -> Body () -> code ()); 
      (*i) -> print ();
    }

    // check for constant objective

    std::set <int> deplist;

    if (0 == (*i) -> Body () -> DepList (deplist,  TAG_AND_RECURSIVE)) {

      // This objective is constant. Store the value and use it in
      // CouenneSolverInterface.

      constObjVal_ = (*i) -> Body () -> Value ();

    } else {

      exprAux *aux = (*i) -> standardize (this);

      if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
	printf (" objective "); (*i) -> print (); 
	if (aux) {printf (" admits aux "); aux -> print ();}
      }

      if (aux) {
	//delete ((*i) -> Body ());
	(*i) -> Body (new exprClone (aux));
      }

      if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
	printf (". New obj: "); (*i) -> print (); printf ("\n");
      }
    }
  }

  // Constraints ----------------------------------------------------------------------------

  // commuted_ is an array with a flag for each original variable,
  // which is true at position i if initially original variable x_i
  // became auxiliary (because of constraint 
  // 
  // x_{k+1} + f(x_1,x_2...,x_k} <=/=/>= 0
  //
  // becoming
  //
  // x_{k+1} <=/=/>= - f(x_1,x_2...,x_k}

  commuted_ = new bool [nVars ()];
  CoinFillN (commuted_, nVars (), false);

  std::vector <std::vector <CouenneConstraint *>::iterator> iters2erase;

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin (); 
       i != constraints_.end (); ++i) {

    if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
      printf ("\nReformulating constraint: "); 
      (*i) -> print ();
    }

    CouNumber 
      conLb = (*((*i) -> Lb ())) (),
      conUb = (*((*i) -> Ub ())) ();

    // sanity check: if (due to bad model or redundancies) the
    // constraint's body is constant, delete it -- check if it's
    // within bounds.

    expression *eBody = (*i) -> Body ();

    if (eBody -> Linearity () <= CONSTANT) {

      CouNumber bodyVal = (*eBody)();

      if ((bodyVal < conLb - COUENNE_BOUND_PREC) ||
	  (bodyVal > conUb + COUENNE_BOUND_PREC)) { // all variables eliminated, but out of bounds
	
	jnlst_ -> Printf (J_SUMMARY, J_PROBLEM, 
			  "Constraint: all variables eliminated, but value %g out of bounds [%g,%g]: ", 			  
			  bodyVal, conLb, conUb);

	if (jnlst_ -> ProduceOutput (J_SUMMARY, J_PROBLEM))
	  (*i) -> print ();

	retval = false;
	break;

      } else {

	iters2erase.push_back (i);
	continue; // all variables eliminated and constraint is redundant
      }
    }

    exprAux *aux = (*i) -> standardize (this);

    if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
      printf (" reformulated: aux w[%d] ", aux ? (aux -> Index ()) : -2); 
      (*i) -> print ();
    }

    if (aux) { // save if standardized

      // this is a top level auxiliary, i.e., an auxiliary whose node
      // in the DAG stands at the maximum level -- no other variable
      // depends on it as it is the lhs of a constraint.

      aux -> top_level () = true;

      // this constraint f(x)<=b is therefore replaced by w<=b, where
      // w is the top-level auxiliary associated with f(x)

      //printf ("delete %x: ", ((*i) -> Body ())); ((*i) -> Body ()) -> print ();
      //printf ("\n"); 
      //delete ((*i) -> Body ());

      (*i) -> Body (new exprClone (aux));

      //      con2.push_back (*i);
    }
    else {

      // left-hand side not reformulated, therefore this is probably
      // an affine expression

      //CouNumber lb, ub;

      //printf ("let's see what's left of the body: "); fflush (stdout);
      //(*i) -> Body () -> print (); printf ("\n");

      // (*i) -> Body () -> getBounds (lb, ub);

      // if ((((*((*i) -> Lb ())) ()) > ub + COUENNE_EPS) ||
      // 	  (((*((*i) -> Ub ())) ()) < lb - COUENNE_EPS)) {

      // 	jnlst_ -> Printf (J_SUMMARY, J_PROBLEM, "found infeasible constraint [%g,%g] vs [%g,%g]\n", 
      // 			  lb, ub,
      // 			  ((*((*i) -> Lb ())) ()),
      // 			  ((*((*i) -> Ub ())) ()));

      // 	if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE))
      // 	  (*i) -> print ();

      // 	retval = false;
      // }

      // delete constraint as no aux associated with it. It is a
      // linear constraint now (after reformulation), so it will be
      // inserted as a linearization of the aux's definition (an
      // affine function).

      iters2erase.push_back (i);
    }

    //(*i) -> Body () -> realign (this);

    if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
      printf (" --> "); (*i) -> print (); printf ("\n\n");
    }

    /*printf ("=== "); fflush (stdout); 
    aux -> print (); printf (" := "); fflush (stdout);
    aux -> Image () -> print (); printf ("\n");*/
  }

  for (unsigned int i = iters2erase.size (); i--;)
    constraints_. erase (iters2erase [i]);

  if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) { 
    // Use with caution. Bounds on auxs are not defined yet, so valgrind complains
    printf ("done with standardization: (careful, bounds below can be nonsense)\n");
    print (); 
  }

  delete auxSet_;

  // Create evaluation order ////////////////////////////////////////////////////

  // reallocate space for enlarged set of variables
  domain_.current () -> resize (nVars ());

  //graph_ -> print ();
  graph_ -> createOrder ();
  //graph_ -> print ();

  assert (graph_ -> checkCycles () == false);

  // Fill numbering structure /////////////////////////////////////////////////

  int n = nVars ();
  numbering_ = new int [n];
  std::set <DepNode *, compNode> vertices = graph_ -> Vertices ();

  for (std::set <DepNode *, compNode>::iterator i = vertices.begin ();
       i != vertices.end (); ++i)

    numbering_ [(*i) -> Order ()] = (*i) -> Index (); 

  //////////////////////////////////////////////////////////////////////////////

  // do initial bound propagation

  for (int i = 0; i < nVars (); i++) {

    int ord = numbering_ [i];

    if (variables_ [ord] -> Type () == AUX) {

      // initial auxiliary bounds are infinite (they are later changed
      // through branching)

      if (variables_ [ord] -> Index () >= nOrigVars_) { // and one that was not an original, originally...

	domain_.lb (ord) = -COIN_DBL_MAX;
	domain_.ub (ord) =  COIN_DBL_MAX;
      }

      //printf ("from "); variables_ [ord] -> Lb    () -> print (); 

      // tighten them with propagated bounds
      variables_ [ord] -> crossBounds ();

      //printf ("to "); variables_ [ord] -> Lb    () -> print (); printf (", now eval\n");

      if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
	printf (":::: %3d %10g [%10g, %10g]", 
		ord, domain_.x (ord), domain_.lb (ord), domain_.ub (ord));
      }

      // and evaluate them
      domain_.x  (ord) = (*(variables_ [ord] -> Image ())) ();
      domain_.lb (ord) = (*(variables_ [ord] -> Lb    ())) ();
      domain_.ub (ord) = (*(variables_ [ord] -> Ub    ())) ();

      // if (DefVarDiffSet.find (ord) != DefVarDiffSet.end ()) {

      // 	domain_.lb (ord) =
      // 	domain_.ub (ord) = 0.;
      // }

      if (jnlst_ -> ProduceOutput (J_ALL, J_REFORMULATE)) {
	printf (" --> %10g [%10g, %10g] [", 
		domain_.x (ord), domain_.lb (ord), domain_.ub (ord));
	variables_ [ord] -> Lb () -> print (); printf (",");
	variables_ [ord] -> Ub () -> print (); printf ("]\n");
      }

      bool integer = variables_ [ord] -> isInteger ();

      if (integer) {
	domain_.lb (ord) = ceil  (domain_.lb (ord) - COUENNE_EPS);
	domain_.ub (ord) = floor (domain_.ub (ord) + COUENNE_EPS);
      }
    }
  }

  // TODO: resolve duplicate index in exprQuad before restoring this

  // remove duplicates

  std::string delete_redund;

  if (bonBase_) bonBase_ -> options () -> GetStringValue ("delete_redundant", delete_redund, "couenne."); 
  else  delete_redund = "yes";

  if (delete_redund == "yes")

    // Look for auxiliaries of the form w:=x and replace each occurrence of w with x
    for (std::vector <exprVar *>::iterator i = variables_.begin (); 
       i != variables_.end (); ++i)

    if (((*i) -> Type () == AUX) && ((*i) -> sign () == expression::AUX_EQ)) {

      int type = (*i) -> Image () -> Type ();

      if ((type == VAR) || (type == AUX)) {

	// found w_k = x_h. 
	// 
	// Check if either is integer, the survivor will be integer too
	// Replace all occurrences of w_k with x_h

	/*printf ("redundancy: "); 
	(*i)             -> print (); printf (" := "); 
	(*i) -> Image () -> print (); printf ("\n");*/

	// use the info on the variable to be discarded: tighter
	// bounds and integrality that the replacement might not have.

	int 
	  indStays  = (*i) -> Image () -> Index (), // index h
	  indLeaves = (*i)             -> Index (); // index k

	if (indStays == indLeaves)  // very strange case, w_h = x_h
	  continue;

	// do not swap! x_h could be in turn an auxiliary...
	//
	//if (indStays > indLeaves) 
	//{int swap = indStays; indStays = indLeaves; indLeaves = swap;} // swap

	exprVar 
	  *varStays  = variables_ [indStays],
	  *varLeaves = variables_ [indLeaves];

	// intersect features of the two variables (integrality, bounds)

	varStays -> lb () = varLeaves -> lb () = CoinMax (varStays -> lb (), varLeaves -> lb ());
	varStays -> ub () = varLeaves -> ub () = CoinMin (varStays -> ub (), varLeaves -> ub ());

	if (varStays  -> isInteger () ||
	    varLeaves -> isInteger ()) {

	  varStays -> lb () = ceil  (varStays -> lb ());
	  varStays -> ub () = floor (varStays -> ub ());

	  if (varStays -> Type () == AUX)
	    varStays -> setInteger (true);
	  else {
	    //expression *old = varStays; // !!! leak
	    variables_ [indStays] = varStays = new exprIVar (indStays, &domain_);
	    auxiliarize (varStays); // replace it everywhere in the problem
	    //delete old;
	  }
	}

	auxiliarize (varLeaves, varStays); // now replace occurrences of w_k with x_h

	//if (varLeaves -> Index () >= nOrigVars_) // why check? It's not there anymore.
	varLeaves -> zeroMult (); // disable this variable
      }
    }

  /// re-check integrality. This is necessary as the initial
  /// integrality check is done on some continuous variables, which
  /// may turn out to be identical to other, integer, variables. See
  /// minlplib/ex1223.nl, where x_29 = x_4^2 and x_4=x_9, with x_4
  /// declared continuous and x_9 integer
  for (int ii=0; ii < nVars (); ii++) {

    int i = numbering_ [ii];

    if ((Var (i) -> Multiplicity () > 0) &&
	(Var (i) -> Type () == AUX) &&
	(Var (i) -> Image () -> isInteger ()) &&
	(Var (i) -> sign () == expression::AUX_EQ))
      Var (i) -> setInteger (true);

    //if (Var (i) -> Multiplicity () == 0)
    //Lb (i) = Ub (i) = 0.;
  }

  // check how many multiplications there are 

//   int nmul = 0;
//   // Look for auxiliaries of the form w:=x and replace each occurrence of w with x
//   for (std::vector <exprVar *>::iterator i = variables_.begin (); 
//        i != variables_.end (); ++i) {

//     if ((*i) -> Type () != AUX ||
// 	(*i) -> Multiplicity () <= 0) 
//       continue;

//     expression *img = (*i) -> Image ();

//     if (img -> code () != COU_EXPRMUL) 
//       continue;

//     expression **args = img -> ArgList ();

//     if ((args [0] -> Type () == AUX ||
// 	 args [0] -> Type () == VAR) &&
// 	(args [1] -> Type () == AUX ||
// 	 args [1] -> Type () == VAR))
//       nmul++;
//   }

//   printf ("\nMULS: %d/%d\n", nmul, variables_.size());

//   exit (-1);

  // TODO: re-compute ranks

  delete [] commuted_;  commuted_ = NULL;
  delete    graph_;     graph_    = NULL;

  return retval;
}
