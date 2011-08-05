/* $Id$
 *
 * Name:    auxiliarize.cpp
 * Author:  Pietro Belotti
 * Purpose: replace occurrences of original variable in a problem with
 *          auxiliary with the same index
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneDepGraph.hpp"

#include <cstdio>
#include <cassert>

using namespace Couenne;

/// replace, in all expressions of the problem (auxiliaries,
/// objectives and constraints) link to an original variable that has
/// gone auxiliary

void CouenneProblem::auxiliarize (exprVar *aux, exprVar *subst) {

  if (graph_ && subst && aux -> Index () != subst -> Index ())
    graph_ -> replaceIndex (aux -> Index (), subst -> Index ());

  if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
    printf ("replacing  "); if (aux)   aux   -> print (); 
    if (subst) {printf (" with "); subst -> print ();}
    printf ("\n");
  }

  bool same_var = (subst == NULL);

  if (!subst) 
    subst = aux;

  // find original variable with index = w -> Index ()

  int index = aux -> Index ();

  assert (index >= 0);

  std::vector <exprVar *>::iterator orig;

  for (orig  = variables_.begin ();
       orig != variables_.end (); ++orig)

    if ((((*orig) -> Type () == VAR) || !same_var) && 
	((*orig) -> Index () == index)) // found it

      break;

  if (orig == variables_ . end ()) {
    printf ("CouenneProblem::auxiliarize: no original variables correspond\n");
    return;
  }

  // all common expressions

  for (std::vector <expression *>::iterator i = commonexprs_.begin ();
       i != commonexprs_.end (); ++i) {

    expression *body = *i;

    if (body) {

      if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	printf ("replacing within common expression [%d]: ", i - commonexprs_.begin ()); fflush (stdout); (*i) -> print (); printf ("\n");
      }

      if ((body -> Type () == VAR) || 
	  (body -> Type () == AUX)) {

	if (body -> Index () == (*orig) -> Index ()) {
      
	  delete body;
	  *i = new exprClone (subst);
	}
      } else body -> replace (*orig, subst);
    }
  }

  // all objectives

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) {

    expression *body = (*i) -> Body ();

    if (body) {

      if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	printf ("replacing within objective: "); fflush (stdout); (*i) -> print (); 
      }

      if ((body -> Type () == VAR) || 
	  (body -> Type () == AUX)) {

	if (body -> Index () == (*orig) -> Index ()) {
      
	  delete body;//(*i) -> Body ();
	  (*i) -> Body (new exprClone (subst));
	}
      } else body -> replace (*orig, subst);
    }
  }

  // and all constraints

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i) {

    expression *body = (*i) -> Body ();

    if (body) {

      if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	printf ("replacing within constraint [%d]: ", i - constraints_.begin ()); fflush (stdout); (*i) -> print (); 
      }

      if ((body -> Type () == VAR) ||
	  (body -> Type () == AUX)) {

	if (body -> Index () == (*orig) -> Index ()) {
      
	  delete body;//(*i) -> Body ();
	  (*i) -> Body (new exprClone (subst));
	}
      } else body -> replace (*orig, subst);
    }
  }

  // substitute it with w in all auxiliaries

  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)

    if (((*i) -> Type () == AUX) &&                  // replace in all aux's image
	((*i) -> Multiplicity () > 0) &&             // this variable is actually used
	((*i) -> Index () != (*orig) -> Index ())) { // skip same variable

      if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	printf ("replacing aux "); fflush (stdout); (*i) -> print (); 
	printf (" := "); fflush (stdout); (*i) -> Image () -> print (); 
	printf ("\n");
      }

      expression *image = (*i) -> Image ();

      if ((image -> Type () == VAR) || 
	  (image -> Type () == AUX)) {

	if (image -> Index () == (*orig) -> Index ()) {

	  delete image;
	  (*i) -> Image (new exprClone (subst));
	  //printf (" ---> "); (*i) -> Image () -> print (); 
	} //else (*i) -> Image () -> replace (*orig, aux);
      } else image  -> replace (*orig, subst);

      //printf ("\n");
    }

  // replace it with new auxiliary

  if (same_var)
    *orig = aux;
}
