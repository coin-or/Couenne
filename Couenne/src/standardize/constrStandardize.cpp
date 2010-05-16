/* $Id$
 *
 * Name:    constrStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization of constraints
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblemElem.hpp"
#include "CouenneProblem.hpp"

#include "exprClone.hpp"
#include "exprAux.hpp"
#include "exprIVar.hpp"
#include "depGraph.hpp"

//#define DEBUG

// replace a variable
void replace (CouenneProblem *p, int wind, int xind);


/// decompose body of constraint through auxiliary variables
exprAux *CouenneConstraint::standardize (CouenneProblem *p) {

  // spot an auxiliary variable in constraint's body w - f(x) and move
  // the explicit w into the vector of auxiliary variables
  //
  // do it only if this is an equality constraint and there is at
  // least one variable that did not show up so far (need a problem
  // structure)

#ifdef DEBUG
  printf ("################################\nStandardizing constraint: "); print ();

  printf (" ["); fflush (stdout); lb_ -> print ();
  printf (","); fflush (stdout);  ub_ -> print (); fflush (stdout);
  /*  printf ("] {with auxset = ");
  for (std::set <exprAux *, compExpr>::iterator i = p -> AuxSet () -> begin ();
       i != p -> AuxSet () -> end (); i++) {
    printf ("<"); (*i) -> print (); 
    printf (","); (*i) -> Image () -> print (); printf ("> ");
    }*/
  printf ("]\n");
#endif

  if (compareExpr (&lb_, &ub_) == 0) {

    // this is an equality constraint, and could be the definition of
    // an auxiliary

    expression *rest;

    // split w from f(x)
    int wind = p -> splitAux ((*lb_) (), body_, rest, p -> Commuted ());

    if (wind >= 0) { // this IS the definition of an auxiliary variable w = f(x)

      // first, simplify expression (you never know)

      expression *restSimple = rest -> simplify ();

      if (restSimple) {
	delete rest;
	rest = restSimple;
      }

      // second, if this is a constraint of the form x=k, reset x's
      // bounds and do nothing else

      if (rest -> code () == COU_EXPRCONST) {

	p -> Var (wind) -> lb () = 
	p -> Var (wind) -> ub () = rest -> Value ();

	delete rest;
	return NULL;
      }

      // now assign a new auxiliary variable (with the same index,
      // "wind", as the old original variable) to the expression
      // contained in "rest"

      p -> Commuted () [wind] = true;

#ifdef DEBUG
      printf ("---> %d & ", wind); fflush (stdout);
      rest -> print (); printf ("\n");
#endif

      assert (p -> Var (wind) -> Type () == VAR);

      // check if constraint now reduces to w_k = x_h, and if so
      // replace all occurrences of x_k with x_h

      int xind = rest -> Index ();

      if (xind >= 0) {

	// create new variable, it has to be integer if original variable was integer
	exprAux *w = new exprAux (new exprClone (p -> Var (xind)), wind, 1 + p -> Var (xind) -> rank (),
				  p -> Var (wind) -> isInteger () ?
				  exprAux::Integer : exprAux::Continuous,
				  p -> domain ());
	p -> auxiliarize (w);
	w -> zeroMult ();
	
	replace (p, wind, xind);

	p -> auxiliarize (p -> Var (wind), p -> Var (xind));
	//p -> Var (wind) -> zeroMult (); // redundant variable is neutralized


	// 	//replace (p, wind, xind);

	// 	p -> auxiliarize (p -> Var (wind), p -> Var (xind));
	// 	p -> Var (wind) -> zeroMult (); // redundant variable is neutralized

      } else {

	// create new variable, it has to be integer if original variable was integer
	exprAux *w = new exprAux (rest, wind, 1 + rest -> rank (),
				  p -> Var (wind) -> isInteger () ? 
				  exprAux::Integer : exprAux::Continuous,
				  p -> domain ());

	std::set <exprAux *, compExpr>::iterator i = p -> AuxSet () -> find (w);

	// no such expression found in the set:
	if (i == p -> AuxSet () -> end ()) {

	  p -> AuxSet      () -> insert (w); // 1) beware of useless copies
	  p -> getDepGraph () -> insert (w); // 2) introduce it in acyclic structure

#ifdef DEBUG
	  printf ("now replacing x [%d] with ", wind); fflush (stdout);
	  w -> print (); printf (" := ");
	  w -> Image () -> print (); printf ("\n");
#endif

	  // replace ALL occurrences of original variable (with index
	  // wind) with newly created auxiliary
	  p -> auxiliarize (w);
	} 

	else {

#ifdef DEBUG
	  printf ("found aux occurrence of "); fflush (stdout);
	  w -> print (); printf (" := ");
	  w -> Image () -> print (); printf (" ... ");
	  (*i) -> print (); printf (" := ");
	  (*i) -> Image () -> print (); printf ("\n");
#endif

	  // if this is an original variable and is competing with an
	  // auxiliary variable, or at least this has a lower index,
	  // we better replace it throughout and eliminate that
	  // aux. See globallib/st_glmp_fp3 for an example where this
	  // otherwise would be a bug (x_1 unlinked from x2-x3 and
	  // leading to unbounded)
	  
	  int xind = (*i) -> Index (), iMax, iMin;

	  if (xind < wind) {
	    iMax = wind;
	    iMin = xind;
	  } else {
	    iMax = xind;
	    iMin = wind;
	  }

	  replace (p, iMax, iMin);

	  p -> auxiliarize (p -> Var (iMax), p -> Var (iMin));
	  p -> Var (iMax) -> zeroMult (); // redundant variable is neutralized
	  p -> auxiliarize (w);
	}
      }

      return NULL;
    }
  }

#ifdef DEBUG
  printf ("\nnormal\n-----------------\n");
#endif

  return body_ -> standardize (p);
}


// Replace a variable ////////////////////////////////
void replace (CouenneProblem *p, int wind, int xind) {

  exprVar 
    *varLeaves = p -> Variables () [wind],
    *varStays  = p -> Variables () [xind];

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
      p -> Variables () [xind] = varStays = new exprIVar (xind, p -> domain ());
      p -> auxiliarize (varStays); // replace it everywhere in the problem
      //delete old;
    }
  }
}
