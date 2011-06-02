/* $Id$
 *
 * Name:    conv-exprMul-reformulate.cpp
 * Author:  Pietro Belotti
 * Purpose: utility methods to reformulate multiplications into bilinear and trilinear terms
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <queue>

#include "CouenneTypes.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprTrilinear.hpp"
#include "CouenneExprBMul.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

/// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {

  register int t1 = v1 -> Type (), t2;
  return (((t1 == VAR) || (t1 == AUX)) &&
	  (((t2 = v2 -> Type ()) == VAR) || (t2 == AUX)) && 
	  (v1 -> Index () == v2 -> Index ()));
}


/// Create standard formulation of this expression

exprAux *exprMul::standardize (CouenneProblem *p, bool addAux) {

  //printf ("standardizing exprMul (addaux=%d) ", addAux); fflush (stdout); print (); printf ("\n");

  CouNumber coeff = 1.;
  std::map <int, CouNumber> indCoe;

  p -> flattenMul (this, coeff, indCoe);

  int nArgs = 0;

  expression **arglist = new expression * [indCoe.size ()];

  for (std::map <int, CouNumber>::iterator i = indCoe.begin (); i != indCoe.end (); ++i)
    if (i -> second == 1.) arglist [nArgs++] = new exprClone (p                -> Var (i -> first));
    else                   arglist [nArgs++] = new exprPow   (new exprClone (p -> Var (i -> first)), new exprConst (i -> second));

  //while (nargs_--)
  //delete arglist_ [nargs_];

  //delete [] arglist_;

  arglist_ = arglist;
  nargs_ = (int) indCoe.size();

  //printf ("new mul [%d]: %g * ", nargs_, coeff); fflush (stdout); print (); printf (" -- ");

  exprOp::standardize (p);

  //printf ("standardized: "); fflush (stdout); print (); printf ("\n");

  if (nargs_ == 1) {

    if (coeff != 1.) {

      // leak: arglist_ should be deleted

      expression *simMul = new exprMul (new exprConst (coeff), new exprClone (arglist_ [0]));
      return (addAux ? (p -> addAuxiliary (simMul)) : new exprAux (simMul, p -> domain ()));

    } else {

      // quick fix for nvs08: expression x_0^.5 * x_0^3 is quite
      // unusual and is translated into an exprMul with a single
      // argument.

      exprAux *var = dynamic_cast <exprAux *> (arglist_ [0] -> Copy ());
      return var;
    }
  }

  //if (nargs_ == 1)  // TODO: what happens really?
      //return NULL;
  /* {
     exprAux *aux = arglist_ [0];
     arglist_ [0] = NULL;
     return aux;
     } */

  // enable this when binary products are available

#if 0
  // check if it is a product of binary variables
  bool isBinProd = true;
  for (int i=nargs_; i--;) {

    expression *arg = arglist_ [i];
    if (arg -> isInteger ()) {

      CouNumber lb, ub;
      arg -> getBounds (lb, ub);
      if ((fabs (ceil  (lb - COUENNE_EPS))      > 0.) ||
	  (fabs (floor (ub + COUENNE_EPS) - 1.) > 0.)) { // if not all conditions hold, 
	isBinProd = false;
	break;
      }
    } else {
      isBinProd = false;
      break;
    }
  }

  if (isBinProd) {
    //printf ("found a BinProd!\n");
  }
#endif

  enum Couenne::TrilinDecompType type = p -> getTrilinDecompType ();

  if (nargs_ <= 2)
    type = Couenne::rAI;

  expression *retExpr;

  switch (type) {

  case Couenne::treeDecomp: {

    //printf ("trying treeDecomp on "); print (); printf ("\n");

    // A more hierarchical decomposition. Example:
    //
    // (x1 x2 x3 x4 x5 x6 x7 x8 x9 x10)
    //
    // becomes
    //
    // ((((x1 x2) (x3 x4)) ((x5 x6) (x7 x8))) (x9 x10))
    //
    // so that x1 .. x10 are the leaves of a binary tree whose
    // non-leaf nodes are auxiliaries
    //
    // use a queue to parse the leaf level first (by pushing all
    // original members) and push the newly added auxs into the queue
    // while extracting _pairs_ of elements

    std::queue <expression *> queue;

    for (int i=0; i<nargs_; i++)
      queue.push (arglist_ [i]);

    expression *aux;

    while (queue.size() > 1) {

      expression *arg1 = queue.front (); queue.pop ();
      expression *arg2 = queue.front (); queue.pop ();

      //printf ("Coupling "); arg1 -> print (); 
      //printf (" with "); arg2 -> print (); 

      if (areSameVariables (arg1, arg2)) aux = new exprPow (new exprClone (arg1), new exprConst (2.));
      else                               aux = new exprMul (new exprClone (arg1), new exprClone (arg2));

      //printf (" --> "); aux -> print (); printf ("\n"); 

      if (!(queue.empty ()))
	aux = p -> addAuxiliary (aux);

      queue.push (aux);
    }

    aux = queue.front (); queue.pop ();

    retExpr = aux; //(addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
  } break;

  // ----------------------------------------------------------------------------------------------

  case Couenne::tri_bi:
  case Couenne::bi_tri: { // the two cases are very similar

    //printf ("trying %s on ", type==tri_bi ? "tri+bi" : "bi+tri"); print (); printf (": "); fflush (stdout);

    // rAI-tre (ok, funny name, but you get the meaning): same as rAI,
    // but with trilinear terms. A product
    //
    // x1 x2 x3... xk
    //
    // is decomposed as
    //
    // (...((x1 x2 x3) x4 x5) x6 x7) ... ) x[k-1] xk)

    expression *aux = arglist_ [0];

    for (int i = 1; i < nargs_;) {

      if (i < nargs_ - 1 && ((type==tri_bi) || (i!=1))) { // this is the only point of departure between tri_bi and bi_tri

	// there are at least two remaining arguments: can create trilinear
	if (areSameVariables (aux, arglist_ [i])) {

	  if (areSameVariables (aux, arglist_ [i+1]))  // this is a term (x_i x_i x_i)

	    if (i == nargs_ - 2) aux =                    new exprPow (new exprClone (aux), new exprConst (3.));
	    else                 aux = p -> addAuxiliary (new exprPow (new exprClone (aux), new exprConst (3.)));

	  else {

	    aux = p -> addAuxiliary (new exprPow (new exprClone (aux), new exprConst (2.)));
	    //aux -> print (); printf (" (tri0) := "); fflush (stdout); aux -> Image () -> print (); printf ("; ");  fflush (stdout);
	    if (i == nargs_ - 2) aux =                    new exprMul (new exprClone (aux), new exprClone (arglist_ [i+1]));
	    else                 aux = p -> addAuxiliary (new exprMul (new exprClone (aux), new exprClone (arglist_ [i+1])));
	  }

	} else 
	  if (areSameVariables (aux, arglist_ [i+1])) {

	    printf ("Couenne, exprTrilinear: bad ordering of factors in product, aborting\n");
	    exit (-1);

	  } else {

	    aux = new exprTrilinear (new exprClone (aux), 
				     new exprClone (arglist_ [i]),
				     new exprClone (arglist_ [i+1]));

	    if (i != nargs_ - 2) 
	      aux = p -> addAuxiliary (aux);
	  }

	//aux -> print (); if (i != nargs_ - 2) {printf (" (tri) := "); aux -> Image () -> print (); printf ("; "); fflush (stdout);}

        i += 2; // covered two variables

      } else {

	if (areSameVariables (aux, arglist_ [i])) aux = new exprPow (new exprClone (aux), new exprConst (2.));
	else                                      aux = new exprMul (new exprClone (aux), new exprClone (arglist_ [i]));

	if (i==1) // must be bi+tri
	  aux = p -> addAuxiliary (aux); // necessary to introduce the auxiliary variable

	//aux -> print (); if (i==1) {printf (" (bi) := "); fflush (stdout); aux -> Image () -> print (); printf ("; "); fflush (stdout);}

	i++; // covered the last variable
      }
    }

    //printf ("\n");

    retExpr = aux; //(addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
  } break;

  // ----------------------------------------------------------------------------------------------

  case Couenne::rAI:
  default:

    //printf ("trying good ol' rAI on "); print (); printf ("\n");

    // rAI -- recursive Arithmetic Intervals (see Ryoo and Sahinidis,
    // JOGO 19 (2001):403-424):
    //
    // All multilinear expressions are decomposed as 
    //
    // (x_1 (x_2 (x_3... (x_{k-1} x_k))...)

    expression *aux = arglist_ [0];

    for (int i = 1; i < nargs_ - 1; i++)
      aux = (areSameVariables (aux, arglist_ [i])) ? 
	(p -> addAuxiliary (new exprPow (new exprClone (aux), new exprConst (2.)))) : 
	(p -> addAuxiliary (new exprMul (new exprClone (aux), new exprClone (arglist_ [i]))));

    if (areSameVariables (aux, arglist_ [nargs_ - 1])) 
      aux    = new exprPow (new exprClone (aux), new exprConst (2.));
    else aux = new exprMul (new exprClone (aux), new exprClone (arglist_ [nargs_ - 1]));

    retExpr = aux;
  }

  if (coeff != 1.) {

    retExpr = p -> addAuxiliary (retExpr);
    retExpr = new exprMul (new exprConst (coeff), new exprClone (retExpr));
  }

  return (addAux ? (p -> addAuxiliary (retExpr)) : new exprAux (retExpr, p -> domain ()));
}
