/* $Id$
 *
 * Name:    conv-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: utility methods to convexify multiplications
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
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

  exprOp::standardize (p);

  if (nargs_ == 1)  // TODO: what happens really?
    return NULL;
  /* {
     exprAux *aux = arglist_ [0];
     arglist_ [0] = NULL;
     return aux;
     } */

  // check if it is a product of binary variables
  bool isBinProd = true;
  for (int i=nargs_; i--;) {

    expression *arg = arglist_ [i];
    if (arg -> isInteger ()) {

      CouNumber lb, ub;
      arg -> getBounds (lb, ub);
      if ((fabs (lb)      > 0.) ||
	  (fabs (ub - 1.) > 0.)) { // if not all conditions hold, 
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

  enum Couenne::TrilinDecompType type = p -> getTrilinDecompType ();

  if (nargs_ <= 2)
    type = Couenne::rAI;

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

    return (addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
  }

    // ----------------------------------------------------------------------------------------------

  case Couenne::tri_bi:
  case Couenne::bi_tri: { // the two cases are very similar

    printf ("trying %s on ", type==tri_bi ? "tri+bi" : "bi+tri"); print (); printf (": "); fflush (stdout);

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
	    aux -> print (); printf (" (tri0) := "); fflush (stdout); aux -> Image () -> print (); printf ("; ");  fflush (stdout);
	    if (i == nargs_ - 2) aux =                    new exprMul (new exprClone (aux), new exprClone (arglist_ [i+1]));
	    else                 aux = p -> addAuxiliary (new exprMul (new exprClone (aux), new exprClone (arglist_ [i+1])));
	  }

	} else 
	  if (areSameVariables (aux, arglist_ [i+1])) {

	    printf ("Couenne, exprTrilinear: bad ordering of factors in product, aborting\n");
	    exit (-1);

	  } else
	    if (i == nargs_ - 2) 
	      aux = new exprTrilinear (new exprClone (aux), 
				       new exprClone (arglist_ [i]),
				       new exprClone (arglist_ [i+1]));
	    else
	      aux = p -> addAuxiliary (new exprTrilinear (new exprClone (aux), 
							  new exprClone (arglist_ [i]),
							  new exprClone (arglist_ [i+1])));

	aux -> print (); if (i != nargs_ - 2) {printf (" (tri) := "); aux -> Image () -> print (); printf ("; "); fflush (stdout);}

        i += 2; // covered two variables

      } else {

	if (areSameVariables (aux, arglist_ [i])) aux = new exprPow (new exprClone (aux), new exprConst (2.));
	else                                      aux = new exprMul (new exprClone (aux), new exprClone (arglist_ [i]));

	if (i==1) // must be bi+tri
	  aux = p -> addAuxiliary (aux); // necessary to introduce the auxiliary variable

	aux -> print (); if (i==1) {printf (" (bi) := "); fflush (stdout); aux -> Image () -> print (); printf ("; "); fflush (stdout);}

	i++; // covered the last variable
      }
    }

    printf ("\n");

    return (addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
  }

    // ----------------------------------------------------------------------------------------------

    // ----------------------------------------------------------------------------------------------

  case Couenne::rAI:
  default:

    //printf ("trying good ol' rAI on "); print (); printf ("\n");

    // rAI (recursive Arithmetic Intervals (see Ryoo and Sahinidis,
    // JOGO 19 (2001):403-424):
    //
    // All multilinear expressions are decomposed as 
    //
    // (x_1 (x_2 (x_3... (x_{k-1} x_k))...)

    //expression *aux = arglist_ [0]; // why not this one?
    expression *aux = new exprClone (arglist_ [0]);

    for (int i = 1; i < nargs_ - 1; i++)
      aux = (areSameVariables (aux, arglist_ [i])) ? 
	(p -> addAuxiliary (new exprPow (new exprClone (aux), new exprConst (2.)))) : 
	(p -> addAuxiliary (new exprMul (new exprClone (aux), new exprClone (arglist_ [i]))));

    if (areSameVariables (aux, arglist_ [nargs_ - 1])) 
      aux    = new exprPow (new exprClone (aux), new exprConst (2.));
    else aux = new exprMul (new exprClone (aux), new exprClone (arglist_ [nargs_ - 1]));

    return (addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
  }
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (expression *&lb, expression *&ub) {

  int i;

  if ((arglist_ [i=0] -> Type () == CONST) ||
      (arglist_ [i=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [i] -> Value ();

    if (!i && (arglist_ [1] -> Type () == CONST)) { 

      // !i means i==0, or the first is constant. If you are here,
      // both are constant, which should not happen...

      CouNumber prod = c * arglist_ [1] -> Value ();

      lb = new exprConst (prod);
      ub = new exprConst (prod);

      return;
    }
    else {

      // expression is of the type c*x

      expression *lbi, *ubi;
      arglist_ [1-i] -> getBounds (lbi, ubi);

      if (c >= 0) {
	lb = new exprMul (new exprConst (c), lbi);
	ub = new exprMul (new exprConst (c), ubi);
      } else {
	lb = new exprMul (new exprConst (c), ubi);
	ub = new exprMul (new exprConst (c), lbi);
      }
    }
  }
  else {

    // expression is of the type x*y

    expression **almin = new expression * [4];
    expression **almax = new expression * [4];

    arglist_ [0] -> getBounds (almin [0], almin [1]);
    arglist_ [1] -> getBounds (almin [2], almin [3]);

    almax [0] = new exprClone (almin [0]);
    almax [1] = new exprClone (almin [1]);
    almax [2] = new exprClone (almin [2]);
    almax [3] = new exprClone (almin [3]);

    lb = new exprLBMul (almin, 4);
    ub = new exprUBMul (almax, 4);
  }
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lb1, ub1, lb2, ub2;

  arglist_ [0] -> getBounds (lb1, ub1);
  arglist_ [1] -> getBounds (lb2, ub2);

  if (ub1 < 0) { // use lb1, dominant
    if      (ub2 < 0) {lb = safeProd(ub1,ub2); ub = safeProd(lb1,lb2);}
    else if (lb2 > 0) {lb = safeProd(lb1,ub2); ub = safeProd(ub1,lb2);}
    else              {lb = safeProd(lb1,ub2); ub = safeProd(lb1,lb2);}
  } else if (lb1 > 0) { // use ub1, dominant
    if      (ub2 < 0) {lb = safeProd(ub1,lb2); ub = safeProd(lb1,ub2);}
    else if (lb2 > 0) {lb = safeProd(lb1,lb2); ub = safeProd(ub1,ub2);}
    else              {lb = safeProd(ub1,lb2); ub = safeProd(ub1,ub2);}
  } else { // there is a zero to consider
    if      (ub2 < 0) {lb = safeProd(ub1,lb2); ub = safeProd(lb1,lb2);}
    else if (lb2 > 0) {lb = safeProd(lb1,ub2); ub = safeProd(ub1,ub2);}
    else              {lb = CoinMin (safeProd(lb1,ub2), safeProd(lb2,ub1));
                       ub = CoinMax (safeProd(lb1,lb2), safeProd(ub1,ub2));}
  }
}
