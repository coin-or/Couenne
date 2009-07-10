/* $Id$ */
/*
 * Name:    flattenMul.cpp
 * Author:  Pietro Belotti
 * Purpose: flatten multiplication expression tree into monomial
 *          c*\Prod_{k\in K} x_{i_k}^{p_k}
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "exprClone.hpp"

//#define DEBUG

/// re-organizes multiplication and stores indices (and exponents) of
/// its variables
void CouenneProblem::flattenMul (expression*& mul, CouNumber &coe, 
				 std::map <int, CouNumber> &indices) {

#ifdef DEBUG
  printf ("flatten %d ---> ", mul -> code ()); mul -> print ();
  printf ("\n");
#endif

  if (mul -> code () != COU_EXPRMUL) {

    exprAux *aux = mul -> standardize (this);

    int ind = (aux) ? aux -> Index () : mul -> Index ();
    if (aux) {
      if (mul->code() == COU_EXPRGROUP ||
          mul->code() == COU_EXPRSUM || 
          mul->code() == COU_EXPRSUB || 
          mul->code() == COU_EXPROPP || 
          mul->code() == COU_EXPRQUAD)
        delete mul;
      mul = new exprClone(aux);
    }

    std::map <int, CouNumber>::iterator 
      where = indices.find (ind);

    if (where == indices.end ()) 
      indices.insert (std::pair <int, CouNumber> (ind, 1));
    else ++ (where -> second);

    return;
  }

  int nargs = mul -> nArgs ();
  expression **al = mul -> ArgList ();

  // for each factor (variable, function, or constant) of the product
  for (int i=0; i < nargs; i++) { 

    expression *arg = al [i];

#ifdef DEBUG
    printf ("  flatten arg %d ---> ", arg -> code ()); arg -> print ();
    printf ("\n");
#endif

    switch (arg -> code ()) {

    case COU_EXPRCONST: // change scalar multiplier

      coe *= arg -> Value ();
      break;

    case COU_EXPRMUL:  // apply recursively

      flattenMul (arg, coe, indices);
      break;

    case COU_EXPRVAR: { // insert index or increment 

      std::map <int, CouNumber>::iterator 
	where = indices.find (arg -> Index ());

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (arg -> Index (), 1));
      else ++ (where -> second);
    } break;

    case COU_EXPROPP: // equivalent to multiplying by -1

      coe = -coe;

      if (arg -> Argument () -> Type () == N_ARY) {
	flattenMul (*arg -> ArgPtr (), coe, indices);
	break;
      } else arg = arg -> Argument ();

    case COU_EXPRPOW: 

      if (arg -> code () == COU_EXPRPOW) { // re-check as it could come from above

      expression *base     = arg -> ArgList () [0],
	         *exponent = arg -> ArgList () [1];

      if (exponent -> Type () == CONST) { // could be of the form k x^2

	double expnum = exponent -> Value ();

	expression *aux = base -> standardize (this);

	if (!aux)
	  aux = base;

	std::map <int, CouNumber>::iterator 
	  where = indices.find (aux -> Index ());

	if (where == indices.end ())
	  indices.insert (std::pair <int, CouNumber> (aux -> Index (), expnum));
	else (where -> second += expnum);

	break;
      }  // otherwise, revert to default
    }

    case COU_EXPRSUM: // well, only if there is one element

      if ((arg -> code  () == COU_EXPRSUM) && // re-check as it could come from above
	  (arg -> nArgs () == 1)) {

	flattenMul (arg, coe, indices);
	break;

      } // otherwise, continue into default case

    default: { // for all other expression, add associated new auxiliary

      exprAux *aux = arg -> standardize (this);

      int ind = (aux) ? aux -> Index () : arg -> Index ();

      if (aux) {
        if (al[i]->code() == COU_EXPROPP) {
          *al[i]->ArgPtr() = NULL;
          delete al[i];
        } else
          assert(arg == al[i]);
        if (arg->code() == COU_EXPRGROUP ||
            arg->code() == COU_EXPRSUM || 
            arg->code() == COU_EXPRSUB || 
            arg->code() == COU_EXPROPP || 
            arg->code() == COU_EXPRQUAD)
        delete arg;
	al[i] = new exprClone(aux);
      }

      std::map <int, CouNumber>::iterator 
	where = indices.find (ind);

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (ind, 1));
      else ++ (where -> second);

    } break;
    }
  }
}
