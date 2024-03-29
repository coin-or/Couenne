/*
 *
 * Name:    splitAux.cpp
 * Author:  Pietro Belotti
 * Purpose: extract auxiliary variable from implicit equality constraint
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprQuad.hpp"

#include "CouenneProblem.hpp"

#include "CoinHelperFunctions.hpp"

#include "CouenneExprAux.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneLQelems.hpp"

using namespace Couenne;

/// given an element of a sum, check if it is a variable (possibly
/// with a coefficient) and return its index (and the coefficient) if
/// it has not been spotted as an auxiliary
void elementBreak (expression *, int &, CouNumber &);


/// split a constraint aw + f(x) >/</= c into w's index (returned)
/// and rest = (-f(x) + c)/a

int CouenneProblem::splitAux (CouNumber rhs, expression *body, expression *&rest,
			      bool *wentAux, enum expression::auxSign &sign) {

  int auxInd = -1,          // index of the auxiliary to be extracted
    code = body -> code (); // type of expression

  expression **alist = body -> ArgList ();

  if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE))
    {printf ("  Splitting expression: "); body -> print (); printf ("\n");}

  switch (code) { // constraint h(x) = 0 may be in different forms:
		  // subtraction, sum, linear group

    /////////////////////////////////////////////////////////////////////////////

  case COU_EXPRSUB: {

    // simplest case, we have w-f(x)=rhs or f(x)-w=rhs or f(x)-g(x)=rhs

    int pos = 0;
    CouNumber coeff = 1;

    auxInd = (*alist) -> Index ();

    if (auxInd < 0) // first argument is not a variable
      elementBreak (*alist, auxInd, coeff); // check first element

    if ((auxInd < 0) ||       // first argument is not a variable
	wentAux [auxInd] ||   // it was, and it already became an auxiliary
	(alist [1] -> dependsOn (auxInd, TAG_AND_RECURSIVE) >= 1)) {
                              // it is not an auxiliary, but the remaining
                              // expression depends on it

      auxInd = (alist [1]) -> Index ();

      if (auxInd < 0)
	elementBreak (alist [1], auxInd, coeff); // check second element
      else coeff = 1;

      if ((auxInd < 0) ||       // no index found
	  wentAux [auxInd] ||   // or, found but invalid
	  (alist [0] -> dependsOn (auxInd, TAG_AND_RECURSIVE) >= 1))
	                        // or, variable depends upon itself
	return -1;

      pos = 1;
    }

    //////////////////////

    // what remains is the "independent" expression
    expression *clone = alist [1 - pos] -> clone (&domain_);

    expression *auxdef;

    if      (coeff ==  1.) auxdef =                                     clone;  // if coefficient is  1, do nothing
    else if (coeff == -1.) auxdef = new exprOpp                        (clone); // if coefficient is -1, return -f(x)
    else                   auxdef = new exprMul (new exprConst (coeff), clone); // else multiply it by coefficient

    rest = (rhs == 0.) ?                                             // no extra constant?
      (auxdef) :                                                     // just put other argument
      (new exprSum (auxdef, new exprConst ((pos==1) ? -rhs : rhs))); // otherwise sum it with \pm rhs

    if (pos==1) {
      if      (sign == expression::AUX_GEQ)  sign = expression::AUX_LEQ;
      else if (sign == expression::AUX_LEQ)  sign = expression::AUX_GEQ;
    }

  } break;

  ////////////////////////////////////////////////////////////////////////////

  case COU_EXPRQUAD:
  case COU_EXPRGROUP:
  case COU_EXPRSUM: {

    // an expr{Sum,Group,Quad}. Decompose the expression and
    // re-assemble (to look for right variable)

    // data structure to be used below if there is a linear term.
    // which specifies position within arrays (negative for linear
    // part of exprGroup, positive for all elements of exprSum)

    int maxindex = -1, nlin = 0, which = 1;
    CouNumber c0 = 0., auxcoe = 1;
    bool which_was_set = false;

    // check indices of linear/quadratic part /////////////////////////////////

    if (code != COU_EXPRSUM) {

      exprGroup *egBody = dynamic_cast <exprGroup *> (body -> isaCopy () ?
						      body -> Copy () :
						      body);
      exprGroup::lincoeff &lcoe = egBody -> lcoeff ();

      // import exprGroup linear structure

      c0 = egBody -> getc0 ();

      for (int i=0, n = lcoe.size (); n--; i++, nlin++) {

	int j = lcoe [i]. first -> Index ();

	//lincoe [i] = lcoe [i]. second;

	// prefer non-integer. If integer, only take it if none chosen yet
	if ((j > maxindex) &&
	    !(wentAux [j]) &&
	    (fabs (lcoe [i]. second) > COUENNE_EPS) &&
	    (!(lcoe [i].first -> isInteger ()) || (which==1))) {

	  // fake cut in linind and check dependence. Only mark if
	  // dependsOn() gives 0

	  exprVar *saveVar = lcoe [i].first;
	  //lcoe [i].first = new exprVar (nVars ()); // better use index -1
	  lcoe [i].first = new exprVar (-1);

	  if (body -> dependsOn (j, TAG_AND_RECURSIVE) == 0) {

	    if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	      printf ("body does not depend on x%d: ", j);
	      body -> print ();
	      printf ("\n");
	    }
	    // mark which with negative number
	    which    = - nlin - 1;
	    auxcoe   = lcoe [i]. second;
	    maxindex = j;
	  }

	  delete lcoe [i].first;
	  lcoe [i].first = saveVar;
	}
      }
    }

    if (which != 1)
      which_was_set = true;
    else which = -1;

    // check indices of elements of (nonlinear) sum /////////////////////////////////

    for (int i = body -> nArgs (); i--;) {

      CouNumber coeff = 1;
      int index = alist [i] -> Index ();

      // if index < 0, this is an expression and not a single variable
      if (index < 0)
	elementBreak (alist [i], index, coeff);

      if ((index > maxindex) &&
	  !(wentAux [index]) &&
	  (fabs (coeff) > COUENNE_EPS)) {

	// fake a cut in the arglist and check

	expression *cut = alist [i];
	alist [i] = new exprConst (0.);

	// not enough... check now linear (and quadratic!) terms

	if (body -> dependsOn (index, TAG_AND_RECURSIVE) == 0) {

	  maxindex = index;
	  which    = i;
	  auxcoe   = coeff;
	}

	delete alist [i];
	alist [i] = cut;
      }
    }

    if (!which_was_set && (which == -1)) // which has not been set
      return -1;

    ///////////////////////////////////////////////////////////////////

    if (maxindex < 0) break; // no substitute found ==> no hidden auxiliary

    // create a new exprGroup, exprQuad, or exprSum with all elements but the
    // extracted auxiliary

    // start with exprSum
    int nargs = body -> nArgs ();
    expression **newarglist;

    if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
      printf (" [[ind %d, coe %.1g, wh %d, nargs %d, nlin %d]] ",
	      maxindex, auxcoe, which, nargs, nlin);
    }

    if (nargs > 0) { // there is an element in the nonlinear sum to be drawn

      int j, mid = (which < 0) ? nargs : which;

      newarglist = new expression * [nargs + 1];

      for (j=0; j<mid;   j++) newarglist [j]   = alist [j];// -> clone (&domain_);
      for (j++; j<nargs; j++) newarglist [j-1] = alist [j];// -> clone (&domain_);

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient

    } else { // no nonlinear arguments, or the only one is the new aux

      nargs++; // !!!!!!!!!!!!!!!!!!!!!!!!!
      newarglist  = new expression *;
      *newarglist = new exprConst (0.);
    }

    // form rhs linear part ////////////////////////////////////////////////////

    int       *linind2 = NULL;
    CouNumber *lincoe2 = NULL;

    // in case this was (and will be) an exprQuad
    int *qindI = NULL,
        *qindJ = NULL;
    CouNumber *qcoe = NULL;

    if (nlin > 0) { // there is an element in the linear sum to be drawn

      exprGroup *egBody = dynamic_cast <exprGroup *> (body -> isaCopy () ?
						      body -> Copy () :
						      body);

      exprGroup::lincoeff &lcoe = egBody -> lcoeff ();

      int mid = (which >= 0) ? nlin : - which - 1;

      linind2 = new int       [nlin + 1];
      lincoe2 = new CouNumber [nlin + 1];

      int j;

      //if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
      //for (j=0; j<mid;  j++) printf ("{%g x%d} ", lincoe [j], linind [j]);
      //for (j++; j<nlin; j++) printf ("{%g x%d} ", lincoe [j], linind [j]);
      //}

      CouNumber divider = -1. / auxcoe;

      for (j=0; j<mid;  j++){linind2[j]  =lcoe[j].first->Index();lincoe2[j]  =divider*lcoe[j].second;}
      for (j++; j<nlin; j++){linind2[j-1]=lcoe[j].first->Index();lincoe2[j-1]=divider*lcoe[j].second;}

      linind2 [j-1] = -1; // terminate list of indices

      if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
	for (j=0; j<mid;  j++) printf ("<%g x%d> ", lincoe2 [j],   linind2 [j]);
	for (j++; j<nlin; j++) printf ("<%g x%d> ", lincoe2 [j-1], linind2 [j-1]);
      }

      if (code == COU_EXPRQUAD) { // copy quadratic elements

	exprQuad *eq = dynamic_cast <exprQuad *> (body -> isaCopy () ?
						  body -> Copy () :
						  body);

	int nqt = eq -> getnQTerms (), j=0;

	qindI = new int [nqt];
	qindJ = new int [nqt];
	qcoe  = new CouNumber [nqt];

	exprQuad::sparseQ &M = eq -> getQ ();

	for (exprQuad::sparseQ::iterator row = M.begin ();
	     row != M.end (); ++row) {

	  int xind = row -> first -> Index ();

	  for (exprQuad::sparseQcol::iterator col = row -> second.begin ();
	       col != row -> second.end (); ++col, ++j) {

	    qindI [j] = xind;
	    qindJ [j] = col -> first -> Index ();
	    qcoe  [j] = divider * col -> second;
	  }
	}
      }

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient
    }

    // the extracted index is one term of...
    if (which >= 0) --nargs; // ...the nonlinear sum
    else            --nlin;  // ...the linear part

    jnlst_ -> Printf (Ipopt::J_ALL, J_REFORMULATE, "\n::: auxcoe %g, rhs %g, lin %d, nl %d\n", auxcoe, rhs, nlin, nargs);

    // all is ready to take the independent stuff to the other side of
    // the inequality.

    if ((code == COU_EXPRQUAD) ||
	((code == COU_EXPRGROUP) && (nlin > 0))) {

      // an exprGroup with at least one linear term left
      //
      // build new vectors for index and coeff. Two cases:
      //
      // 1)  f(x) + c0 -  w = rhs   =====>   w =       f(x) + c0 - rhs
      // 2)  f(x) + c0 + aw = rhs   =====>   w = -1/a (f(x) + c0 - rhs), a != -1

      if (fabs (auxcoe + 1) < COUENNE_EPS) {

	//std::vector <std::pair <exprVar *, CouNumber> >
	exprGroup::lincoeff lcoeff;
	indcoe2vector (linind2, lincoe2, lcoeff);

	if (code == COU_EXPRGROUP)
	  rest = new exprGroup (c0-rhs, lcoeff, newarglist, nargs);

	else {
	  std::vector <quadElem> qcoeff;
	  indcoe2vector (qindI, qindJ, qcoe, qcoeff);
	  rest = new exprQuad  (c0-rhs, lcoeff, qcoeff, newarglist, nargs);
	}
      }
      else {

	expression **mullist = new expression * [1];

	// only nl term is constant
	if ((nargs <= 1) && ((*newarglist) -> Linearity () <= CONSTANT)) {
	  *mullist = new exprConst ((*newarglist) -> Value ());
	  //delete *newarglist;
	  delete [] newarglist;
	}
	else if ((nargs <= 1) &&
		 ((*newarglist) -> code () == COU_EXPROPP) &&
		 (fabs (auxcoe - 1) < COUENNE_EPS))
	  //*mullist = new exprSum (newarglist, nargs);
	  *mullist = (*newarglist) -> Argument ();//, nargs);
	else { // multiple nonlinear terms, multiply them by -1/auxcoe

	  if (nargs <= 1)
	    *mullist = new exprMul (new exprConst (-1. / auxcoe),
				    *newarglist); // !!! one more leak?
	  else
	    *mullist = new exprMul (new exprConst (-1. / auxcoe),
				    new exprSum (newarglist, nargs));
	}

	std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
	indcoe2vector (linind2, lincoe2, lcoeff);

	// final outcome: -1/a (f(x) + c0 - rhs)
	if (code == COU_EXPRGROUP)
	  rest = new exprGroup ((rhs - c0) / auxcoe, lcoeff,        mullist,1);

	else {
	  std::vector <quadElem> qcoeff;
	  indcoe2vector (qindI, qindJ, qcoe, qcoeff);
	  rest = new exprQuad  ((rhs - c0) / auxcoe, lcoeff, qcoeff, mullist,1);
	}
      }
    }
    else { // simple exprSum

      if (fabs (rhs) > COUENNE_EPS) { // have to add constant to exprSum

	if ((nargs == 1) && ((*newarglist) -> Type () == CONST)) {

	  CouNumber val = (*newarglist) -> Value () - rhs;
	  //delete *newarglist;
	  *newarglist = new exprConst (val);
	}
	else newarglist [nargs++] = new exprConst (-rhs);
      }

      // now exprSum is complete with -rhs. Send it to right hand side

      expression *auxDef;

      if (nargs==1) {
	auxDef = *newarglist;
	delete [] newarglist;
      } else auxDef = new exprSum (newarglist, nargs);

      if (auxcoe == -1.)
	rest = auxDef;
      else if ((auxcoe == 1.) &&
	       (auxDef -> code () == COU_EXPROPP))
	rest = auxDef -> Argument ();
      else
	if (auxDef -> Type () == CONST)

	  rest = new exprConst (- auxDef -> Value () / auxcoe);
	else {

	  // TODO: check if auxdef is an exprMul (k*f(x)) and
	  // -1/auxcoe simplifies

	  if (auxDef -> code () == COU_EXPROPP) {

	    auxDef = auxDef -> Argument ();
	    auxcoe = -auxcoe;
	  }

	  // TODO: if it's a c*(ax+b) return new exprGroup with (c*a, c*b)

	  rest = new exprMul (new exprConst (-1./auxcoe),
			      ((auxDef -> Type () == AUX) ||
			       (auxDef -> Type () == VAR)) ? new exprClone (auxDef) : auxDef);
	}
    }

    if (linind2) {
      delete [] linind2;
      delete [] lincoe2;
    }

    if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
      printf ("  Generated "); rest -> print (); printf ("\n");
    }

    auxInd = maxindex;

    if (auxcoe < 0) {
      if      (sign == expression::AUX_GEQ)  sign = expression::AUX_LEQ;
      else if (sign == expression::AUX_LEQ)  sign = expression::AUX_GEQ;
    }

  } break;

  default: break;
  }

  if ((auxInd < 0) || (wentAux [auxInd]))
    return -1;

  // we have a variable, meaning the constraint body is of the form
  // w +/- f(x) or f(x) +/- w

  // standardize remaining of the expression

  if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
    printf   ("  Standardize rest (2nd level). Rest: "); fflush (stdout); rest -> print ();
    printf ("\n                                Body: "); fflush (stdout); body -> print ();
  }

  // second argument is set to false so as to instruct standardize not
  // to create a new auxiliary variable (we are creating it ourselves,
  // it's aux)

  exprAux *aux = rest -> standardize (this, false);

  if (aux && jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
      printf ("  Aux: "); aux -> print ();
      printf (" := ");    aux -> Image () -> print ();
  }

  if (aux) {
    //delete rest;
    //rest = aux -> Image () -> clone (&domain_);
    rest = aux -> Image (); // -> clone (&domain_);
    aux -> Image (NULL);
    //delete aux; // segfaults on simplified expressions? (see test/simple.nl)
  }

  if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE))
    {printf (" ==> (%d)", auxInd); rest -> print (); printf ("\n");}

  return auxInd;
}
