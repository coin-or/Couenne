/* $Id$
 *
 * Name:    writeLP.cpp
 * Author:  Pietro Belotti
 * Purpose: save a problem in LP format (MIQCQP only)
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <fstream>
#include <iomanip> // to use the setprecision manipulator

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"

#include "CouenneExprGroup.hpp"
#include "CouenneExprQuad.hpp"

using namespace Couenne;

// print linear or quadratic expression (only used here for LP file
// saving). Returns (unprintable) constant terms

double printLPquad (std::ofstream &f, const expression *ex, double mult) {

  expression *mutex = const_cast <expression *> (ex);

  exprGroup *e = dynamic_cast <exprGroup *> (mutex);

  double sign = 1., retval = 0;

  if (!e) {
    if (mutex -> code () == COU_EXPROPP) {
      ex = ex -> Argument () -> Original ();
      mutex = const_cast <expression *> (ex);
      sign = -1.;
      e = dynamic_cast <exprGroup *> (mutex);
    }
  }

  int wrapCount = 0;

  if (e) {

    // print linear part

    if (e -> getc0 () != 0)
      retval -= sign * e -> getc0 ();

    for (std::vector <std::pair <exprVar *, CouNumber> >::iterator i = e -> lcoeff () . begin (); 
	 i != e -> lcoeff (). end (); ++i) {

      CouNumber c = sign * i -> second;

      if ((c > 0) &&
	  ((i != e -> lcoeff (). begin ()) || (e -> getc0 () != 0)))
	f << '+';

      f << c << ' ';

      i -> first -> print (f);

      if (!(++wrapCount % 10))
	f << std::endl;
      else f << ' ';
    }
  } else {

    // assume the expression is a constant, take its value and print
    // it

    retval -= ex -> Value ();
  }

  // print quadratic part, if any
  exprQuad *q = dynamic_cast <exprQuad *> (mutex);

  // todo: an exprGroup has nonlinear components in argList_

  if (q) {

    std::vector <std::pair <exprVar *, exprQuad::sparseQcol> > &Q = q -> getQ ();

    f << " + [ ";

    for (std::vector <std::pair <exprVar *, exprQuad::sparseQcol> >::iterator i = Q. begin ();
	 i != Q . end (); ++i) {

      for (std::vector <std::pair <exprVar *, CouNumber> >::iterator j = i -> second. begin ();
	   j != i -> second . end (); ++j) {

	CouNumber c = mult * sign * j -> second;

	if (c > 0)
	  f << " +";

	f << c << " ";

	if (i -> first -> Index () == 
	    j -> first -> Index ()) {
	  i -> first -> print (f);
	  f << "^2";
	} else {
	  i -> first -> print (f); f << " * ";
	  j -> first -> print (f);
	}

	if (!(++wrapCount % 10))
	  f << std::endl;
	else f << ' ';
      }
    }

    f << " ] ";
    if (mult != 1) 
      f << "/ " << mult;
  }

  if ((mutex -> nArgs () == 1) && 
      (mutex -> ArgList () [0] -> Type () == CONST)) {

    retval -= sign * mutex -> ArgList () [0] -> Value ();

  } else {

    f << " + [ ";

    expression *arg = NULL, *arg2 = NULL;
    double signI;

    for (int i=0; i<mutex -> nArgs ();) {

      if (arg2) {

	arg = arg2;
	arg2 = NULL;
	signI = -signI;

      } else {

	signI = sign;

	arg = mutex -> ArgList () [i];

	if (arg -> code () == COU_EXPROPP) {

	  arg = arg -> Argument ();
	  signI = - sign;
	}

	if (arg -> code () == COU_EXPRSUB) {

	  arg2 = arg -> ArgList () [1];
	  arg  = arg -> ArgList () [0];
	}
      }

      if (arg -> code () == COU_EXPROPP) {

	arg = arg -> Argument ();
	signI = - signI;
      }

      if        (arg -> Type () == CONST)          retval -= signI * arg -> Value (); // constant

      else  if ((arg -> code () == COU_EXPRPOW) &&
		(arg -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> Type () == CONST)) {                       // xi^2

	f << ((signI > 0) ? " + " : " - ");
	if (mult != 1) f << mult << " ";
	arg -> ArgList () [0] -> print (f);
	f << "^2";
      }

      else  if ((arg -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [0] -> Type () == CONST) &&
		(arg -> ArgList () [1] -> code () == COU_EXPRPOW) &&
		(arg -> ArgList () [1] -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> ArgList () [1] -> Type () == CONST)) {     // alpha xi^2

	double c = mult * signI * arg -> ArgList () [0] -> Value ();
	f << ((c > 0) ? '+' : ' ') << c << " "; 
	arg -> ArgList () [1] -> ArgList () [0] -> print (f);
	f << "^2";
      }

      else  if ((arg -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> Type () == VAR) &&
		(arg -> ArgList () [0] -> Index () !=
		 arg -> ArgList () [1] -> Index ())) {                               // x1 * x2

	f << ((signI > 0) ? " + " : " - "); 
	if (mult != 1) f << mult << " ";
	arg -> ArgList () [0] -> print (f); f << " * ";
	arg -> ArgList () [1] -> print (f);
      }

      else  if ((arg -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> Type () == VAR) &&
		(arg -> ArgList () [0] -> Index () ==
		 arg -> ArgList () [1] -> Index ())) {                               // x1 * x1

	f << ((signI > 0) ? " + " : " - ");
	if (mult != 1) f << mult << " ";
	arg -> ArgList () [1] -> print (f);
	f << "^2 ";
      }

      else  if ((arg -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [0] -> Type () == CONST) &&
		(arg -> ArgList () [1] -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [1] -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> ArgList () [1] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> ArgList () [0] -> Index () !=
		 arg -> ArgList () [1] -> ArgList () [1] -> Index ())) {     // alpha x1 * x2

	double c = mult * signI * arg -> ArgList () [0] -> Value ();
	f << ((c > 0) ? '+' : ' ') << c << " "; 
	arg -> ArgList () [1] -> ArgList () [0] -> print (f); f << " * ";
	arg -> ArgList () [1] -> ArgList () [1] -> print (f);
      }

      else  if ((arg -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [0] -> Type () == CONST) &&
		(arg -> ArgList () [1] -> code () == COU_EXPRMUL) &&
		(arg -> ArgList () [1] -> ArgList () [0] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> ArgList () [1] -> Type () == VAR) &&
		(arg -> ArgList () [1] -> ArgList () [0] -> Index () !=
		 arg -> ArgList () [1] -> ArgList () [1] -> Index ())) {     // alpha x1 * x1

	double c = mult * signI * arg -> ArgList () [0] -> Value ();
	f << ((c > 0) ? '+' : ' ') << c << " "; 
	arg -> ArgList () [1] -> ArgList () [0] -> print (f);
	f << "^2 ";

      } else {

	printf ("Can't recognize expression (code: %d), exiting: ", arg -> code ());
	arg -> print ();
	printf ("\nExpression: ");
	ex -> print ();
	printf ("\n");
	exit (-1);
      }

      if (!(++wrapCount % 10))
	f << std::endl;
      else f << ' ';	

      if (!arg2) ++i;
      else if (arg -> code () == COU_EXPROPP)
	signI = - signI;
    }

    f << " ] ";
    if (mult != 1) 
      f << "/ " << mult;
  }

  return retval;
}


// store problem in a .mod file (AMPL)

void CouenneProblem::writeLP (const std::string &fname) { /// name of the mod file

  // only write LP problems from original. 

  for (int i=0; i < nVars (); i++)
    if (variables_ [i] -> Type () == AUX) {
    printf ("Auxiliary variables not supported in writeLP yet, bailing out\n");
    return;
  }

  if (objectives_ [0] -> Body () -> Linearity () > QUADRATIC) {
    printf ("Objective is nonlinear and not quadratic, bailing out\n");
    return;
  }

  for (int i=0; i < nCons (); i++) {
    if (constraints_ [i] -> Body () -> Linearity () > QUADRATIC) {
      printf ("Constraint %d is nonlinear and not quadratic, bailing out\n", i);
      return;
    }
  }

  std::ofstream f (fname.c_str ());

  f << std::setprecision (15);

  // original variables, integer and non //////////////////////////////////////////////

  f << "\\ Problem name (saved using Couenne): " << problemName_ << std::endl << std::endl;

  // objective function /////////////////////////////////////////////////////////////

  //if (objectives_ [0] -> Sense () == MINIMIZE) 
  f << "minimize obj: ";

  double objConst = printLPquad (f, objectives_ [0] -> Body (), 2);
  if (objConst != 0.) f << ((objConst > 0) ? " + " : " ") << objConst;
  f << std::endl << std::endl << "Subject To" << std::endl << std::endl;

  // // defined (aux) variables, with formula ///////////////////////////////////////////

  // if (aux) {

  //   f << std::endl << "# aux. variables defined" << std::endl << std::endl;

  //   for (int i=0; i < nVars (); i++)

  //     if ((variables_ [i] -> Type () == AUX) && 
  // 	  (variables_ [i] -> Multiplicity () > 0)) {

  // 	f << "aux_" << i << ": "; variables_ [i] -> print (f, false);
  // 	f << " = ";  

  // 	variables_ [i] -> Image () -> print (f, false);
  // 	f << std::endl;
  //     }
  // }

  // write constraints //////////////////////////////////////////////////////////////

  // f << std::endl << "# constraints" << std::endl << std::endl;

  // if (!aux) // print h_i(x,y) <= ub, >= lb
  //   for (std::vector <exprVar *>::iterator i = variables_.begin ();
  // 	 i != variables_.end ();
  // 	 ++i) 

  //     if (((*i) -> Type () == AUX) && 
  // 	  ((*i) -> Multiplicity () > 0)) {
	
  // 	CouNumber bound;

  // 	if ((bound = (*i) -> lb ()) > - COUENNE_INFINITY) {
  // 	  f << "conAuxLb" << (*i) -> Index () << ": ";
  // 	  (*i) -> print (f, true);
  // 	  f << ">= " << bound << std::endl;
  // 	}

  // 	if ((bound = (*i) -> ub ()) <   COUENNE_INFINITY) {
  // 	  f << "conAuxUb" << (*i) -> Index () << ": ";
  // 	  (*i) -> print (f, true);
  // 	  f << "<= " << bound << std::endl;
  // 	}
  //     }

  for (int i=0; i < nCons (); i++) {

    // get numerical value of lower, upper bound
    CouNumber lb = (constraints_ [i] -> Lb () -> Value ()),
              ub = (constraints_ [i] -> Ub () -> Value ());

    f << "con_" << i << ": ";
    double extraTerm = printLPquad (f, constraints_ [i] -> Body (), 1);

    lb += extraTerm;
    ub += extraTerm;

    if (lb > - COUENNE_INFINITY) {
      f << ' ';
      if (fabs (ub-lb) > COUENNE_EPS) 
	f << '>';
      f << "= " << lb << std::endl;
    }
    else f << " <= " << ub << std::endl;

    // if range constraint, print it once again

    if ((   lb > - COUENNE_INFINITY + 1) 
	&& (ub <   COUENNE_INFINITY - 1)
	&& (fabs (ub-lb) > COUENNE_EPS)) {

      f << "con_" << i << "_rng: ";
      printLPquad (f, constraints_ [i] -> Body (), 1);
      f << " <= " << ub << std::endl;
    }
  }

  f << "Bounds" << std::endl << std::endl;

  for (int i=0; i < nVars (); i++) {

    if ((Lb (i) == 0) && (Ub (i) >= COUENNE_INFINITY/2))
      continue;

    if (Lb (i) != 0) f << Lb (i) << " <= ";
    variables_ [i] -> print (f);
    if (Ub (i) < + COUENNE_INFINITY/2) f << " <= " << Ub (i);
    f << std::endl;
  }

  f << "Generals" << std::endl << std::endl;

  int wrapcount = 0;

  for (int i=0; i < nVars (); i++)

    if (variables_ [i] -> isInteger ()) {
      variables_ [i] -> print (f);
      if (!(++wrapcount % 10))
	f << std::endl;
      else 
	f << " ";
    }

  f << std::endl << std::endl << "End" << std::endl;

  f.close ();
}
