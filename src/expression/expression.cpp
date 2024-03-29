/*
 *
 * Name:    expression.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <iostream>

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprOp.hpp"
#include "CouenneExprUnary.hpp"
#include "CouenneExprStore.hpp"

using namespace Couenne;

// empty differentiation method
expression *expression::differentiate (int)
{return new exprConst (0.);}


// Get lower and upper bound of a generic expression
void expression::getBounds (expression *&lb, expression *&ub) {

  lb = new exprConst (- COIN_DBL_MAX);
  ub = new exprConst (  COIN_DBL_MAX);
}


/// Get lower and upper bound of an expression (if any) -- real values
void expression::getBounds (CouNumber &lb, CouNumber &ub) {

  expression *le, *ue;
  getBounds (le, ue);

  lb = (*le) ();
  ub = (*ue) ();

  delete le;
  delete ue;
}


// generate one cut for a constant
void exprConst::generateCuts (expression *w, //const OsiSolverInterface &si,
			      OsiCuts &cs, const CouenneCutGenerator *cg,
			      t_chg_bounds *chg, int,
			      CouNumber, CouNumber) {
  if (cg -> isFirst ())
    cg -> createCut (cs, value_, 0, w -> Index (), 1.);
}


/// compare generic expression with other expression
int expression::compare (expression &e1) {

  int c0 =     code (),
               c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;

  // same code, check arguments

  if (c0 >= COU_EXPRUNARY) { // both are exprUnary's. COU_EXPRUNARY > COU_EXPROP, so if this is run then the next is not

    exprUnary *ne0 = dynamic_cast <exprUnary *> (const_cast <expression *> (this->Original()));
    exprUnary *ne1 = dynamic_cast <exprUnary *> (const_cast <expression *> (e1.Original()));

    return ne0 -> compare (*ne1);
  }

  if (c0 >= COU_EXPROP) { // both are exprOp's

    exprOp *ne0 = dynamic_cast <exprOp *> (const_cast <expression *> (this->Original()));
    exprOp *ne1 = dynamic_cast <exprOp *> (const_cast <expression *> (e1.Original()));

    return ne0 -> compare (*ne1);
  }

  // expressions are both variables or constants

  {
    int
      i0 =     Index (),
      i1 = e1. Index ();

    if (i0 < i1) return -1;
    if (i0 > i1) return  1;
    if (i0 >= 0) return  0; // same variables
  }

  // both are numbers
  {
    CouNumber
      v0 =     Value (),
      v1 = e1. Value ();

    if (v0 < v1) return -1;
    if (v0 > v1) return  1;
  }

  return  0;
}


/// Copy constructor
exprCopy::exprCopy (const exprCopy &e, Domain *d) {

  //expression *copy = getOriginal (e.copy_);

  //if ((copy -> Type () == VAR) ||
  //    (copy -> Type () == AUX))

    // if this is a variable, don't copy, just refer to the same
    // object (otherwise many useless copies are created)

    //  copy_ = copy;// -> clone (d);
  //else copy_ = copy -> clone (d);

  copy_ = e.copy_ -> Original () -> clone (d); // !!! REMOVE " -> Original ()"
  //else copy_ = e.copy_ -> clone (d);

  value_ = e.value_;
}


/// compare expressions (used in bsearch within CouenneProblem::standardize)
int expression::compare (exprCopy &c)
{return compare (const_cast <expression &> (*(c. Original ())));}


/// replace occurrence of a variable with another variable
void exprCopy::replace (exprVar *orig, exprVar *aux) {

  if (!aux)
    aux = orig;

  enum nodeType copyType = copy_ -> Type ();

  if ((copyType == VAR) ||
      (copyType == AUX)) {

    if (copy_ -> Index () == orig -> Index ()) {
      if (copy_ -> isaCopy ())  // !!! REMOVE
	delete copy_;           // !!! REMOVE
      copy_ = aux;
    }

  } else copy_ -> replace (orig, aux);

  //copy_ -> replace (orig, aux);

  /*if ((copy_ -> Index () == orig -> Index ()) && (copy_ -> Type  () != AUX)) {
    delete copy_;
    copy_ = new exprClone (aux);
    }*/
}


/// dependence on variable set: return cardinality of subset of the
/// set of indices in first argument which occur in expression.
int expression::dependsOn (int *ind, int n, enum dig_type type) {

  std::set <int>
    indlist (ind, ind + n),
    deplist;
  //intersectn;

  /*printf (":::::: indlist = {");
  for (std::set <int>::iterator i=indlist.begin (); i != indlist.end(); ++i)
    printf ("%d ", *i);
    printf ("} -- deplist = {");*/

  DepList (deplist, type);

  /*for (std::set <int>::iterator i=deplist.begin (); i != deplist.end(); ++i)
    printf ("%d ", *i);
    printf ("} -- intersection: \n");*/

  for (std::set <int>::iterator
	 i = deplist.begin (),
	 j = indlist.begin ();
       (i != deplist.end ()) &&
	 (j != indlist.end ());) {

    if (*i == *j) return 1;

    if (*i > *j) ++j;
    else         ++i;
  }

  return 0;
}


/// redirect variables to proper variable vector
void exprCopy::realign (const CouenneProblem *p) {

  if (((copy_ -> Type () == VAR) ||
       (copy_ -> Type () == AUX)) &&
      (copy_ -> Original () != p -> Var (copy_ -> Index ()))) {

    /*printf ("exprCopy::realign replaces %x with %x (",
	    copy_ -> Original (), p -> Var (copy_ -> Index ()));
    copy_ -> Original () -> print (); printf (" --> ");
    p -> Var (copy_ -> Index ()) -> print (); printf (")\n");*/

    expression *trash = copy_;

    copy_ = p -> Var (copy_ -> Index ());
    //delete trash; // Otherwise ticket 13
  } else copy_ -> realign (p);
}

/// printing method for copy expressions
void exprCopy::print (std::ostream &out, bool descend) const
{copy_ -> Original () -> print (out, descend);}
//{out << "["; copy_ -> print (out, descend); out << "]<" << copy_ << ">"; } // Must go


/// printing method
void exprClone::print (std::ostream &out, bool descend) const
{copy_ -> Original () -> print (out, descend);}
//{out << "{"; copy_ -> print (out, descend); out << "}<" << copy_ << ">"; } // Must go


/// printing method
void exprStore::print (std::ostream &out, bool descend) const
{copy_ -> Original () -> print (out, descend);}
//{out << "<"; copy_ -> print (out, descend); out << "><" << copy_ << ">"; }


// /// redirect variables to proper variable vector
// void exprClone::realign (const CouenneProblem *p) {

//   if (((copy_ -> Type () == VAR) ||
//        (copy_ -> Type () == AUX)) &&
//       (copy_ -> Original () != p -> Var (copy_ -> Index ()))) {

//     expression *trash = copy_;

//     copy_ = p -> Var (copy_ -> Index ());
//     printf ("deleting "); trash -> print (); printf ("\n");
//     delete trash;
//   } else {printf ("not deleted "); print (); printf ("\n");}
// }


/// Closest feasible points in function in both directions
void expression::closestFeasible (expression *varind, expression *vardep,
				  CouNumber& left, CouNumber& right) const
{
  assert(isBijective());
  CouNumber inv = inverse(vardep);
  CouNumber curr = (*varind) ();
  if (curr > inv) {
    left  = inv;
    right = curr;
  }
  else {
    left  = curr;
    right = inv;
  }
}
