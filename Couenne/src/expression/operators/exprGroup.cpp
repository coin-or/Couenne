/* $Id$
 *
 * Name:    exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprGroup
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneExprConst.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneProblem.hpp"

#include <cassert>

using namespace Couenne;

namespace Couenne {
class Domain;
}

// eliminates elements with zero coefficients
void cleanZeros (std::vector <std::pair <exprVar *, CouNumber> > &lcoeff) {

  std::vector <std::pair <exprVar *, CouNumber> >::iterator i = lcoeff.begin ();

  int    ind  = 0;
  size_t size = lcoeff.size ();
  
  while (size-- > 0) {
    if ((i -> second ==  0.) || 
	(i -> second == -0.)) {
      lcoeff.erase (i);
      i = lcoeff.begin () + ind;
    } else {
      ++i;
      ++ind;
    }
  }
}



/// Generalized (static) constructor: check parameters and return a
/// constant, a single variable, or a real exprGroup
expression *exprGroup::genExprGroup (CouNumber c0,
				     std::vector <std::pair <exprVar *, CouNumber> > &lcoeff, 
				     expression **al, 
				     int n) {
  size_t nl = lcoeff.size ();
  expression *ret = NULL;

  cleanZeros (lcoeff);

  // a constant
  if ((n==0) && (nl==0))
    ret = new exprConst (c0); // a constant auxiliary? FIX!

  else if ((n==0) && (fabs (c0) < COUENNE_EPS) && (nl==1)) { // a linear monomial, cx

    if (fabs (lcoeff[0]. second - 1) < COUENNE_EPS)
      ret    = new exprClone (lcoeff[0]. first);
    else ret = new exprMul (new exprConst (lcoeff[0]. second), new exprClone (lcoeff[0]. first));
    
  } else ret = new exprGroup (c0, lcoeff, al, n);

  return ret;
}


/// Constructor
exprGroup::exprGroup (CouNumber c0,
		      std::vector <std::pair <exprVar *, CouNumber> > &lcoeff, 
		      expression **al, 
		      int n):
  exprSum  (al, n),
  lcoeff_  (lcoeff),
  c0_      (c0) {

  cleanZeros (lcoeff_);
}


/// copy constructor
exprGroup::exprGroup  (const exprGroup &src, Domain *d): 
  exprSum   (src.clonearglist (d), src.nargs_),
  c0_       (src.c0_) {

  for (lincoeff::iterator i = src.lcoeff_.begin (); i != src.lcoeff_.end (); ++i)

    lcoeff_ . push_back (std::pair <exprVar *, CouNumber> 
			 //(dynamic_cast <exprVar *> (i -> first -> clone (d)), i -> second));
			 (new exprVar (i -> first -> Index (), d), i -> second));
}


/// Destructor -- check if there are exprBounds and delete them 
exprGroup::~exprGroup () {

  for (lincoeff::iterator i = lcoeff_.begin (); i != lcoeff_.end (); ++i) {
    enum expr_type code = i -> first -> code ();
    if ((code == COU_EXPRLBOUND) || 
	(code == COU_EXPRUBOUND))
      delete i -> first;
  }
}


/// I/O
void exprGroup::print (std::ostream &out, bool descend) const {

  //if (code () == COU_EXPRGROUP)
  if (lcoeff_.size () > 0)
    out << '(';
  
  bool nzNL = nargs_ && ((nargs_ > 1) ||
			 ((*arglist_) -> Type () != CONST) ||
			 (fabs ((*arglist_) -> Value ()) > COUENNE_EPS));
    
  if (nzNL)
    exprSum::print (out, descend);

  if      (c0_ >   0.) {if (nzNL) out << '+'; out << c0_;}
  else if (c0_ < - 0.)                        out << c0_;

  for (size_t n = lcoeff_.size (), i=0; n--; i++) {

    CouNumber coeff = lcoeff_ [i]. second;

    if      (coeff >   0.) { if (i || (c0_ != 0.) || nzNL) out << '+'; if (coeff !=  1.) out <<  coeff << "*";}
    else if (coeff < - 0.) {                               out << '-'; if (coeff != -1.) out << -coeff << "*";}

    lcoeff_ [i]. first -> print (out, descend);
    if (!((i + 1) % MAX_ARG_LINE) && n)
      out << std::endl;
  }

  //if (code () == COU_EXPRGROUP)
  if (lcoeff_.size () > 0)
    out << ')';
}


/// differentiation
expression *exprGroup::differentiate (int index) {

  expression **arglist = new expression * [nargs_ + 1];

  CouNumber totlin=0;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    if (el -> first -> Index () == index)
      totlin += el -> second;

  int nargs = 0;

  if (fabs (totlin) > COUENNE_EPS)
    arglist [nargs++] = new exprConst (totlin);

  for (int i = 0; i < nargs_; i++) 
    if (arglist_ [i] -> dependsOn (&index, 1))
      arglist [nargs++] = arglist_ [i] -> differentiate (index);

  if ((nargs == 0) ||
      ((nargs == 1) && (fabs (totlin) > COUENNE_EPS))) {
    delete [] arglist;
    return new exprConst (totlin);
  }
  else return new exprSum (arglist, nargs);
}


/// get a measure of "how linear" the expression is:
int exprGroup::Linearity () {

  int 
    nllin = exprSum::Linearity (),    // linearity of nonlinear part
    llin  = (lcoeff_.size () == 0) ?  //              linear part
    ((fabs (c0_) < COUENNE_EPS) ? ZERO : CONSTANT) : 
    LINEAR;

  return (llin > nllin) ? llin : nllin;
}


/// compare affine terms
int exprGroup::compare (exprGroup &e) {

  // !!! why?

  //int sum = exprSum::compare (e);

  //if (sum != 0) 
  //return sum;

  if (c0_ < e.c0_ - COUENNE_EPS) return -1;
  if (c0_ > e.c0_ + COUENNE_EPS) return  1;

  if (lcoeff_.size () < e.lcoeff_.size ()) return -1;
  if (lcoeff_.size () > e.lcoeff_.size ()) return  1;

  for (lincoeff::iterator 
	 el1 =   lcoeff_.begin (),
	 el2 = e.lcoeff_.begin ();
       el1 != lcoeff_.end (); 
       ++el1, ++el2) {

    int 
      ind1 = el1 -> first -> Index (),
      ind2 = el2 -> first -> Index ();

    CouNumber 
      coe1 = el1 -> second,
      coe2 = el2 -> second;

    if (ind1 < ind2) return -1;
    if (ind1 > ind2) return  1;

    if (coe1 < coe2 - COUENNE_EPS) return -1;
    if (coe1 > coe2 + COUENNE_EPS) return  1;
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprGroup::rank () {

  int maxrank = exprOp::rank ();

  if (maxrank < 0) 
    maxrank = 0;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    int r = el -> first -> rank ();
    if (r > maxrank)
      maxrank = r;
  }

  return maxrank;
}


/// update dependence set with index of this variable
void exprGroup::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {

  exprOp::fillDepSet (dep, g);

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    dep -> insert (g -> lookup (el -> first -> Index ()));
}


/// fill in the set with all indices of variables appearing in the
/// expression
int exprGroup::DepList (std::set <int> &deplist,
			enum dig_type type) {

  int deps = exprOp::DepList (deplist, type);

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    /*printf ("before ["); 
    el -> first -> print ();
    printf ("]: {");
    for (std::set <int>::iterator i=deplist.begin (); i != deplist.end(); ++i)
    printf ("%d ", *i);*/

    deps += el -> first -> DepList (deplist, type);

    /*printf ("}, after: {");
    for (std::set <int>::iterator i=deplist.begin (); i != deplist.end(); ++i)
      printf ("%d ", *i);
      printf ("}\n");*/
  }

  return deps;
}


/// is this linear term integer?
bool exprGroup::isInteger () {

  if (!(::isInteger (c0_)) ||
      !(exprOp::isInteger ()))
    return false;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    CouNumber coe = el -> second;

    bool
      intCoe = ::isInteger (coe),
      intVar = el -> first -> isInteger ();

    if (intCoe && intVar)
      continue;

    CouNumber 
      lb = el -> first -> lb (), 
      ub = el -> first -> ub ();

    // check var fixed and product is integer
    if ((fabs (lb - ub) < COUENNE_EPS) &&
	(::isInteger (lb * coe) ||
	 (intCoe && ::isInteger (lb)))) 
      continue;

    return false;
  }

  return true;
}


/// replace variable x with new (aux) w
void exprGroup::replace (exprVar *x, exprVar *w) {

  exprOp::replace (x, w);

  int 
    xind = x -> Index (),
    wind = w -> Index ();

  // find occurrences of x and w in vector of variabls

  lincoeff::iterator x_occur = lcoeff_.begin ();

  // Do not assume index vector is sorted in ascending order
  // w.r.t. (*iterator) -> first () -> Index()
  while ((x_occur != lcoeff_.end ()) && 
	 (x_occur -> first -> Index () != xind))
    ++x_occur;

  if (x_occur == lcoeff_ .end ()) // not found, bail out
    return;

  if (xind == wind)
    x_occur -> first = w;
  else { // replacing two variables, not the features of one

    lincoeff::iterator w_occur = lcoeff_.begin ();

    // Do not assume index vector is sorted in ascending order
    // w.r.t. (*iterator) -> first () -> Index()
    while ((w_occur != lcoeff_.end ()) &&
	   (w_occur -> first -> Index () != wind))
      ++w_occur;

    if (w_occur == lcoeff_ . end ()) // not found w, simply substitute 
      x_occur -> first = w;
    else {
		if ((w_occur -> second += x_occur -> second) == 0.) { // add coefficients
	lcoeff_.erase (w_occur);                // if they cancel out, delete w as well

	// under Microsoft, x_occur may have been invalidated by removing w_occur from lcoeff_, so we search for it again
	for( x_occur = lcoeff_.begin (); x_occur -> first -> Index () != xind; ++x_occur )
		assert(x_occur != lcoeff_ .end ()); // it was found before, so it should be still found
	}
      lcoeff_.erase   (x_occur);                // delete entry of x
    }
  }
}


/// return l-2 norm of gradient at given point. Not needed for now, as
/// we only use it with nonlinear operators
CouNumber exprGroup::gradientNorm (const double *x) {

  CouNumber retval = 0;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    retval += el -> second * el -> second;

  return sqrt (retval);
}


/// Redirect variables to proper variable vector
void exprGroup::realign (const CouenneProblem *p) {

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    exprVar *var = el -> first;

    if (((var -> Type () == VAR) ||  
	 (var -> Type () == AUX)) &&
	(var -> Original () != p -> Var (var -> Index ()))) {

      expression *trash = var;
      el -> first = p -> Var (var -> Index ());
      delete trash;
    }
  }
}


/// simplification
expression *exprGroup::simplify () {
  exprOp::simplify (); 
  //if (lcoeff_. size () <= 0) // this is just a constant
  //return new exprConst (c0_);
  //else
  return NULL;
}

