/*
 *
 * Name:    exprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprQuad
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneLQelems.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

using namespace Couenne;

namespace Couenne {
class Domain;
}

//#define DEBUG

struct cmpVar {
  bool operator() (const exprVar* v1, const exprVar* v2) const
  {return (v1 -> Index () < v2 -> Index ());}
};

/// Constructor
exprQuad::exprQuad (CouNumber c0,
		    std::vector <std::pair <exprVar *, CouNumber> > &lcoeff,
		    std::vector <quadElem> &qcoeff,
		    expression **al,
		    int n):

  exprGroup (c0, lcoeff, al, n) {

  nqterms_ = 0;

  typedef std::map <exprVar *, CouNumber, cmpVar> rowMap;
  typedef std::map <exprVar *, rowMap,    cmpVar> matrixMap;

  matrixMap qMap;

  for (std::vector <quadElem>::iterator qel = qcoeff.begin (); qel != qcoeff.end (); ++qel) {

    CouNumber coe = qel -> coeff ();

    exprVar
      *varI = qel -> varI (),
      *varJ = qel -> varJ ();

    if (varI -> Index () != varJ -> Index ()) {

      //coe /= 2.;

      // pick smaller index as row reference
      if (varI -> Index () > varJ -> Index ()) {

	exprVar *swap = varJ;
	varJ = varI;
	varI = swap;
      }
    }

    matrixMap::iterator rowp = qMap.find (varI);

    if (rowp == qMap.end ()) { // add new row

      std::pair <exprVar *, CouNumber> newcell (varJ, coe);
      rowMap rmap;
      rmap.insert (newcell);

      std::pair <exprVar *, rowMap>    newrow  (varI, rmap);
      qMap.insert (newrow);

    } else { // insert element into row

      rowMap::iterator cell = rowp -> second.find (varJ);

      if (cell == rowp -> second.end ()) { // normal case, add entry

	std::pair <exprVar *, CouNumber> newcell (varJ, coe);
	rowp -> second.insert (newcell);

      } else { // strange, but add coefficient

	if (fabs (cell -> second += coe) < COUENNE_EPS)
	  // eliminate element of map if null coefficient
	  rowp -> second.erase (cell);
      }
    }
  }

  // transform maps into vectors

  for (matrixMap::iterator row = qMap.begin (); row != qMap.end (); ++row) {

    sparseQcol line;

    // insert first element in bound map
    if (bounds_.find (row -> first) == bounds_.end ()) {

      std::pair <CouNumber, CouNumber> newbound (-COIN_DBL_MAX, COIN_DBL_MAX);
      std::pair <exprVar *, std::pair <CouNumber, CouNumber> > newvar (row -> first, newbound);
      bounds_.insert (newvar);
    }

    for (rowMap::iterator cell = row -> second.begin (); cell != row -> second.end (); ++cell) {

      line.push_back (std::pair <exprVar *, CouNumber> (*cell));

      // insert second element in bound map
      if (bounds_.find (cell -> first) == bounds_.end ()) {

	std::pair <CouNumber, CouNumber> newbound (-COIN_DBL_MAX, COIN_DBL_MAX);
	std::pair <exprVar *, std::pair <CouNumber, CouNumber> > newvar (cell -> first, newbound);
	bounds_.insert (newvar);
      }
    }

    matrix_.push_back (std::pair <exprVar *, sparseQcol> (row -> first, line));
    nqterms_ += (int) (line.size ());
  }
}


/// copy constructor
exprQuad::exprQuad (const exprQuad &src, Domain *d):
  exprGroup (src, d),
  bounds_   (src.bounds_),
  nqterms_  (src.nqterms_) {

  for (sparseQ::iterator row = src.matrix_.begin (); row != src.matrix_ . end (); ++row) {

    sparseQcol column;

    for (sparseQcol::iterator i = row -> second. begin (); i != row -> second. end (); ++i)
      column.push_back (std::pair <exprVar *, CouNumber>
			//(dynamic_cast <exprVar *> (i -> first -> clone (d)), i -> second));
			(new exprVar (i -> first -> Index (), d), i -> second));

    matrix_.push_back (std::pair <exprVar *, sparseQcol>
		       //dynamic_cast <exprVar *> (row -> first -> clone (d)), column));
		       (new exprVar (row -> first -> Index (), d), column));
  }

  //////////////////////////////////////////////////////////////////////////////

  std::vector
    <std::pair <CouNumber, std::vector
    <std::pair <exprVar *, CouNumber> > > >::iterator row;

  for (row = src.eigen_ . begin ();
       row != src.eigen_ . end (); ++row) {

    std::vector <std::pair <exprVar *, CouNumber> > eigVec;

    for (std::vector <std::pair <exprVar *, CouNumber> >::iterator
	   i = row -> second. begin ();
	 i != row -> second. end (); ++i)
      eigVec.push_back (std::pair <exprVar *, CouNumber>
			(dynamic_cast <exprVar *> (i -> first -> clone (d)), i -> second));

    eigen_.push_back (std::pair <CouNumber, std::vector
		       <std::pair <exprVar *, CouNumber> > > (row -> first, eigVec));
  }
}


/// I/O
void exprQuad::print (std::ostream &out, bool descend) const {

  //if (code () == COU_EXPRQUAD)
  if (matrix_.size () > 0)
    out << '(';

  // print linear and nonquadratic part
  exprGroup::print (out, descend);

  int noperands = 0;

  for (size_t n = matrix_.size (), i=0; n--; i++) {
    //sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    int xind = matrix_ [i].first -> Index ();
    const sparseQcol row = matrix_ [i].second;

    for (int m = row.size (), j=0; m--; j++) {
      //sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      if (fabs (row [j]. second - 1.) > COUENNE_EPS) {
	if (fabs (row [j]. second + 1.) < COUENNE_EPS) out << "- ";
	else {
	  if (row [j]. second > 0.) out << '+';
	  out << row [j]. second << "*";
	}
      } else out << '+';

      if (row [j].first -> Index () == xind) {
	matrix_ [i]. first -> print (out, descend);
	out << "^2";
      } else {
	matrix_ [i]. first -> print (out, descend);
	out << '*';
	row [j]. first -> print (out, descend);
      }

      if (!((noperands + 1) % MAX_ARG_LINE))
	out << std::endl;
    }
  }

  //if (code () == COU_EXPRGROUP)
  if (matrix_.size () > 0)
    out << ')';
}


/// differentiation
expression *exprQuad::differentiate (int index) {

  std::map <exprVar *, CouNumber> lmap;

  CouNumber c0 = 0;

  // derive linear part (obtain constant)
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    c0 += el -> second;

  // derive quadratic part (obtain linear part)
  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    int xind = row -> first -> Index ();

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      int yind = col -> first -> Index ();

      CouNumber coe = col -> second;
      exprVar *var = col -> first;

      if      (xind == index)
	if    (yind == index) {var = col -> first; coe *= 2;}
	else                   var = col -> first;
      else if (yind == index)  var = row -> first;
      else continue;

      std::map <exprVar *, CouNumber>::iterator i = lmap.find (var);

      if (i != lmap.end()) {
	if (fabs (i -> second += coe) < COUENNE_EPS)
	  lmap.erase (i);
      } else {
	std::pair <exprVar *, CouNumber> npair (var, coe);
	lmap.insert (npair);
      }
    }
  }

  // derive nonlinear sum
  expression **arglist = new expression * [nargs_ + 1];
  int nargs = 0;

  for (int i = 0; i < nargs_; i++)
    if (arglist_ [i] -> dependsOn (index))
      arglist [nargs++] = arglist_ [i] -> differentiate (index);

  // special cases

  // 1) no linear part
  if (lmap.empty ()) {

    // and no nonlinear part either
    if (!nargs) {
      delete arglist;
      return new exprConst (c0);
    }

    if (fabs (c0) > COUENNE_EPS)
      arglist [nargs++] = new exprConst (c0);

    return new exprSum (arglist, nargs);
  }

  lincoeff coe;

  for (std::map <exprVar *, CouNumber>::iterator i = lmap.begin (); i != lmap.end (); ++i)
    coe.push_back (std::pair <exprVar *, CouNumber> (i -> first, i -> second));

  return new exprGroup (c0, coe, arglist, nargs);
}


/// compare quadratic terms

int exprQuad::compare (exprQuad &e) {

  int sum = exprGroup::compare (e);

  if (sum != 0)
    return sum;

  if (matrix_.size() < e.matrix_.size()) return -1;
  if (matrix_.size() > e.matrix_.size()) return  1;

  for (sparseQ::iterator
	 row1 =   matrix_.begin (),
	 row2 = e.matrix_.begin ();
       row1 != matrix_.end ();
       ++row1, ++row2) {

    if (row1 -> first -> Index () < row2 -> first -> Index ()) return -1;
    if (row1 -> first -> Index () > row2 -> first -> Index ()) return  1;

    if (row1 -> second.size () < row2 -> second.size ()) return -1;
    if (row1 -> second.size () > row2 -> second.size ()) return  1;

    //    if (matrix_.size() > e.matrix_.size()) return  1;
    //    int xind = row -> first -> Index ();
    //    CouNumber x = (*(row -> first)) ();

    for (sparseQcol::iterator
	   col1 = row1 -> second.begin (),
	   col2 = row2 -> second.begin ();
	 col1 != row1 -> second.end ();
	 ++col1, ++col2) {

      if (col1 -> first -> Index () < col2 -> first -> Index ()) return -1;
      if (col1 -> first -> Index () > col2 -> first -> Index ()) return  1;

      if (col1 -> second < col2 -> second - COUENNE_EPS) return -1;
      if (col1 -> second > col2 -> second + COUENNE_EPS) return  1;
    }
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprQuad::rank () {

  int maxrank = exprGroup::rank ();

  if (maxrank < 0)
    maxrank = 0;

  int r;

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    if ((r = row -> first -> rank ()) > maxrank) maxrank = r;

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col)
      if ((r = col -> first -> rank ()) > maxrank) maxrank = r;
  }

  return maxrank;
}


/// update dependence set with index of this variable
void exprQuad::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {

  exprGroup::fillDepSet (dep, g);

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    dep -> insert (g -> lookup (row -> first -> Index ()));

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col)
      dep -> insert (g -> lookup (col -> first -> Index ()));
  }
}


/// fill in the set with all indices of variables appearing in the
/// expression
int exprQuad::DepList (std::set <int> &deplist,
		       enum dig_type type) {

  int deps = exprGroup::DepList (deplist, type);

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {
    deps += row -> first -> DepList (deplist, type);
    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col)
      deps += col -> first -> DepList (deplist, type);
  }

  return deps;
}


/// is this quadratic expression integer?
bool exprQuad::isInteger () {

  if (!(exprGroup::isInteger ()))
    return false;

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    bool intI = row -> first -> isInteger ();

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      CouNumber coe = col -> second;

      bool
	intCoe = ::isInteger (coe),
	intJ   = row -> first -> isInteger ();

      if (intI && intJ && intCoe)
	continue;

      if (!intCoe  // coefficient fractional, check all is fixed and product is integer
	  && row -> first -> isFixed ()
	  && col -> first -> isFixed ()
	  && ::isInteger (coe *
			  row -> first -> lb () *
			  col -> first -> lb ()))
	continue;

      if (!intI && (row -> first -> isFixed ()) && ::isInteger ((*(row -> first)) ())) continue;
      if (!intJ && (col -> first -> isFixed ()) && ::isInteger ((*(col -> first)) ())) continue;

      //if (!intI && !intJ &&  intCoe) ; // check x y fixed int
      //if (!intI &&  intJ &&  intCoe) ; // check x   fixed int
      //if ( intI && !intJ &&  intCoe) ; // check   y fixed int

      return false;
    }
  }

  return true;
}


/// replace variable x with new (aux) w
void exprQuad::replace (exprVar *x, exprVar *w) {

  exprGroup::replace (x, w);
  int xind = x -> Index ();
  int wind = w -> Index ();

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    exprVar * &vr = row -> first;
    if ((vr -> Index () == xind)) {

      //fprintf (stderr, "Didn't fix exprQuad::replace() yet");
      //exit (-1);
      vr = w;
    }

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      exprVar * &vc = col -> first;
      if ((vc -> Index () == wind)) {

	//fprintf (stderr, "Didn't fix exprQuad::replace() yet");
	//exit (-1);
	vc = w;
      }
    }
  }
}


/// Redirect variables to proper variable vector
void exprQuad::realign (const CouenneProblem *p) {

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    exprVar * &vr = row -> first;
    int indVar;

    // substitute variable representing this row with its newest version

    if (((vr -> Type () == VAR) ||
	 (vr -> Type () == AUX)) &&
	(vr -> Original () != p -> Var (indVar = vr -> Index ()))) {

      expression *trash = vr;
      row -> first = p -> Var (indVar);
      delete trash;
    }

    // substitute each variable of this row with its newest version

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      exprVar * &vc = col -> first;
      int indVar;

      // substitute variable representing this row with its newest version

      if (((vc -> Type () == VAR) ||
	   (vc -> Type () == AUX)) &&
	  (vc -> Original () != p -> Var (indVar = vc -> Index ()))) {

	expression *trash = vc;
	col -> first = p -> Var (indVar);
	delete trash;
      }
    }
  }
}


/// compute $y^{lv}$ and $y^{uv}$ for Violation Transfer algorithm
void exprQuad::closestFeasible (expression *varind,
				expression *vardep,
				CouNumber &left,
				CouNumber &right) const {

  fprintf (stderr, "exprQuad::closestFeasible() not available for VT\n");
  exit (-1);
}


/// return l-2 norm of gradient at given point
CouNumber exprQuad::gradientNorm (const double *x) {

  CouNumber grad = 0.;

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    CouNumber gradEl = 0.;
    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col)
      gradEl += col -> second * x [col -> first -> Index ()];

    grad += gradEl * gradEl;
  }

  return sqrt (grad);
}

/// Simplify expression
expression *exprQuad::simplify () {
  exprOp::simplify ();
  return NULL;
}

