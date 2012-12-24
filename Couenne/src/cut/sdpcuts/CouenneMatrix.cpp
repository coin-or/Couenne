/* $Id$
 *
 * Name:    CouenneMatrix.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of the class of expression matrices
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

#include <stdio.h>

#include "CouenneMatrix.hpp"
#include "CouenneExprConst.hpp"

using namespace Couenne;


// copy constructor
CouenneSparseVector::CouenneSparseVector (const CouenneSparseVector &rhs) {

  for (register std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::iterator 
	 i  = rhs. elem_. begin ();
       i   != rhs. elem_. end   (); ++i)

    elem_. insert (new CouenneScalar (**i));
}


// assignment operator
CouenneSparseVector &CouenneSparseVector::operator= (const CouenneSparseVector &rhs) {

  for (register std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::iterator 
	 i  = rhs. elem_. begin ();
       i   != rhs. elem_. end   (); ++i)

    elem_. insert (new CouenneScalar (**i));

  return *this;
}

#define copy_vectors(from, to) {                                                                                         \
  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator                  \
	 rowIt  = from. begin ();                                                                                        \
       rowIt   != from. end   (); ++rowIt) {                                                                             \
    to . insert (std::pair <int, CouenneSparseVector *> (rowIt -> first, new CouenneSparseVector (*(rowIt -> second)))); \
  }                                                                                                                      \
}

/// copy constructor
CouenneExprMatrix::CouenneExprMatrix (const CouenneExprMatrix &rhs):
  varIndices_ (rhs.varIndices_) {

  copy_vectors(rhs.row_, row_);
  copy_vectors(rhs.col_, col_);
}

/// Assignment operator
CouenneExprMatrix &CouenneExprMatrix::operator= (const CouenneExprMatrix &rhs) {

  varIndices_ = rhs.varIndices_;

  copy_vectors(rhs.row_, row_);
  copy_vectors(rhs.col_, col_);

  return *this;
}

void CouenneScalar::print () const {
  printf ("[%d,", index_); 
  if (elem_) 
    elem_ -> print (); 
  printf ("]");
}


/// Insertion into vector
void CouenneSparseVector::add_element (int index, expression *elem) {

  CouenneScalar *element = new CouenneScalar (index, elem);
  elem_ . insert (element);
}


/// return size of (square sub-) matrix
long unsigned int CouenneExprMatrix::size () 
{return CoinMax (row_ . size (), col_ . size ());}


/// Insertion into matrix
void CouenneExprMatrix::add_element (int rowInd, int colInd, expression *elem) {

  // don't duplicate code, macroize

#define check_and_insert(indMaj,indMin,vecMaj,elem)                                                  \
  {									        	 	     \
    std::          pair <int, CouenneSparseVector *> findme (indMaj, NULL);  	        	     \
    std::set <std::pair <int, CouenneSparseVector *>,                                                \
	      CouenneExprMatrix::compare_pair_ind>::const_iterator check = vecMaj.find (findme);     \
                                                                                                     \
    if (check == vecMaj. end ()) {					        	 	     \
      std::pair <int, CouenneSparseVector *> new_vector (indMaj, new CouenneSparseVector);           \
      new_vector.second -> add_element (indMin, elem);			        	 	     \
      vecMaj. insert (new_vector);					        	 	     \
    } else check -> second -> add_element (indMin, elem);	                                     \
  }

  check_and_insert (rowInd, colInd, row_, elem);
  if (elem -> code () == COU_EXPRCONST) 
    elem = new exprClone (elem);
  check_and_insert (colInd, rowInd, col_, elem);
}


/// Dot product
inline double CouenneSparseVector::operator * (const CouenneSparseVector &v2) const 
{return multiply_thres (v2, COIN_DBL_MAX);}


/// Threshold dot product
double CouenneSparseVector::multiply_thres (const CouenneSparseVector &v2, double thres) const {

  double prod = 0.;

  for (register std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::iterator 
	 i1 =    elem_. begin (),
	 i2 = v2.elem_. begin ();

	 ((i1 !=    elem_.end ()) && 
	  (i2 != v2.elem_.end ()));) {

    while ((i1 !=    elem_.end ()) && ((*i1) -> getIndex () < (*i2) -> getIndex ())) ++i1; if (i1 ==     elem_. end ()) return prod;
    while ((i2 != v2.elem_.end ()) && ((*i2) -> getIndex () < (*i1) -> getIndex ())) ++i2; if (i2 == v2. elem_. end ()) return prod;

    prod += 
      (*((*i1) -> getElem ())) () * 
      (*((*i2) -> getElem ())) ();

    if (prod > thres) 
      break;
  }

  return prod;
}


/// vector * matrix
CouenneSparseVector &CouenneSparseVector::operator * (const CouenneExprMatrix &post) const {

  CouenneSparseVector *product = new CouenneSparseVector;

  const std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind> &columns = post. getCols ();

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator colIt = columns. begin (); 
       colIt != columns. end (); ++colIt) {

    double single = operator* (*(colIt -> second));

    if (single != 0.)
      product -> add_element (colIt -> first, new exprConst (single));
  }

  return *product;
}


/// matrix * vector
CouenneSparseVector &CouenneExprMatrix::operator * (const CouenneSparseVector &post) const {

  CouenneSparseVector *product = new CouenneSparseVector;

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator rowIt = row_. begin (); 
       rowIt != row_. end (); ++rowIt) {

    double single = post. operator* (*(rowIt -> second));

    if (single != 0.)
      product -> add_element (rowIt -> first, new exprConst (single));
  }

  return *product;
}


/// matrix * matrix
CouenneExprMatrix &CouenneExprMatrix::operator * (const CouenneExprMatrix &post) const {

  CouenneExprMatrix *product = new CouenneExprMatrix;
  return *product;
}

/// Destructor
CouenneExprMatrix::~CouenneExprMatrix () {

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator 
	 i  = row_ . begin ();
       i   != row_ . end   (); ++i)

    delete i -> second;

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator 
	 i  = col_ . begin ();
       i   != col_ . end   (); ++i)

    delete i -> second;
}


/// Destructor
CouenneSparseVector::~CouenneSparseVector () {

  for (register std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::iterator 
	 i  = elem_. begin ();
       i   != elem_. end   (); ++i)
    delete (*i);
}


CouenneScalar::~CouenneScalar () {

  if (delete_) // only delete constants
    delete elem_;
} 

#define WRAP 20

/// Pretty print
void CouenneSparseVector::print () const {

  int cnt=0;

  printf ("Vector (%ld) (", elem_ .  size ());

  for (std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::iterator i = elem_ . begin (); 
       i != elem_ . end (); ++i) {

    if (i != elem_ . begin ())
      printf (",");

    (*i) -> print ();

    if (!(++cnt) % WRAP)
      printf ("\n   ");
  }

  printf (")");
}


/// Pretty print
void CouenneExprMatrix::print () const {

  printf ("Matrix (%ld x %ld):\n", 
	  row_ . size (), 
	  col_ . size ());

  // print rows

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator 
	 i  = row_ . begin ();
       i   != row_ . end (); ++i) {

    printf ("Row [%d]: ", (*i) . first);
    (*i) . second -> print ();
    printf ("\n");
  }

  // print columns

  for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::iterator 
	 i  = col_ . begin ();
       i   != col_ . end (); ++i) {

    printf ("Col [%d]: ", (*i) . first);
    (*i) . second -> print ();
    printf ("\n");
  }

  if (varIndices_ . size () > 0) {
    printf ("varIndices: (");
    for (std::vector <expression *>::const_iterator 
	   i  = varIndices_ . begin (); 
	 i   != varIndices_ . end   (); ++i) {
      if (i != varIndices_ . begin ())
	printf (",");
      (*i) -> print ();
    }
    printf (")\n");
  }
}
