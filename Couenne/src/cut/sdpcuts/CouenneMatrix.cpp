/* $Id$
 *
 * Name:    CouenneMatrix.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of the class of expression matrices
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinFinite.hpp"

#include "CouenneMatrix.hpp"
#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"

using namespace Couenne;

/// Insertion into vector
void CouenneSparseVector::add_element (int index, expression *elem) {

  CouenneScalar *element = new CouenneScalar (index, elem);
  elem_ . insert (element);
}


/// Insertion into matrix
void CouenneSparseMatrix::add_element (int rowInd, int colInd, expression *elem) {

  // don't duplicate, macroize

#define check_and_insert(indMaj,indMin,vecMaj)                                                       \
  {									        	 	     \
    std::pair <int, CouenneSparseVector *> findme (indMaj, NULL);  	        	 	     \
    std::set <std::pair <int, CouenneSparseVector *> >::const_iterator check = vecMaj.find (findme); \
                                                                                                     \
    if (check == vecMaj. end ()) {					        	 	     \
      std::pair <int, CouenneSparseVector *> new_vector (indMaj, new CouenneSparseVector);           \
      new_vector.second -> add_element (indMin, elem);			        	 	     \
      vecMaj. insert (new_vector);					        	 	     \
    } else check -> second -> add_element (indMin, elem);	                                     \
  }
  
  check_and_insert (rowInd, colInd, row_);
  check_and_insert (colInd, rowInd, col_);
}


/// Dot product
inline double CouenneSparseVector::operator * (const CouenneSparseVector &v2) const 
{return multiply_thres (v2, COIN_DBL_MAX);}


/// Threshold dot product
double CouenneSparseVector::multiply_thres (const CouenneSparseVector &v2, double thres) const {

  double prod = 0.;

  for (register std::set <CouenneScalar *>::iterator 
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
CouenneSparseVector &CouenneSparseVector::operator * (const CouenneSparseMatrix &post) const {

  CouenneSparseVector *product = new CouenneSparseVector;

  const std::set <std::pair <int, CouenneSparseVector *> > &columns = post. getColumns ();

  for (std::set <std::pair <int, CouenneSparseVector *> >::iterator colIt = columns. begin (); 
       colIt != columns. end (); ++colIt) {

    double single = operator* (*(colIt -> second));

    if (single != 0.)
      product -> add_element (colIt -> first, new exprConst (single));
  }

  return *product;
}


/// matrix * vector
CouenneSparseVector &CouenneSparseMatrix::operator * (const CouenneSparseVector &post) const {

  CouenneSparseVector *product = new CouenneSparseVector;

  for (std::set <std::pair <int, CouenneSparseVector *> >::iterator rowIt = row_. begin (); 
       rowIt != row_. end (); ++rowIt) {

    double single = post. operator* (*(rowIt -> second));

    if (single != 0.)
      product -> add_element (rowIt -> first, new exprConst (single));
  }

  return *product;
}


// /// matrix * matrix
// CouenneSparseVector <T> &CouenneSparseVector::operator * (CouenneSparseMatrix <T> &post) const {
//   CouenneSparseVector <T> *product = new CouenneSparseVector <T>;
//   product -> add_element ();
// }


#define WRAP 20

// /*****************************************************************************************/
// template <class T> void CouenneSparseVector::print () {

//   int i;

//   if (!v) return;

//   printf ("Vector(%d): ", v -> n);

//   for (i=0; i<v -> n; i++) {
//     printf ("[%d %g]", v -> indices [i], v -> values [i]);
//     if (!((i+1) % WRAP))
//       printf ("\n   ");
//   }

//   printf ("\n");
// }


// /*****************************************************************************************/
// template <class T> void CouenneSparseMatrix::print () {

//   int i;

//   printf ("matrix (%d,%d):\nRows\n", m -> nrow, m -> ncol);

//   for (i=0; i<m -> nrow; i++)
//     print (m -> row [i]);

//   printf ("\nColumns\n");

//   for (i=0; i<m -> ncol; i++)
//     m -> col [i] -> print ();

//   printf ("\n");
// }
