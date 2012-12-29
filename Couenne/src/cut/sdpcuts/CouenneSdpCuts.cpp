/* $Id$
 *
 * Name:    CouenneSdpCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: wrapper for Couenne to insert sdpcuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "IpOptionsList.hpp"

#include "CouenneJournalist.hpp"
#include "CouenneMatrix.hpp"
#include "CouennePSDcon.hpp"
#include "CouenneSdpCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

#include "CoinTime.hpp"

using namespace Couenne;

//#define DEBUG

/// Constructor
CouenneSdpCuts::CouenneSdpCuts (CouenneProblem *p,
				JnlstPtr jnlst,
				const Ipopt::SmartPtr <Ipopt::OptionsList> options):
  problem_  (p),
  doNotUse_ (false) {

  std::string s;

  options -> GetIntegerValue ("sdp_cuts_num_ev", numEigVec_, "couenne.");
  options -> GetStringValue  ("sdp_cuts_neg_ev", s,          "couenne.");
  onlyNegEV_ = (s == "yes");

  CouenneExprMatrix *cauldron = new CouenneExprMatrix;

  // 1) Construct matrix with entries x_i, x_j

  for (std::vector <exprVar *>::iterator 
	 i  = p -> Variables (). begin (); 
       i   != p -> Variables (). end   (); ++i)

    if ((*i) -> Type () == AUX) {

      expression *image = (*i) -> Image ();

      /// either it is x_i x_j
      if ((image -> code () == COU_EXPRMUL) && 
	  (image -> ArgList () [0] -> Type () != CONST) && 
	  (image -> ArgList () [1] -> Type () != CONST)) {

	int 
	  index0 = image -> ArgList () [0] -> Index (),
	  index1 = image -> ArgList () [1] -> Index ();

	if ((index0 >= 0) &&
	    (index1 >= 0) &&
	    ((*i) -> Index () >= 0))

	  // cauldron -> add_element (CoinMin (index0, index1),
	  //  			      CoinMax (index0, index1), (*i));

	  cauldron -> add_element (index0, index1, (*i));
	  cauldron -> add_element (index1, index0, (*i));
      }

      /// or it is x_i ^ 2
      if ((image -> code () == COU_EXPRPOW) && 
	  (image -> ArgList () [0] -> Type  () != CONST) && 
	  (image -> ArgList () [1] -> Value () == 2.)) {

	int index0 = image -> ArgList () [0] -> Index ();

	if ((index0 >= 0) &&
	    ((*i) -> Index () >= 0))

	  cauldron -> add_element (index0, index0, (*i));
      }
    }

#ifdef DEBUG
  printf ("Cauldron so far\n");
  cauldron -> print ();
#endif

  // TODO

  // 2) Block-partition it (optional), obtain matrices. Replace line
  //    below to decompose cauldron

  minors_ . push_back (cauldron);

  // 3) Bottom-right border each block with a row vector, a column vector,
  // and the constant 1

  for (std::vector <CouenneExprMatrix *>::iterator 
	 i  = minors_ . begin ();
       i   != minors_ . end   (); ++i) {

    int size = problem_ -> nVars ();

#ifdef DEBUG
    printf ("minor has %ld rows and %ld columns\n",
    	    (*i) -> getRows () . size (),
    	    (*i) -> getCols () . size ());
#endif

    for (std::set <std::pair <int, CouenneSparseVector *> >::iterator 
	   j  = (*i) -> getCols () . begin (); 
	 j   != (*i) -> getCols () . end   (); ++j)

      (*i) -> varIndices () . push_back (problem_ -> Var (j -> first));

    for (std::vector <expression *>::iterator 
	   j  = (*i) -> varIndices () . begin ();
	 j   != (*i) -> varIndices () . end   (); ++j) {

      int indexVar = (*j) -> Index ();
#ifdef DEBUG
      printf ("adding at [%d,%d] and viceversa\n", indexVar, size);
#endif
      (*i) -> add_element (indexVar, size, *j); // note: problem_ -> Var (indexVar) = (*j)
      (*i) -> add_element (size, indexVar, *j);
    }

    (*i) -> add_element (size, size, new exprConst (1.));

#ifdef DEBUG
    (*i) -> print ();
#endif
  }

  // 0) Search for X \succeq 0 constraints, if any, then add matrix to
  //    minors for each such constraint

  if (p -> ConstraintClass ("PSDcon"))
    for (std::vector <CouenneConstraint *>::iterator 
	   i  = p -> ConstraintClass ("PSDcon") -> begin (); 
	 i   != p -> ConstraintClass ("PSDcon") -> end   (); ++i) {

      CouennePSDcon *con = dynamic_cast <CouennePSDcon *> (*i);

      if (!con) 
	continue;

      minors_ . push_back (con -> getX ());
    }
}


/// Destructor
CouenneSdpCuts::~CouenneSdpCuts () {

  for (std::vector <CouenneExprMatrix *>::iterator 
  	 i  = minors_ . begin ();
       i   != minors_ . end   (); ++i)

    delete (*i);

  // Destroy matrix structures
}


/// Copy constructor
CouenneSdpCuts::CouenneSdpCuts (const CouenneSdpCuts &rhs):

  problem_   (rhs. problem_),
  numEigVec_ (rhs. numEigVec_),
  onlyNegEV_ (rhs. onlyNegEV_) {

  for (std::vector <CouenneExprMatrix *>::const_iterator 
  	 i  = rhs.minors_ . begin ();
       i   != rhs.minors_ . end   (); ++i)

    minors_ . push_back (new CouenneExprMatrix (**i));
}


/// Assignment
CouenneSdpCuts &CouenneSdpCuts::operator= (const CouenneSdpCuts &rhs) {

  problem_   = rhs. problem_;
  numEigVec_ = rhs. numEigVec_;
  onlyNegEV_ = rhs. onlyNegEV_;

  for (std::vector <CouenneExprMatrix *>::const_iterator 
  	 i  = rhs.minors_ . begin ();
       i   != rhs.minors_ . end   (); ++i)

    minors_ . push_back (new CouenneExprMatrix (**i));

  return *this;
}


/// clone constructor
CglCutGenerator *CouenneSdpCuts::clone () const
{return new CouenneSdpCuts (*this);}


/// Add list of options to be read from file
void CouenneSdpCuts::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("sdp_cuts",
     "The frequency (in terms of nodes) at which Couenne SDP cuts are generated.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. \
Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. \
A negative number -n means that generation should be attempted at the root node, and if successful it can be repeated at every n nodes, otherwise it is stopped altogether."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("sdp_cuts_num_ev",
     "The number of eigenvectors of matrix X to be used to create sdp cuts.",
     -1, -1,
     "Set to -1 to indicate that all n eigenvectors should be used. Eigenvalues are \
sorted in non-decreasing order, hence selecting 1 will provide cuts on the most negative eigenvalue"
    );

  roptions -> AddStringOption2
    ("sdp_cuts_neg_ev",
     "Only use negative eigenvalues to create sdp cuts.",
     "yes", 
     "no", "use all eigenvalues regardless of their sign",
     "yes", "exclude all non-negative eigenvalues"
    );

  /////////////////////////////////////////

  roptions -> AddStringOption2
    ("sdp_cuts_sparsify",
     "Make cuts sparse by greedily reducing X one column at a time before extracting eigenvectors.",
     "yes", 
     "no", "",
     "yes", ""
    );

#if 0
  roptions -> AddStringOption2
    ("sdp_cuts_neg_ev",
     "Only use negative eigenvalues to create sdp cuts.",
     "yes", 
     "no", "use all eigenvalues regardless of their sign",
     "yes", "exclude all non-negative eigenvalues"
    );
#endif
}
