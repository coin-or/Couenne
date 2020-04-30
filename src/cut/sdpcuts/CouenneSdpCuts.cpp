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
#include "CouenneExprAux.hpp"
#include "operators/CouenneExprPow.hpp"
#include "operators/CouenneExprMul.hpp"

#include "CoinTime.hpp"

#include <set>

using namespace Couenne;

//#define DEBUG

/// Constructor
CouenneSdpCuts::CouenneSdpCuts (CouenneProblem *p,
				JnlstPtr jnlst,
				const Ipopt::SmartPtr <Ipopt::OptionsList> options):
  problem_  (p),
  doNotUse_ (false) {

  std::string s;

  options -> GetIntegerValue ("sdp_cuts_num_ev",      numEigVec_, "couenne.");
  options -> GetStringValue  ("sdp_cuts_neg_ev",      s,          "couenne."); onlyNegEV_        = (s == "yes");
  options -> GetStringValue  ("sdp_cuts_sparsify",    s,          "couenne."); useSparsity_      = (s == "yes");
  options -> GetStringValue  ("sdp_cuts_fillmissing", s,          "couenne."); fillMissingTerms_ = (s == "yes");

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
	    ((*i) -> Index () >= 0)) {

	  // cauldron -> add_element (CoinMin (index0, index1),
	  //  			      CoinMax (index0, index1), (*i));

	  cauldron -> add_element (index0, index1, (*i));
	  cauldron -> add_element (index1, index0, (*i));
        }
      }

      /// or it is x_i ^ 2
      if ((image -> code () == COU_EXPRPOW) && 
	  (image -> ArgList () [0] -> Type  () != CONST) && 
	  (image -> ArgList () [1] -> Value () == 2.)) {

	int index0 = image -> ArgList () [0] -> Index ();

	if (        (index0   >= 0) &&
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

  // 2.5) If option says so, add fictitious auxiliary variables if not there

  if (fillMissingTerms_) {

    for (std::vector <CouenneExprMatrix *>::iterator 
	   i  = minors_ . begin ();
	 i   != minors_ . end   (); ++i) {

      // First: construct (possibly sparse) index set, to check
      // against each row (whose index set is a SUBSET, if non-proper,
      // of varNumIndices). Do not use varIndices, as it is empty for
      // now.

      std::set <int> varNumIndices;

      for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::const_iterator 
	     rowIt  = (*i) -> getRows (). begin ();
	   rowIt   != (*i) -> getRows (). end   (); ++rowIt) {

	varNumIndices. insert (rowIt -> first);

	for (std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::const_iterator 
	       elemIt  = rowIt -> second -> getElements () . begin ();
	     elemIt   != rowIt -> second -> getElements () . end   (); ++elemIt)

	  varNumIndices. insert ((*elemIt) -> getIndex ());
      }

      // Second: check every row for elements (i,j) not in this row by
      // parallel scanning of varNumINdices

      for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::const_iterator 
	     rowIt  = (*i) -> getRows (). begin (); 
	   rowIt   != (*i) -> getRows (). end   (); ++rowIt) {

	int rowInd = rowIt -> first;

	std::set <int>::iterator vniIt = varNumIndices . begin ();

	for (std::set <CouenneScalar *, CouenneSparseVector::compare_scalars>::const_iterator 
	       elemIt  = rowIt -> second -> getElements () . begin ();
	     elemIt   != rowIt -> second -> getElements () . end   (); ++elemIt) {

	  int colInd = (*elemIt) -> getIndex ();

	  while ((vniIt != varNumIndices . end ()) && (*vniIt < colInd)) {

	    if (rowInd <= *vniIt) {
	      //printf ("missing term: %d, %d\n", rowInd, *vniIt);

	      expression *image;

	      if (rowInd == *vniIt) image = new exprPow (new exprClone (problem_ -> Var (rowInd)), new exprConst (2.));
	      else                  image = new exprMul (new exprClone (problem_ -> Var (rowInd)), new exprClone (problem_ -> Var (*vniIt)));

	      exprAux *yIJ = problem_ -> addAuxiliary (image);

	      // seek expression in the set
	      if (problem_ -> AuxSet () -> find (yIJ) == 
		  problem_ -> AuxSet () -> end ()) {

	        // no such expression found in the set, create entry therein
	        problem_ -> Variables () . push_back (yIJ);
	        problem_ -> AuxSet () -> insert (yIJ); // insert into repetition checking structure
	      }

	      (*i) -> add_element (rowInd, *vniIt, yIJ);
	      (*i) -> add_element (*vniIt, rowInd, yIJ);
	    }

	    ++vniIt;
	  }

	  if (vniIt == varNumIndices . end ())
	    break;
	  else ++vniIt;
	}
      }
    }

    // post-rescan: update 
    //
    // numbering_
    // domain_
    // commuted_
    // optimum_
    // integerRank_
    // unusedOriginalsIndices_
    //
    // since there are new variables
  }

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

    for (std::set <std::pair <int, CouenneSparseVector *>, CouenneExprMatrix::compare_pair_ind>::const_iterator 
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
      (*i) -> add_element (indexVar, size, *j); // note: problem_ -> Var (indexVar) == (*j)
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
}


/// Copy constructor
CouenneSdpCuts::CouenneSdpCuts (const CouenneSdpCuts &rhs):

  problem_          (rhs. problem_),
  doNotUse_         (rhs. doNotUse_),
  numEigVec_        (rhs. numEigVec_),
  onlyNegEV_        (rhs. onlyNegEV_),
  useSparsity_      (rhs. useSparsity_),
  fillMissingTerms_ (rhs. fillMissingTerms_) {

  for (std::vector <CouenneExprMatrix *>::const_iterator 
  	 i  = rhs.minors_ . begin ();
       i   != rhs.minors_ . end   (); ++i)

    minors_ . push_back (new CouenneExprMatrix (**i));
}


/// Assignment
CouenneSdpCuts &CouenneSdpCuts::operator= (const CouenneSdpCuts &rhs) {

  problem_          = rhs. problem_;
  doNotUse_         = rhs. doNotUse_;
  numEigVec_        = rhs. numEigVec_;
  onlyNegEV_        = rhs. onlyNegEV_;
  useSparsity_      = rhs. useSparsity_;
  fillMissingTerms_ = rhs. fillMissingTerms_;

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
sorted in non-decreasing order, hence selecting 1 will provide cuts on the most negative eigenvalue."
    );

  roptions -> AddStringOption2
    ("sdp_cuts_neg_ev",
     "Only use negative eigenvalues to create sdp cuts.",
     "yes", 
     "no", "use all eigenvalues regardless of their sign.",
     "yes", "exclude all non-negative eigenvalues."
    );

  /////////////////////////////////////////

  roptions -> AddStringOption2
    ("sdp_cuts_sparsify",
     "Make cuts sparse by greedily reducing X one column at a time before extracting eigenvectors.",
     "no", 
     "no", "",
     "yes", ""
    );

  roptions -> AddStringOption2
    ("sdp_cuts_fillmissing",
     "Create fictitious auxiliary variables to fill non-fully dense minors. Can make a difference when Q has at least one zero term.",
     "no", 
     "no", "Do not create auxiliaries and simply use Fourier-Motzkin to substitute a missing auxiliary y_ij with inequalities that use bounds and the definition y_ij = x_i x_j \
Advantage: limits the creation of auxiliaries, reformulation stays small. Default.",
     "yes", "Create (at the beginning) auxiliaries that are linearized (through McCormick) and used within an SDP cut. This allows tighter cuts although it increases the size \
of the reformulation and hence of the linear relaxation."
    );

#if 0
  roptions -> AddStringOption2
    ("sdp_cuts_",
     "",
     "yes", 
     "no", "",
     "yes", ""
    );
#endif
}
