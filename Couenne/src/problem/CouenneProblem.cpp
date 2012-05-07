/* $Id$
 *
 * Name:    CouenneProblem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CbcBranchActual.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprIVar.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprOpp.hpp"

#include "CouenneProblem.hpp"
#include "CouenneGlobalCutOff.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneLQelems.hpp"

#include "CouenneJournalist.hpp"
#include "CouenneObject.hpp"

using namespace Couenne;

/// methods to add objective function. 

void CouenneProblem::addObjective (expression *newobj, const std::string &sense) {
  objectives_ . push_back
    (new CouenneObjective ((sense == "min") ? 
			   newobj : new exprOpp (new exprClone (newobj))));
}


/// methods to add nonlinear constraints:

/// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs) {

  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

/// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (COUENNE_INFINITY)));
}

/// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (-COUENNE_INFINITY), rhs));
}

/// Add (non linear) objective function
void CouenneProblem::setObjective (int indObj, expression * newObj, const std::string &sense) {
  objectives_ [indObj] = (new CouenneObjective ((sense == "min") ? 
						newObj : new exprOpp (new exprClone (newObj))));
}


/// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb, expression *ub) {
  if (!lb) lb = new exprConst (0.);
  if (!ub) ub = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}



/// add variable to the problem -- check whether it is integer (isDiscrete)

expression *CouenneProblem::addVariable (bool isDiscrete, Domain *d) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size (), d)) :
    (new exprVar  (variables_ . size (), d));

  variables_ . push_back (var);

  if (isDiscrete) 
    nIntVars_++;

  nOrigVars_++;

  return var;
}


/// add auxiliary variable and associate it with pointer to expression
/// given as argument
exprAux *CouenneProblem::addAuxiliary (expression *symbolic) {

  // check if image is already in the expression database auxSet_
  std::set <exprAux *, compExpr>::iterator i;

  int var_ind = variables_ . size ();
  domain_. current () -> resize (var_ind + 1);

  symbolic -> getBounds (domain_. lb (var_ind),
			 domain_. ub (var_ind));

  // create new aux associated with that expression
  exprAux *w = new exprAux (symbolic,
			    var_ind,
			    1 + symbolic -> rank (),
			    exprAux::Unset, 
			    &domain_);
  //symbolic -> isInteger () ? exprAux::Integer : exprAux::Continuous);

  //  w -> linkDomain (&domain_);

  // seek expression in the set
  if ((i = auxSet_ -> find (w)) == auxSet_ -> end ()) {

    // no such expression found in the set, create entry therein
    variables_ . push_back (w);
    auxSet_ -> insert (w); // insert into repetition checking structure
    graph_  -> insert (w); // insert into acyclic structure

  } else {  // otherwise, just return the entry's pointer

    w -> Image (NULL); // otherwise "delete w" will also delete user given expression "symbolic"
    delete w;
    w = *i;
    (*i) -> increaseMult ();
  }

  return w;
}


/// translates pair (indices, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexL, 
				    CouNumber *coeff,
				    std::vector <std::pair <exprVar *, CouNumber> > &lcoeff) {
  // TODO: sort

  for (int i=0; indexL [i] >= 0; i++)
    lcoeff.push_back (std::pair <exprVar *, CouNumber> (Var (indexL [i]), coeff [i]));
}


/// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexI,
				    int *indexJ,
				    CouNumber *coeff,
				    std::vector <quadElem> &qcoeff) {
  // TODO: sort

  for (int i=0; indexI [i] >= 0; i++)
    qcoeff.push_back (quadElem (Var (indexI [i]), Var (indexJ [i]), coeff [i]));
}


/// fill in the integerRank_ array
void CouenneProblem::fillIntegerRank () const {

  if (integerRank_)
    return;

  int nvars = nVars ();

  integerRank_ = new int [nvars];

  // 0: fractional
  // 1: integer
  // k: integer,    depending on at least one integer with associated value k-1, or
  // k: fractional, depending on at least one integer with associated value k

  for (int ii = 0; ii < nvars; ii++) {

    int index = numbering_ [ii];

    if (Var (index) -> Multiplicity () <= 0) {
      integerRank_ [index] = 0;
      continue;
    }

    bool isInt = Var (index) -> isDefinedInteger ();

    // sets all originals to either 0 (fractional) or 1 (integer)
    integerRank_ [index] = (isInt) ? 1 : 0;

    if (Var (index) -> Type () == AUX) {

      std::set <int> deplist;

      if (Var (index) -> Image () -> DepList (deplist, STOP_AT_AUX) != 0) // depends on something
	for (std::set <int>::iterator i = deplist.begin (); i != deplist.end (); ++i) {

	  int token = integerRank_ [*i];
	  if (isInt) token++;

	  if (token > integerRank_ [index]) // there's a free integer below us
	    integerRank_ [index] = token;
	}
    }
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "Free (original) integers\n");
  for (int i=0; i<nOrigVars_; i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, integerRank_ [i]);

  // fill in numberInRank_
  for (int i=0; i<nOrigVars_ - ndefined_; i++)
    if ((variables_ [i] -> isDefinedInteger ()) &&
	(variables_ [i] -> Multiplicity () > 0)) {

    int rank = integerRank_ [i];

    if (numberInRank_.size () <= (unsigned int) rank)
      for (int j=numberInRank_.size (); j <= rank; j++)
	numberInRank_ .push_back (0);

    numberInRank_ [rank] ++;
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "numInteger [neglect non-originals]\n");
  for (unsigned int i=0; i<numberInRank_.size(); i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, numberInRank_ [i]);
}


/// Get cutoff
CouNumber CouenneProblem::getCutOff () const
{return pcutoff_ -> getCutOff ();}

/// Get cutoff solution
CouNumber *CouenneProblem::getCutOffSol () const
{return pcutoff_ -> getCutOffSol ();}

/// Provide Journalist
ConstJnlstPtr CouenneProblem::Jnlst () const 
{return ConstPtr (jnlst_);}

// set lastPrioSort_
void CouenneProblem::setLastPrioSort(int givenLastPS) {
  lastPrioSort_ = givenLastPS;
}
