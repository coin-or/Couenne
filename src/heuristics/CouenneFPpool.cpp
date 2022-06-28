/*
 *
 * Name:    CouenneFPpool.cpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Pool of MILP- (and why not? NLP-) feasible solutions for
 *          restart use in the Feasibility Pump
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneFPpool.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprAux.hpp"

static const double COMP_TOLERANCE = COUENNE_EPS; // not global

using namespace Couenne;

/// CouenneProblem-aware constructor
CouenneFPsolution::CouenneFPsolution (CouenneProblem *p, CouNumber *x, bool copied):

  x_        (NULL),
  n_        (p -> nVars ()),
  nNLinf_   (0),
  nIinf_    (0),
  objVal_   (0.),
  maxNLinf_ (0.),
  maxIinf_  (0.),
  copied_   (copied),
  problem_  (p) {

  /// NOTE: the copied flag means we won't use this solution for
  /// anything but testing it against a tabu list, so we don't need to
  /// perform all checks

  if (copied_) {
    x_ = x;
    return;
  }

  x_ = CoinCopyOfArray (x, p -> nVars ());

  for (std::vector <exprVar *>::iterator i = p -> Variables (). begin (); 
       i != p -> Variables (). end ();
       ++i) {

    if ((*i) -> Multiplicity () <= 0)
      continue;

    CouNumber 
      vval = (**i) ();

    if ((*i) -> isInteger ()) {

      double inf = CoinMax (vval - floor (vval + COUENNE_EPS),
			    ceil (vval - COUENNE_EPS) - vval);

      if (inf > COUENNE_EPS) {

	++nIinf_;
 
	if (inf > maxIinf_) 
	  maxIinf_ = inf;
      }
    }

    if (((*i) -> Type () == AUX) &&
	((*i) -> Image () -> Linearity () > LINEAR)) {

      double
	diff = 0.,
	fval = (*((*i) -> Image ())) ();

      if      ((*i) -> sign () != expression::AUX_GEQ) diff = CoinMax (diff, vval - fval);
      else if ((*i) -> sign () != expression::AUX_LEQ) diff = CoinMax (diff, fval - vval);

      if (diff > COUENNE_EPS) {

	++nNLinf_;

	if (diff > maxNLinf_)
	  maxNLinf_ = diff;
      }
    }
  }
}


/// copy constructor
CouenneFPsolution::CouenneFPsolution (const CouenneFPsolution &src):
  x_          (src.x_ ? CoinCopyOfArray (src.x_, src.n_) : NULL),
  n_          (src.n_),
  nNLinf_     (src.nNLinf_),
  nIinf_      (src.nIinf_),
  objVal_     (src.objVal_),
  maxNLinf_   (src.maxNLinf_),
  maxIinf_    (src.maxIinf_),
  copied_     (false),
  problem_    (src.problem_) {}


/// assignment
CouenneFPsolution &CouenneFPsolution::operator= (const CouenneFPsolution &src) {

  x_         = src.x_ ? CoinCopyOfArray (src.x_, src.n_) : NULL;
  n_         = src.n_;
  nNLinf_    = src.nNLinf_;
  nIinf_     = src.nIinf_;
  objVal_    = src.objVal_;
  maxNLinf_  = src.maxNLinf_;
  maxIinf_   = src.maxIinf_;
  copied_    = false;
  problem_   = src.problem_;

  return *this;
}

/// destructor
CouenneFPsolution::~CouenneFPsolution () {

  if (x_ && !copied_)
    delete [] x_;
}

/// basic comparison procedure -- what to compare depends on user's choice
bool CouenneFPsolution::compare (const CouenneFPsolution &other, enum what_to_compare comparedTerm) const {

  switch (comparedTerm) {

  case SUM_INF:      if (maxNLinf_ + maxIinf_ < other.maxNLinf_ + other.maxIinf_)                                             return true;
  case OBJVAL:       if (objVal_              < other.objVal_ - COUENNE_EPS * CoinMax (1., CoinMax (objVal_, other.objVal_))) return true;
  case SUM_NINF:     if (nNLinf_   + nIinf_   < other.nNLinf_   + other.nIinf_)                                               return true;
    // so that if objective is close these two solutions will be deemed equal

  case INTEGER_VARS: {

    // lexicographical comparison: unless the two solutions have the
    // same integer subvector, comparison will tell them apart

    for (std::vector <exprVar *>::iterator i = problem_ -> Variables (). begin (); 
	 i != problem_ -> Variables (). end ();
	 ++i)

      if (((*i) -> Multiplicity () > 0) && 
	  (*i) -> isInteger ()) {

	int indVar = (*i) -> Index ();

	// LEXICOGRAPHICAL ORDERING: return true upon finding the
	// first element with lower element

#if 0
	if (0) {
	  printf (" (%d,%g,%g,%g)", indVar, x_ [indVar], other.x_ [indVar], x_ [indVar] - other.x_ [indVar]);

	  if (x_ [indVar] < other.x_ [indVar] - COUENNE_EPS)
	    printf ("\n-----\n");
	}
#endif
	if (x_ [indVar] < other.x_ [indVar] - COUENNE_EPS)
	  return true;
      }

    //printf ("\n####\n");
    return false;
  }

  case ALL_VARS: {
    // lexicographical comparison: unless the two solutions have the
    // same subvector, comparison will tell them apart

    for (std::vector <exprVar *>::iterator i = problem_ -> Variables (). begin (); 
	 i != problem_ -> Variables (). end ();
	 ++i) 

      if ((*i) -> Multiplicity () > 0) {

	int indVar = (*i) -> Index ();

	if (x_ [indVar] < other.x_ [indVar] + COUENNE_EPS)
	  return true;
      }

    return false;
  }
  }

  printf ("CouenneFPsolution::compare: bad compared term\n");
  abort ();
}


/// copy constructor
CouenneFPpool::CouenneFPpool (const CouenneFPpool &src):
  set_     (src.set_),
  problem_ (src.problem_) {}


/// assignment
CouenneFPpool &CouenneFPpool::operator= (const CouenneFPpool &src) {

  set_     = src.set_;
  problem_ = src.problem_;

  return *this;
}

/// compare, base version
bool compareSol::operator() (const CouenneFPsolution &one, 
			     const CouenneFPsolution &two) const {

  return one. compare (two, comparedTerm_);

  // const double
  //   *x1 = one.x (),
  //   *x2 = two.x ();

  // int n = one.n ();

  // while (n--)
  //   if ((*x1++ - *x2++) <= COMP_TOLERANCE)
  //     return true;

  // return false;
}

/// finds, in pool, solution x closest to nSol; removes it from the
/// pool and overwrites it to sol
void CouenneFPpool::findClosestAndReplace (double *&sol, const double *nSol, int nvars)  {

   double bestdist = COIN_DBL_MAX;
   std::set <CouenneFPsolution, compareSol>::iterator bestsol = set_. end ();

   if( nSol )
   {
     for (std::set <CouenneFPsolution, compareSol>::iterator i = set_. begin (); 
           i != set_. end (); ++i)
      {
	//compute distance of pool solution and NLP solution

	double 
	   dist = 0.0,
	   delta;

	const double 
	  *x = i -> x (),
	  *s = nSol;

	bool move_on = false;

	for (int j = nvars, k=0; j--; ++k) {

	  delta = *x++ - *s++;

	  /// forget about this variable if eliminated by
	  /// reformulation
	  if (problem_ -> Var (k) -> Multiplicity () <= 0)
	    continue;

	  dist += delta * delta;

	  if (dist >= bestdist) { // interrupt check of this solution
				  // if already above current best
	    move_on = true;
	    break;
	  }
	}

	if (move_on) 
	  continue;

	//update best solution
	if( dist < bestdist )
         {
            bestdist = dist;
            bestsol = i;
         }
      }
   }
   else 
      bestsol = set_. begin ();

   if( bestsol != set_. end () )
   {
     delete [] sol;
     sol = CoinCopyOfArray ((*bestsol).x(), nvars); 
     set_. erase(bestsol);
   }   
}
