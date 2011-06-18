/* $Id$
 *
 * Name:    CouenneFPpool.hpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Pool of MILP- (and why not? NLP-) feasible solutions for
 *          restart use in the Feasibility Pump and sets of solutions
 *          to be used as tabu list
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneFPpool_hpp
#define CouenneFPpool_hpp

#include <set>

#include "CouenneTypes.hpp"
#include "CoinFinite.hpp"

namespace Couenne {

  class CouenneProblem;
  class CouenneFeasPump;

  /// what term to compare: the sum of infeasibilities, the sum of
  /// numbers of infeasible terms, or the objective function
  static enum what_to_compare {SUM_NINF = 0, SUM_INF, OBJVAL} comparedTerm_;

  /// Class containing a solution with infeasibility evaluation
  class CouenneFPsolution {

  protected:

    CouNumber *x_;        ///< solution
    int        n_;        ///< number of variables (for independence from CouenneProblem
    int        nNLinf_;   ///< number of NL      infeasibilities
    int        nIinf_;    ///< number of integer infeasibilities
    CouNumber  objVal_;   ///< objective function value
    CouNumber  maxNLinf_; ///< maximum NL      infeasibility
    CouNumber  maxIinf_;  ///< maximum integer infeasibility

    /// This is a temporary copy, not really a solution holder. As a
    /// result, all the above members are meaningless for copied
    /// solutions

    bool       copied_;   

  public:

    CouenneFPsolution (CouenneProblem *p, CouNumber *x, bool copied = false); ///< CouenneProblem-aware constructor

    /// independent constructor --- must provide other data as no
    /// CouenneProblem to compute them
    CouenneFPsolution (CouNumber *x,
		       int n,
		       int nNLinf         = 1,
		       int nIinf          = 1,
		       CouNumber objVal   = COIN_DBL_MAX,
		       CouNumber maxNLinf = COIN_DBL_MAX,
		       CouNumber maxIinf  = COIN_DBL_MAX);

    CouenneFPsolution (const CouenneFPsolution &src); ///< copy constructor

    CouenneFPsolution &operator= (const CouenneFPsolution &src); ///< assignment

    ~CouenneFPsolution (); ///< destructor

    /// returns size
    const int n () const {return n_;}

    /// returns vector
    const double *x () const {return x_;}

    /// basic comparison procedure -- what to compare depends on user's choice
    bool compare (const CouenneFPsolution &other, enum what_to_compare comparedTerm) const;
  };


  /// compare, base version
  inline bool operator< (const CouenneFPsolution &one, 
			 const CouenneFPsolution &two)
  {return one.compare (two, comparedTerm_);}


  /// class for comparing solutions (used in tabu list)
  class compareSol {

  public:
    bool operator () (const CouenneFPsolution &one, 
		      const CouenneFPsolution &two) const;
  };


  /// Pool of solutions
  class CouenneFPpool {

  protected:

    /// Pool
    std::set <CouenneFPsolution> set_;

  public:

    /// simple constructor (empty pool)
    CouenneFPpool (enum what_to_compare c)
    {comparedTerm_ = c;}

    /// copy constructor
    CouenneFPpool (const CouenneFPpool &src);

    /// assignment
    CouenneFPpool &operator= (const CouenneFPpool &src);

    /// return the main object in this class
    std::set <CouenneFPsolution> &Set ()
    {return set_;}

    /// finds, in pool, solution x closest to sol; removes it from the
    /// pool and overwrites it to sol
    void findClosestAndReplace (double *sol, double *nSol, int nvars) ;
  };
}

#endif
