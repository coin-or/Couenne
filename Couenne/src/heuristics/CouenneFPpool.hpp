/* $Id$
 *
 * Name:    CouenneFPpool.hpp
 * Authors: Pietro Belotti
 *          Timo Berthold, ZIB Berlin
 * Purpose: Pool of MILP- (and why not? NLP-) feasible solutions for
 *          restart use in the Feasibility Pump
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneFeasPump.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneTNLP.hpp"

namespace Couenne {

  enum what_to_compare {SUM_NINF, SUM_INF, OBJVAL};

  /// what term to compare: the sum of infeasibilities, the sum of
  /// numbers of infeasible terms, or the objective function
  static enum what_to_compare comparedTerm_;

  /// Class containing a solution with infeasibility evaluation
  class CouenneFPsolution {

  protected:

    CouNumber *x_;        ///< solution
    int        n_;        ///< number of variables (for independence from CouenneProblem
    int        nNLPinf_;  ///< number of NLP     infeasibilities
    int        nIinf_;    ///< number of integer infeasibilities
    CouNumber  objVal_;   ///< objective function value
    CouNumber  maxNLinf_; ///< maximum NLP     infeasibility
    CouNumber  maxIinf_;  ///< maximum integer infeasibility

  public:

    CouenneFPsolution (CouenneProblem *p, CouNumber *x); ///< CouenneProblem-aware constructor

    /// independent constructor --- must provide other data as no
    /// CouenneProblem to compute them
    CouenneFPsolution (CouNumber *x, 
		       int n, 
		       int nNLPinf        = 1, 
		       int nIinf          = 1, 
		       CouNumber objVal   = COIN_DBL_MAX,
		       CouNumber maxNLinf = COIN_DBL_MAX,
		       CouNumber maxIinf  = COIN_DBL_MAX);

    /// basic comparison procedure -- what to compare depends on user's choice
    int compare (const CouenneFPsolution &other, enum what_to_compare comparedTerm) const;
  };

  /// compare, base version
  inline int operator< (const CouenneFPsolution &one, 
			const CouenneFPsolution &two)
  {return one.compare (two, comparedTerm_);}

  /// Pool of solutions
  class CouenneFPpool {

  protected:

    /// Pool
    std::priority_queue <CouenneFPsolution> queue_;

  public:

    /// simple constructor (empty pool)
    CouenneFPpool (enum what_to_compare c)
    {comparedTerm_ = c;}

    /// return the main object in this class
    std::priority_queue <CouenneFPsolution> &Queue ()
    {return queue_;}
  };
}
