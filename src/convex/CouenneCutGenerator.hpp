/*
 *
 * Name:    CouenneCutGenerator.hpp
 * Author:  Pietro Belotti
 * Purpose: a convexification cut generator for MINLP problems
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_CUT_GENERATOR_HPP
#define COUENNE_CUT_GENERATOR_HPP

//#include "BonRegisteredOptions.hpp"

#include "BonOaDecBase.hpp"
#include "CglConfig.h"
#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
#include "BonAuxInfos.hpp"
#include "BonBabInfos.hpp"
#include "OsiSolverInterface.hpp"
#include "CouenneConfig.h"
#include "CouenneJournalist.hpp"
#include "CouenneTypes.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
  class Journalist;
}

namespace Bonmin {
  class RegisteredOptions;
  class BabInfo;
  class OsiTMINLPInterface;
  class BabSetupBase;
}

struct ASL;

namespace Couenne {

class CouenneProblem;
class funtriplet;

/// Cut Generator for linear convexifications

class COUENNELIB_EXPORT CouenneCutGenerator: public CglCutGenerator {

 protected:

  /// True if no convexification cuts have been generated yet for this
  /// problem
  mutable bool firstcall_;

  /// True if we should add the violated cuts only, false if all of
  /// them should be added
  mutable bool addviolated_;

  /// what kind of sampling should be performed?
  enum conv_type convtype_;

  /// how many cuts should be added for each function?
  int nSamples_;

  /// pointer to symbolic repr. of constraint, variables, and bounds
  CouenneProblem *problem_;

  /// number of cuts generated at the first call
  mutable int nrootcuts_;

  /// total number of cuts generated
  mutable int ntotalcuts_;

  /// separation time (includes generation of problem)
  mutable double septime_;

  /// Record obj value at final point of CouenneConv.
  mutable double objValue_;

  /// nonlinear solver interface as used within Bonmin (used at first
  /// Couenne pass of each b&b node
  Bonmin::OsiTMINLPInterface *nlp_;

  /// pointer to the Bab object (used to retrieve the current primal
  /// bound through bestObj())
  Bonmin::Bab *BabPtr_;

  /// signal infeasibility of current node (found through bound tightening)
  mutable bool infeasNode_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// Time spent at the root node
  mutable double rootTime_;

  /// Check all generated LPs through an independent call to
  /// OsiClpSolverInterface::initialSolve()
  bool check_lp_;

  /// Take advantage of OsiClpSolverInterface::tightenBounds (), known
  /// to have caused some problems some time ago
  bool enable_lp_implied_bounds_;

  /// Running count of printed info lines
  mutable int lastPrintLine;

 public:

  /// constructor
  CouenneCutGenerator  (Bonmin::OsiTMINLPInterface * = NULL,
			Bonmin::BabSetupBase *base = NULL,
			CouenneProblem * = NULL,
			struct ASL * = NULL);

  /// copy constructor
  CouenneCutGenerator  (const CouenneCutGenerator &);

  /// destructor
  ~CouenneCutGenerator ();

  /// clone method (necessary for the abstract CglCutGenerator class)
  CouenneCutGenerator *clone () const
  {return new CouenneCutGenerator (*this);}

  /// return pointer to symbolic problem
  inline CouenneProblem *Problem () const
  {return problem_;}

  /// return pointer to symbolic problem
  inline void setProblem (CouenneProblem *p)
  {problem_ = p;}

  /// total number of variables (original + auxiliary)
  int getnvars () const;

  /// has generateCuts been called yet?
  inline bool isFirst () const
  {return firstcall_;}

  /// should we add the violated cuts only (true), or all of them (false)?
  inline bool addViolated () const
  {return addviolated_;}

  /// get convexification type (see CouenneTypes.h)
  inline enum conv_type ConvType () const
  {return convtype_;}

  /// get number of convexification samples
  inline int nSamples () const
  {return nSamples_;}

  /// the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &,
		     OsiCuts &,
		     const CglTreeInfo = CglTreeInfo ())
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
   const
#endif
  ;

  /// create cut and check violation. Insert and return status
  int createCut (OsiCuts &, // cutset to insert
		 CouNumber, // lb
		 CouNumber, // ub
		            // index, coeff  (index -1: "don't care")
		 int,    CouNumber,    // of first  term
		 int=-1, CouNumber=0., // of second term
		 int=-1, CouNumber=0., // of third  term
		 bool = false) const;  // is it a global cut? No, by default

  /// create cut and check violation. Other version with only one bound
  int createCut (OsiCuts &, // cutset to insert
		 CouNumber, // rhs
		 int,       // sign: -1: <=, 0: =, +1: >=
		            // index, coeff  (index -1: "don't care")
		 int,    CouNumber,    // of first  term
		 int=-1, CouNumber=0., // of second term
		 int=-1, CouNumber=0., // of third  term
		 bool = false) const;  // is it a global cut? No, by default

  /// Add general linear envelope to convex function, given its
  /// variables' indices, the (univariate) function and its first
  /// derivative
  void addEnvelope (OsiCuts &,
		    int,
		    unary_function, unary_function,
		    int, int,
		    CouNumber, CouNumber, CouNumber,
		    t_chg_bounds * = NULL,
		    bool = false) const;

  /// Add general linear envelope to convex function, given its
  /// variables' indices, the (univariate) function and its first
  /// derivative
  void addEnvelope (OsiCuts &,
		    int,
		    funtriplet *,
		    int, int,
		    CouNumber, CouNumber, CouNumber,
		    t_chg_bounds * = NULL,
		    bool = false) const;

  /// Add half-plane through (x1,y1) and (x2,y2) -- resp. 4th, 5th,
  /// 6th, and 7th argument
  int addSegment (OsiCuts &, int, int,
		  CouNumber, CouNumber,
		  CouNumber, CouNumber, int) const;

  /// add tangent at given poing (x,w) with given slope
  int addTangent (OsiCuts &, int, int,
		  CouNumber, CouNumber,
		  CouNumber, int) const;

  /// Method to set the Bab pointer
  void setBabPtr (Bonmin::Bab *p)
  {BabPtr_ = p;}

  /// Get statistics
  void getStats (int &nrc, int &ntc, double &st) {
    nrc = nrootcuts_;
    ntc = ntotalcuts_;
    st  = septime_;
  }

  /// Allow to get and set the infeasNode_ flag (used only in generateCuts())
  bool &infeasNode () const
  {return infeasNode_;}

  /// generate OsiRowCuts for current convexification
  void genRowCuts (const OsiSolverInterface &, OsiCuts &cs,
		   int, int *, t_chg_bounds * = NULL) const;

  /// generate OsiColCuts for improved (implied and propagated) bounds
  void genColCuts (const OsiSolverInterface &, OsiCuts &, int, int *) const;

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  /// print node, depth, LB/UB/LP info
  void printLineInfo() const;

  /// Provide Journalist
  inline ConstJnlstPtr Jnlst() const
  {return ConstPtr (jnlst_);}

  void setJnlst (JnlstPtr jnlst__)
  { jnlst_ = jnlst__; }

  /// Time spent at root node
  double &rootTime ()
  {return rootTime_;}

  /// return check_lp flag (used in CouenneSolverInterface)
  bool check_lp () const
  {return check_lp_;}

  /// returns value of enable_lp_implied_bounds_
  bool enableLpImpliedBounds () const
  {return enable_lp_implied_bounds_;}
};


/// translate sparse to dense vector (should be replaced)
COUENNELIB_EXPORT
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);

}

#endif
