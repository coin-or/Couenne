/* $Id$
 *
 * Name:    CouenneCrossConv.hpp
 * Author:  Pietro Belotti
 * Purpose: Convexification cuts on redundant relationships between auxiliaries
 *
 * (C) Pietro Belotti, 2010-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNECROSSCONV_HPP
#define COUENNECROSSCONV_HPP

#include "BonRegisteredOptions.hpp"

#include "CglConfig.h"
#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
#include "CouenneJournalist.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
}

namespace Couenne {

  class CouenneProblem;

  /// Base class definition for relations between auxiliaries

  class AuxRelation {

  public:

    virtual int findRelations () = 0;

    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;
  protected:

  };
  
  /// Identifies 5-ples of variables of the form
  ///
  /// x_3 := log x_1
  /// x_4 := log x_2
  /// x_5 := x_1 x_2 in [l,u]
  ///
  /// and generates a cut 
  ///
  /// x_3 + x_4 in [max {0, log l}, max {0, log u}].
  ///
  /// This has to be repeatedly generated, even when l=u (l and/or u
  /// could change in other nodes).

  class SumLogAuxRel: public AuxRelation {

  public:

    virtual int findRelations ();

    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;
  };


  /// Identifies 5-ples of variables of the form
  ///
  /// x_k := x_i x_j
  /// x_l := x_i x_p
  ///
  /// x_q := x_k x_p     OR     x_q := x_k / x_j
  /// x_r := x_k x_j            x_r := x_l / x_p
  ///
  /// and generates, ONLY ONCE, a cut 
  ///
  /// x_q = x_r (in both cases).

  class MultiProdRel: public AuxRelation {

  public:

    virtual int findRelations ();

    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;
  };

  /// Identifies 5-tuple of the form
  /// 
  /// x_j := x_i / x_k
  /// x_p := x_i / x_q
  ///
  /// x_l := x_j / x_p   OR   x_l := x_j x_k
  /// x_m := x_q / x_k        x_m := x_p x_q
  ///
  /// and generates, ONLY once, a cut
  ///
  /// x_l = x_m (in both cases).

  class BiProdDivRel: public AuxRelation {

  public:

    virtual int findRelations ();

    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;
  };

  /// Identifies 5-tuple of the form
  /// 
  /// x_j := x_i ^ alpha
  /// x_p := x_i ^ beta
  ///
  /// and generates cuts based on the relation
  ///
  /// x_p = x_j ^ {beta/alpha}

  class PowRel: public AuxRelation {

  public:

    virtual int findRelations ();

    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;
  };


  /// Cut Generator that uses relationships between auxiliaries

  class CouenneCrossConv: public CglCutGenerator {

  public:

    /// constructor
    CouenneCrossConv (CouenneProblem *,
		      JnlstPtr,
		      const Ipopt::SmartPtr <Ipopt::OptionsList>);

    /// copy constructor
    CouenneCrossConv  (const CouenneCrossConv &);

    /// destructor
    virtual ~CouenneCrossConv ();

    /// clone method (necessary for the abstract CglCutGenerator class)
    virtual CouenneCrossConv *clone () const
    {return new CouenneCrossConv (*this);}

    /// the main CglCutGenerator
    virtual void generateCuts (const OsiSolverInterface &, 
		       OsiCuts &, 
		       const CglTreeInfo = CglTreeInfo ())
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
    const
#endif
    ;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

    /// Set up data structure to detect redundancies
    virtual void setup ();

  protected:

    /// Journalist
    JnlstPtr jnlst_;

    /// pointer to the CouenneProblem representation
    CouenneProblem *problem_;
  };
}

#endif
