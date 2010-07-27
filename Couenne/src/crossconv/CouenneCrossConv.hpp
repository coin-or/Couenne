/* $Id$
 *
 * Name:    CouenneCrossConv.hpp
 * Author:  Pietro Belotti
 * Purpose: Convexification cuts on redundant relationships between auxiliaries
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECROSSCONV_HPP
#define COUENNECROSSCONV_HPP

#include "BonRegisteredOptions.hpp"

//#include "BonOaDecBase.hpp"
#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
//#include "OsiSolverInterface.hpp"
#include "CouenneJournalist.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
}

namespace Couenne {

  class CouenneProblem;

  /// Cut Generator for FBBT fixpoint

  class CouenneCrossConv: public CglCutGenerator {

  public:

    /// constructor
    CouenneCrossConv (CouenneProblem *,
		      JnlstPtr,
		      const Ipopt::SmartPtr <Ipopt::OptionsList>);

    /// copy constructor
    CouenneCrossConv  (const CouenneCrossConv &);

    /// destructor
    ~CouenneCrossConv ();

    /// clone method (necessary for the abstract CglCutGenerator class)
    CouenneCrossConv *clone () const
    {return new CouenneCrossConv (*this);}

    /// the main CglCutGenerator
    void generateCuts (const OsiSolverInterface &, 
		       OsiCuts &, 
		       const CglTreeInfo = CglTreeInfo ()) const;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

    /// Set up data structure to detect redundancies
    void setup (CouenneProblem *p);

  protected:

    /// Journalist
    JnlstPtr jnlst_;

    /// pointer to the CouenneProblem representation
    CouenneProblem *problem_;
  };
}

#endif
