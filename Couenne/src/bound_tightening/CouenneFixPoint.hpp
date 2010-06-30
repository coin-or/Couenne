/* $Id$
 *
 * Name:    CouenneFixPoint.hpp
 * Author:  Pietro Belotti
 * Purpose: A bound tightener based on fixpoint computation
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEFIXPOINT_HPP
#define COUENNEFIXPOINT_HPP

#include "BonRegisteredOptions.hpp"

#include "BonOaDecBase.hpp"
#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
#include "OsiSolverInterface.hpp"

namespace Couenne {

  /// Cut Generator for FBBT fixpoint

  class CouenneFixPoint: public CglCutGenerator {

  public:

    /// constructor
    CouenneFixPoint (const Ipopt::SmartPtr<Ipopt::OptionsList> );

    /// copy constructor
    CouenneFixPoint  (const CouenneFixPoint &);

    /// destructor
    ~CouenneFixPoint ();

    /// clone method (necessary for the abstract CglCutGenerator class)
    CouenneFixPoint *clone () const
    {return new CouenneFixPoint (*this);}

    /// the main CglCutGenerator
    void generateCuts (const OsiSolverInterface &, 
		       OsiCuts &, 
		       const CglTreeInfo = CglTreeInfo ()) const;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  private:

    /// should we use an extended model or a more compact one?
    bool extendedModel_;
  };
}

#endif
