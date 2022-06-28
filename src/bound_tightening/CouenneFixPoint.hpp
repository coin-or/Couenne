/*
 *
 * Name:    CouenneFixPoint.hpp
 * Author:  Pietro Belotti
 * Purpose: A bound tightener based on fixpoint computation
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEFIXPOINT_HPP
#define COUENNEFIXPOINT_HPP

#include "BonRegisteredOptions.hpp"

#include "BonOaDecBase.hpp"
#include "CglConfig.h"
#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
#include "OsiSolverInterface.hpp"

#include "CouenneBTPerfIndicator.hpp"

namespace Couenne {

  class CouenneProblem;

  /// Cut Generator for FBBT fixpoint

  class COUENNELIB_EXPORT CouenneFixPoint: public CglCutGenerator {

  public:

    /// constructor
    CouenneFixPoint (CouenneProblem *,
		     const Ipopt::SmartPtr<Ipopt::OptionsList>);

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
		       const CglTreeInfo = CglTreeInfo ())
#if CGL_VERSION_MAJOR == 0 && CGL_VERSION_MINOR <= 57
    const
#endif
    ;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  protected:

    /// should we use an extended model or a more compact one?
    bool extendedModel_;

    /// pointer to the CouenneProblem representation
    CouenneProblem *problem_;

    /// Is this the first call?
    mutable bool firstCall_;

    /// CPU time
    mutable double CPUtime_;

    /// Number of bounds tightened
    mutable int nTightened_;

    /// Level at which to stop
    int levelStop_;

    /// Create a single cut
    void createRow (int, int,
		    int,
		    OsiSolverInterface *,
		    const int    *,
		    const double *,
		    const double,
		    const int,
		    bool,
		    int, int) const;

    /// Performance indicator
    CouenneBTPerfIndicator perfIndicator_;
  };
}

#endif
