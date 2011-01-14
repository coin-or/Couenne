/* $Id$
 *
 * Name:    CouenneAggrProbing.hpp
 * Author:  Giacomo Nannicini
 * Purpose: A bound tightener based on aggressive probing
 *
 * (C) Giacomo Nannicini, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNEAGGRPROBING_HPP
#define COUENNEAGGRPROBING_HPP

#include "BonRegisteredOptions.hpp"

#include "BonOaDecBase.hpp"
#include "CglCutGenerator.hpp"
#include "OsiColCut.hpp"
#include "OsiSolverInterface.hpp"
#include "CouenneProblem.hpp"
#include "BonCouenneSetup.hpp"

namespace Couenne {

  /// Cut Generator for aggressive BT; i.e., an aggressive probing.
  /// This probing strategy is very expensive and was initially
  /// developed to be run in parallel; hence, the user can choose to
  /// probe just a particular variable, without adding this cut
  /// generator to the list of cut generators normally employed by
  /// Couenne. However, it can also be used in the standard way; in
  /// that case, it chooses automatically the variables to probe (in a
  /// very naive way, for the moment).
  /// TODO: Implement some way to automatically choose the variables
  /// TODO: Implement the generateCuts method, for use in Branch-and-Bound

  class CouenneAggrProbing: public CglCutGenerator {

  public:

    /// Constructor
    CouenneAggrProbing(CouenneSetup* couenne,
		       const Ipopt::SmartPtr<Ipopt::OptionsList> options);

    /// Copy constructor
    CouenneAggrProbing(const CouenneAggrProbing& rhs);

    /// Destructor
    ~CouenneAggrProbing();

    /// Clone method (necessary for the abstract CglCutGenerator class)
    CouenneAggrProbing *clone () const
    {return new CouenneAggrProbing (*this);}

    /// The main CglCutGenerator; not implemented yet
    void generateCuts(const OsiSolverInterface & solver,
		      OsiCuts& cuts, 
		      const CglTreeInfo = CglTreeInfo ()) const;

    /// Probe one variable (try to tigthen the lower or the upper
    /// bound, depending on the value of the second argument), so that
    /// we can generate the corresponding column cut. This runs the
    /// main algorithm. It returns the new bound (equal to the initial
    /// one if we could not tigthen)
    double probeVariable(int index, bool probeLower);

    /// Alternative probing algorithm. This one does not work yet!
    /// Do not use, will probably segfault.
    double probeVariable2(int index, bool lower);

    /// Add list of options to be read from file
    static void registerOptions(Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

    /// Set/get maximum time to probe one variable
    void setMaxTime(double value);
    double getMaxTime() const;

    /// Set/get maximum number of failed steps
    void setMaxFailedSteps(int value);
    int getMaxFailedSteps() const;

    /// Set/get maximum number of nodes to probe one variable
    void setMaxNodes(int value);
    int getMaxNodes() const;

    /// Set/get restoreCutoff parameter (should we restore the initial
    /// cutoff value after each probing run?)
    void setRestoreCutoff(bool value);
    bool getRestoreCutoff() const;

  protected:    

    /// Pointer to the CouenneProblem representation
    CouenneSetup* couenne_;

    /// Number of columns (want to have this handy)
    int numCols_;

    /// Maximum time to probe one variable
    double maxTime_;

    /// Maximum number of failed iterations
    int maxFailedSteps_;

    /// Maximum number of nodes in probing
    int maxNodes_;

    /// Restore initial cutoff (value and solution)?
    bool restoreCutoff_;

    /// Initial cutoff
    double initCutoff_;

  };
}

#endif
