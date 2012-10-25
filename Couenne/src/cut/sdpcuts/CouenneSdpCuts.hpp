/* $Id$
 *
 * Name:    CouenneSdpCuts.hpp
 * Author:  Pietro Belotti
 * Purpose: wrapper for Couenne to insert sdpcuts
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef CouenneSdpCuts_hpp
#define CouenneSdpCuts_hpp

#include "BonRegisteredOptions.hpp"

#include "CglCutGenerator.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
}

namespace Couenne {

  class CouenneProblem;
  class CouenneMatrix;

  ///
  /// These are cuts of the form
  ///
  /// a' X a >= 0
  ///
  /// where X = / 1 x' \   <--- where x' stands for "x transposed"
  ///           \ x X0 /
  ///
  /// and X0 = (x_ij)_{i,j in N}, and x_ij is the auxiliary variable
  /// for x_i * x_j. After reformulation, matrices like X0 arise
  /// naturally and can be used to separate cuts that help strengthen
  /// the lower bound. See Sherali and Fraticelli for the base idea,
  /// and Qualizza, Belotti and Margot for an efficient rework and its
  /// implementation. Andrea Qualizza's code has been made open source
  /// and is used here (thanks Andrea!).
  ///

  class CouenneSdpCuts: public CglCutGenerator {

  protected:

    bool doNotUse_; ///< after construction, true if there are enough
		    ///< product terms to justify application. If not,
		    ///< do not add this cut generator

    std::vector <CouenneMatrix *> minors_; ///< minors on which to apply cuts

  public:

    CouenneSdpCuts  (CouenneProblem *);                  ///< Constructor
    ~CouenneSdpCuts ();                                  ///< Destructor
    CouenneSdpCuts  &operator= (const CouenneSdpCuts &); ///< Assignment
    CouenneSdpCuts             (const CouenneSdpCuts &); ///< Copy constructor
    virtual CglCutGenerator *clone () const;             ///< Cloning constructor

    const bool doNotUse () const {return doNotUse_;}

    /// The main CglCutGenerator
    virtual void generateCuts (const OsiSolverInterface &, 
			       OsiCuts &, 
			       const CglTreeInfo = CglTreeInfo ()) const;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);
  };
}

#endif
