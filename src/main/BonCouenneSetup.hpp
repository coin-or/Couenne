/* */
// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007

#ifndef BonCouenneSetup_H
#define BonCouenneSetup_H

#include "BonBabSetupBase.hpp"
#include "BonBonminSetup.hpp"
#include "CbcFeasibilityBase.hpp"

#include "CouenneConfig.h"
#include "CouenneTNLP.hpp"

struct ASL;

namespace Bonmin{
  class Bab;
}

namespace Couenne {

  class CouenneInterface;
  class CouenneCutGenerator;
  class CouenneProblem;
  class CouenneTNLP;

  class COUENNELIB_EXPORT SmartAsl : public Ipopt::ReferencedObject{
  public:
    ASL * asl;
    SmartAsl():
      Ipopt::ReferencedObject(),
      asl(NULL)
    {}
    virtual ~SmartAsl();
  };

  class COUENNELIB_EXPORT CouenneSetup : public Bonmin::BonminSetup{
  public:
    /** Default constructor*/
    CouenneSetup():
    BonminSetup(),
    aslfg_(NULL),
    CouennePtr_ (NULL),
    displayStats_ (false),
    couenneProb_ (NULL),
    couenneProb_is_own_(true) {}

    /** Copy constructor.*/
    CouenneSetup(const CouenneSetup& other):
      BonminSetup(other),
      aslfg_(NULL),
      displayStats_ (other.displayStats_),
      couenneProb_ (other.couenneProb_) {}

    /** virtual copy constructor.*/
    virtual Bonmin::BabSetupBase * clone () const
    {return new CouenneSetup (*this);}

    /// destructor
    virtual ~CouenneSetup();

    /** Initialize from command line arguments. */
    bool InitializeCouenne(char ** argv = NULL,
			   CouenneProblem *couenneProb = NULL,
			   Ipopt::SmartPtr<Bonmin::TMINLP> tminlp = NULL,
			   CouenneInterface *ci = NULL,
			   Bonmin::Bab *bb = NULL);

    /** the options */
    virtual void registerOptions();
    /** Register all Couenne options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

    /** Get the basic options if don't already have them.*/
    virtual void readOptionsFile(){
      if (readOptions_) return;
      Bonmin::BabSetupBase::readOptionsFile ("couenne.opt");
    }

    /// return pointer to cut generator (used to get pointer to problem)
    CouenneCutGenerator *couennePtr () const
    {return CouennePtr_;}

    /// true if one wants to display statistics at the end of program
    bool displayStats ()
    {return displayStats_;}

    /// add cut generators
    void addMilpCutGenerators ();

    /// modify parameter (used for MaxTime)
    inline void setDoubleParameter (const DoubleParameter &p, const double val)
    {doubleParam_ [p] = val;}

    /// modify parameter (used for MaxTime)
    inline double getDoubleParameter (const DoubleParameter &p) const
    {return doubleParam_ [p];}

    void setNodeComparisonMethod (Bonmin::BabSetupBase::NodeComparison c)
    {nodeComparisonMethod_ = c;}

private:
    Ipopt::SmartPtr<SmartAsl> aslfg_;

    /// hold a reference to Couenne cut generator to delete it at
    /// last. The alternative would be to clone it every time the
    /// CouenneSolverInterface is cloned (that is, at each call of
    /// Optimality-based bound tightening).
    CouenneCutGenerator *CouennePtr_;

    /// true if one wants to display statistics at the end of program
    bool displayStats_;

    /// MINLP formulation
    CouenneProblem *couenneProb_;

    /// whether the couenneProb_ has been created by Couenne, and thus will be deleted by Couenne
    bool couenneProb_is_own_;
  };
}

#endif
