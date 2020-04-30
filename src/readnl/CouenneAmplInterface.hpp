// $Id$
//
// (C) Copyright XXX 2009
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pietro Belotti, Lehigh University
// Stefan Vigerske, Humboldt University
//
// Date : 07/18/2009

#ifndef COUENNEAMPLINTERFACE_HPP_
#define COUENNEAMPLINTERFACE_HPP_

#include "CouenneUserInterface.hpp"
#include "BonRegisteredOptions.hpp"

struct ASL;
struct expr;

namespace Couenne {

class expression;

class CouenneAmplInterface : public CouenneUserInterface {
private:
	CouenneProblem*                  problem;
	Ipopt::SmartPtr<Bonmin::TMINLP>  tminlp;
	
	Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions;
	struct ASL*                      asl;
	
	bool readASLfg();
	bool readnl();
	expression* nl2e(expr *e);
	
public:
	static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
	
	CouenneAmplInterface(Ipopt::SmartPtr<Ipopt::OptionsList> options_, Ipopt::SmartPtr<Ipopt::Journalist> jnlst_)
	: CouenneUserInterface(options_, jnlst_), problem(NULL), asl(NULL)
	{ }
	
	~CouenneAmplInterface();
	
	CouenneProblem* getCouenneProblem();
	
	Ipopt::SmartPtr<Bonmin::TMINLP> getTMINLP();
	
	bool writeSolution(Bonmin::Bab& bab);
	
	void setRegisteredOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions_)
	{ roptions = roptions_; }
	
};

}

#endif /*COUENNEAMPLINTERFACE_HPP_*/
