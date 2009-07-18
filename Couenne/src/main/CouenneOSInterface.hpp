// $Id$
//
// (C) Copyright XXX 2009
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belotti, Lehigh University
// Stefan Vigerske, Humboldt University
//
// Date : 07/18/2009

#ifndef COUENNEOSINTERFACE_HPP_
#define COUENNEOSINTERFACE_HPP_

#include "CouenneUserInterface.hpp"
#include "BonRegisteredOptions.hpp"

class OSInstance;

class CouenneOSInterface : public CouenneUserInterface {
private:
	CouenneProblem*                  problem;
	Ipopt::SmartPtr<Bonmin::TMINLP>  tminlp;
	
	OSInstance*                      osinstance;
	
public:
	static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
	
	CouenneOSInterface(Ipopt::SmartPtr<Ipopt::OptionsList> options_, Ipopt::SmartPtr<Ipopt::Journalist> jnlst_)
	: CouenneUserInterface(options_, jnlst_), problem(NULL), osinstance(NULL)
	{ }
	
	~CouenneOSInterface();
	
	CouenneProblem* getCouenneProblem();
	
	Ipopt::SmartPtr<Bonmin::TMINLP> getTMINLP();
	
	bool writeSolution(Bonmin::Bab& bab);
};

#endif /*COUENNEOSINTERFACE_HPP_*/
