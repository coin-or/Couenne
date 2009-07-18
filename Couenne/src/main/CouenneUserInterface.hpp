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

#ifndef COUENNEUSERINTERFACE_HPP_
#define COUENNEUSERINTERFACE_HPP_

#include "CouenneConfig.h"
#include "IpOptionsList.hpp"
#include "IpJournalist.hpp"

class CouenneProblem;
class CouenneBaB;
namespace Bonmin {
class TMINLP;
class Bab;
}

/* abstract base class of an interface for Couenne users
 */ 
class CouenneUserInterface {
protected:
	Ipopt::SmartPtr<Ipopt::OptionsList> options;
	Ipopt::SmartPtr<Ipopt::Journalist>  jnlst;
	
public:
	CouenneUserInterface(Ipopt::SmartPtr<Ipopt::OptionsList> options_, Ipopt::SmartPtr<Ipopt::Journalist> jnlst_)
	: options(options_), jnlst(jnlst_)
	{ }
	
	virtual ~CouenneUserInterface() { }
	
	virtual bool setupJournals() {
		Ipopt::SmartPtr<Ipopt::Journal> stdout_jrnl = jnlst->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
		stdout_jrnl->SetPrintLevel(Ipopt::J_DBG, Ipopt::J_NONE);
		return true;
	}
	
	virtual CouenneProblem* getCouenneProblem() = 0;
	
	virtual Ipopt::SmartPtr<Bonmin::TMINLP> getTMINLP() = 0;
	
  virtual bool addBabPlugins(Bonmin::Bab& bab) {
  	//  CutGenerator1   myCutGenerator   (problem, optionList);
  	//  BoundTightener1 myBoundTightener (problem, optionList);
  	//  Heuristic       myHeuristic      (problem, optionList);

  	//  bab.addCutGenerator   (myCutGenerator);
  	//  bab.addBoundTightener (myBoundTightener);
  	//  bab.addHeuristic      (myHeur1);
  	//  bab.addJournalist     (myJour);
  	return true;
  }

	virtual bool writeSolution(Bonmin::Bab& bab) {
		return true;
	}
};


#endif /*COUENNEUSERINTERFACE_HPP_*/
