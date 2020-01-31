// $Id$
//
// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pietro Belotti, Lehigh University
// Stefan Vigerske, Humboldt University
//
// Date : 07/06/2009

//#include "CouenneConfig.h"

#include <cstdlib>

#include "CoinPragma.hpp"
#include "CoinError.hpp"
#include "CoinTime.hpp"

#include "CouenneUserInterface.hpp"
#ifdef COIN_HAS_ASL
#include "CouenneAmplInterface.hpp"
#endif
#ifdef COIN_HAS_OS
#include "CouenneOSInterface.hpp"
#endif

#include "BonRegisteredOptions.hpp"
#include "BonCbc.hpp"

#include "BonCouenneSetup.hpp"
#include "BonCouenneInterface.hpp"

// for printing of statistics
#include "CbcCutGenerator.hpp"      
#include "CouenneCutGenerator.hpp" 
#include "CouenneProblem.hpp"

using namespace Couenne;

// the maximum difference between a printed optimum and a CouNumber
#define PRINTED_PRECISION 1e-5

using Ipopt::SmartPtr;

static const int infeasible = 1;

bool parseCommandLine(int argc, char* argv[], Ipopt::SmartPtr<Ipopt::OptionsList> options) {
	assert(IsValid(options));
	
	if (argc==3 && strcmp(argv[1], "-AMPL")==0)
		options->SetStringValue("nlfile", argv[2]);

	if (argc==3 && strcmp(argv[1], "-OSIL")==0)
		options->SetStringValue("osilfile", argv[2]);

	return true;
}

int main (int argc, char *argv[]) {
  WindowsErrorPopupBlocker();

  double time_start = CoinCpuTime();
  
	// register options to prepare for parsing the command line
	SmartPtr<Bonmin::RegisteredOptions> roptions = new Bonmin::RegisteredOptions();
	Bonmin::CouenneSetup::registerAllOptions(roptions);
#ifdef COIN_HAS_ASL
	CouenneAmplInterface::registerOptions(roptions);
#endif
#ifdef COIN_HAS_OS
	CouenneOSInterface::registerOptions(roptions);
#endif
	
	SmartPtr<Ipopt::Journalist> jnlst = new Ipopt::Journalist();
	// do not add journals yet, maybe the user wants to do so; but what if parsing the command line gives errors?
	
	SmartPtr<Ipopt::OptionsList> options = new Ipopt::OptionsList(GetRawPtr(roptions), jnlst);
	if (!parseCommandLine(argc, argv, options))
		return EXIT_FAILURE;
	
	

	CouenneUserInterface* userinterface = NULL;
	
	std::string dummy;
#ifdef COIN_HAS_ASL
	if (!userinterface && options->GetStringValue("nlfile", dummy, "")) {
		userinterface = new CouenneAmplInterface(options, jnlst);
		((CouenneAmplInterface*)userinterface) -> setRegisteredOptions(roptions); // for some reason the TMINLP constructor needs the registered options
	}
#endif
#ifdef COIN_HAS_OS
	if (!userinterface && options->GetStringValue("osilfile", dummy, "")) {
		userinterface = new CouenneOSInterface();
	}
#endif
	
	if (!userinterface) {
		fprintf(stderr, "Error: No input file given.\n");
		return EXIT_FAILURE;
	}
	
	if (!userinterface->setupJournals())
		return EXIT_FAILURE;
	
	CouenneProblem* problem = userinterface->getCouenneProblem();
	if (!problem)
		return EXIT_FAILURE;
	problem->initOptions(options);
	
	SmartPtr<Bonmin::TMINLP> tminlp = userinterface->getTMINLP();
	if (Ipopt::IsNull(tminlp))
		return EXIT_FAILURE;

	try {
    Bonmin::Bab bb;
    bb.setUsingCouenne (true);

    Bonmin::CouenneSetup couenne;
    couenne.setOptionsAndJournalist(roptions, options, jnlst);
    if (!couenne.InitializeCouenne (NULL, problem, tminlp))
      throw infeasible;

    double timeLimit = 0;
    options -> GetNumericValue ("time_limit", timeLimit, "couenne.");
    couenne.setDoubleParameter (Bonmin::BabSetupBase::MaxTime, timeLimit - (time_start = (CoinCpuTime () - time_start)));
  
    if (!userinterface->addBabPlugins(bb))
    	return EXIT_FAILURE;

    bb (couenne); // do branch and bound
    
    // retrieve test value to check
    double global_opt;
    options -> GetNumericValue ("couenne_check", global_opt, "couenne.");

    if (global_opt < COUENNE_INFINITY) { // some value found in couenne.opt
      double opt = bb.model (). getBestPossibleObjValue ();

      jnlst -> Printf(Ipopt::J_SUMMARY, J_PROBLEM, "Global Optimum Test on %-40s %s\n", 
	      problem -> problemName ().c_str (), 
	      (fabs (opt - global_opt) / 
	       (1. + CoinMax (fabs (opt), fabs (global_opt))) < PRINTED_PRECISION) ? 
	      "OK" : "FAILED");

    } else if (couenne.displayStats ()) { // print statistics

      int nr=-1, nt=-1;
      double st=-1;

      CouenneCutGenerator* cg = NULL;
      if (bb.model (). cutGenerators ())
        cg = dynamic_cast <CouenneCutGenerator *>	(bb.model (). cutGenerators () [0] -> generator ());
      if (cg) cg -> getStats (nr, nt, st);
      else jnlst -> Printf(Ipopt::J_WARNING, J_PROBLEM, "Warning: Could not get pointer to CouenneCutGenerator\n");

    	jnlst -> Printf(Ipopt::J_SUMMARY, J_PROBLEM, "Stats: %-15s %4d [var] %4d [int] %4d [con] %4d [aux] "
    			"%6d [root] %8d [tot] %6g [sep] %8g [time] %8g [bb] "
    			"%20e [lower] %20e [upper] %7d [nodes]\n",// %s %s\n",
    			problem -> problemName ().c_str (),
    			problem -> nOrigVars   (), 
    			problem -> nOrigIntVars(), 
    			problem -> nOrigCons   (),
    			problem -> nVars       () - problem -> nOrigVars (),
    			nr, nt, st, 
    			CoinCpuTime () - time_start,
    			cg ? (CoinCpuTime () - cg -> rootTime ()) : CoinCpuTime (),
    			bb.model (). getBestPossibleObjValue (),
    			bb.model (). getObjValue (),
    			//bb.bestBound (),
    			//bb.bestObj (),
    			bb.numNodes ()
    			//bb.iterationCount (),
    			//status.c_str (), message.c_str ()
    	);
    }    

    if (!userinterface->writeSolution(bb))
    	return EXIT_FAILURE;
 	
	} catch(Bonmin::TNLPSolver::UnsolvedError *E) {
     E->writeDiffFiles();
     E->printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
     
  } catch(Bonmin::OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
    
  } catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
    
  } catch (Ipopt::OPTION_INVALID &E) {
  	std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
   
  } catch (int generic_error) {
    if (generic_error == infeasible)
      jnlst->Printf(Ipopt::J_SUMMARY, J_PROBLEM, "problem infeasible\n");
  }
  
  delete userinterface;
  
  return EXIT_SUCCESS;
}
