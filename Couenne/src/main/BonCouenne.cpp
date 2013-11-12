// $Id$
//
// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pietro Belotti, Carnegie Mellon University,
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iomanip>
#include <fstream>

#include <stdlib.h>

#include "CoinTime.hpp"
#include "CoinError.hpp"
#include "BonCouenneInterface.hpp"

#include "BonCouenneSetup.hpp"

#include "BonCbc.hpp"
#include "CouenneBab.hpp"

#include "CbcCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

// the maximum difference between a printed optimum and a CouNumber
#define PRINTED_PRECISION 1e-5

#include "CouenneExprVar.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneProblem.hpp"
#include "CouenneJournalist.hpp"

#ifdef COIN_HAS_NTY
// FIXME: horrible global variable. Brrr.
#include "CouenneBranchingObject.hpp"
#endif

#include "CoinSignal.hpp"

#if 0
extern "C" {

  static int nInterrupts = 0;
  static void signal_handler (int sig) {

    if (!nInterrupts) {
      std::cerr << "[BREAK]" << std::endl;
      abort ();
    }
    return;
  }
}
#endif

//#define FM_FRES

int main (int argc, char *argv[]) {

  //CoinSighandler_t saveSignal = SIG_DFL;
  //saveSignal = signal (SIGINT, signal_handler);

    printf ("Couenne %s --  an Open-Source solver for Mixed Integer Nonlinear Optimization\n\
Mailing list: couenne@list.coin-or.org\n\
Instructions: http://www.coin-or.org/Couenne\n", 
	    strcmp (COUENNE_VERSION, "trunk") ? COUENNE_VERSION : "");

  WindowsErrorPopupBlocker();
  using namespace Ipopt;

  char * pbName = NULL;

  const int infeasible = 1;

  try {

    CouenneBab bb;
    //    bb.setUsingCouenne (true);

    CouenneProblem *p = NULL;
    CouenneInterface *ci = NULL;

#if 0
    //ci = new CouenneInterface;
    p = new CouenneProblem;

    p -> addVariable (false, p -> domain ());
    p -> addVariable (false, p -> domain ());
    p -> addVariable (false, p -> domain ());
    p -> addVariable (false, p -> domain ());

    p -> addObjective    (new exprSum (new exprClone (p->Var (1)), new exprClone (p->Var (2))), "min");
    p -> addLEConstraint (new exprSum (new exprClone (p->Var (0)), new exprClone (p->Var (2))), new exprConst (1));
    p -> addEQConstraint (new exprSum (new exprClone (p->Var (1)), new exprClone (p->Var (2))), new exprConst (1));
    p -> addEQConstraint (new exprSum (new exprClone (p->Var (1)), new exprClone (p->Var (3))), new exprConst (1));
    p -> addEQConstraint (new exprSum (new exprClone (p->Var (2)), new exprClone (p->Var (3))), new exprConst (1));
#endif

    CouenneSetup couenne;

    if (!couenne.InitializeCouenne (argv, p, NULL, ci, &bb))
      throw infeasible;

    //////////////////////////////
    CouenneCutGenerator *cg = NULL;

    // there is only one CouenneCutGenerator object; scan array until
    // dynamic_cast returns non-NULL

    if (couenne. cutGenerators () . size () > 0) {

      for (std::list <Bonmin::BabSetupBase::CuttingMethod>::iterator 
	           i  = couenne.cutGenerators () . begin ();
	   !cg && (i != couenne.cutGenerators () . end ()); 
	   ++i) 

	cg = dynamic_cast <CouenneCutGenerator *> (i -> cgl);
    }

    // This assumes that first cut generator is CouenneCutGenerator
    // CouenneCutGenerator *cg = dynamic_cast<CouenneCutGenerator *> 
    //   (couenne.cutGenerators().begin()->cgl);

    if (cg)

      cg -> setBabPtr (&bb);

    else {
      printf ("main(): ### ERROR: Can not get CouenneCutGenerator\n");
      exit (1);
    }

    // initial printout

    ConstJnlstPtr jnlst = couenne. couennePtr () -> Jnlst ();

    CouenneProblem *prob = couenne. couennePtr () -> Problem ();

    bb. setProblem (prob);

    jnlst -> Printf (J_ERROR, J_COUENNE, "\
Loaded instance \"%s\"\n\
Constraints:     %8d\n\
Variables:       %8d (%d integer)\n\
Auxiliaries:     %8d (%d integer)\n\n",
		     prob -> problemName ().c_str (),
		     prob -> nOrigCons (),
		     prob -> nOrigVars (),
		     prob -> nOrigIntVars (),
		     prob -> nVars    () - prob -> nOrigVars (),
		     CoinMax (0, prob -> nIntVars () - prob -> nOrigIntVars ()));

    double time_start = CoinCpuTime();

#if 0
    CouenneFeasibility feasibility;
    bb.model().setProblemFeasibility (feasibility);
#endif

    /// update time limit (read/preprocessing might have taken some)
    double timeLimit = 0;
    couenne.options () -> GetNumericValue ("time_limit", timeLimit, "couenne.");
    couenne.setDoubleParameter (Bonmin::BabSetupBase::MaxTime, 
				CoinMax (1., timeLimit - time_start));

    //jnlst -> Printf (J_ERROR, J_COUENNE, "Starting branch-and-bound\n");

    //////////////////////////////////

    bb (couenne); // do branch and bound

#ifdef COIN_HAS_NTY
    if (CouenneBranchingObject::nOrbBr)
      printf ("%d orbital nontrivial branchings\n", CouenneBranchingObject::nOrbBr);
#endif

    std::cout.precision (10);

    int nr=-1, nt=-1;
    double st=-1;

    if (cg) cg -> getStats (nr, nt, st);
    else printf ("Warning, could not get pointer to CouenneCutGenerator\n");

    CouenneProblem *cp = cg ? cg -> Problem () : NULL;

#if defined (FM_TRACE_OPTSOL) || defined (FM_FRES)
    double cbcLb = bb.model ().getBestPossibleObjValue();
    double printObj = 0;
    bool foundSol = false;
#endif

#ifdef FM_TRACE_OPTSOL

    //FILE *fSol = fopen ("bidon.sol", "w");

    // if(fSol == NULL) {
    //   printf("### ERROR: can not open bidon.sol\n");
    //   exit(1);
    // }

    if(cp != NULL) {
      double cbcObjVal = bb.model().getObjValue();
      int modelNvars = bb.model().getNumCols();

      CouenneRecordBestSol *rs = cp->getRecordBestSol(); 
      const double *cbcSol = bb.model().getColSolution();
      double *modCbcSol = new double[modelNvars];
      double modCbcSolVal= 1e100, modCbcSolMaxViol = 0;
      bool cbcSolIsFeas = false;

      if(modelNvars != cp->nVars()) {
	printf("### ERROR: modelNvars: %d nVars: %d\n", 
	       modelNvars, cp->nVars());
	exit(1);
      }

      if(cbcObjVal < 1e49) {

#ifdef FM_CHECKNLP2
	int cMS = rs->getCardModSol();
	cbcSolIsFeas = cp->checkNLP2(cbcSol, 0, false, // do not care about obj
				     false, // do not stop at first viol 
				     true, // checkAll 
				     cp->getFeasTol());
	CoinCopyN(rs->getModSol(cMS), cMS, modCbcSol);
	modCbcSolVal = rs->getModSolVal();
	modCbcSolMaxViol = rs->getModSolMaxViol();
#else /* not FM_CHECKNLP2 */
	int cMS = cp->nVars();
	cbcSolIsFeas = cp->checkNLP(cbcSol, modCbcSolVal, true);
	CoinCopyN(cbcSol, cMS, modCbcSol);
	modCbcSolMaxViol = cp->getFeasTol();
#endif /* not FM_CHECKNLP2 */
	foundSol = true;
      }

      const double *couenneSol = rs->getSol();
      double *modCouenneSol = new double[modelNvars];
      double modCouenneSolVal= 1e100, modCouenneSolMaxViol = 0;
      bool couenneSolIsFeas = false;

      if(couenneSol != NULL) {
#ifdef FM_CHECKNLP2
	couenneSolIsFeas = cp->checkNLP2(couenneSol, 0, false, 
					 false, true, 
					 cp->getFeasTol());
	int cMS = rs->getCardModSol();
	CoinCopyN(rs->getModSol(cMS), cMS, modCouenneSol);
	modCouenneSolVal = rs->getModSolVal();
	modCouenneSolMaxViol = rs->getModSolMaxViol();
#else /* not FM_CHECKNLP2 */
	couenneSolIsFeas = cp->checkNLP(couenneSol, modCouenneSolVal, true);
	int cMS = cp->nVars();
	CoinCopyN(couenneSol, cMS, modCouenneSol);
	modCouenneSolMaxViol = cp->getFeasTol();
#endif /* not FM_CHECKNLP2 */
	foundSol = true;
      }

      int retcomp = rs->compareAndSave(modCbcSol, modCbcSolVal, 
				       modCbcSolMaxViol, cbcSolIsFeas,
				       modCouenneSol, modCouenneSolVal, 
				       modCouenneSolMaxViol, couenneSolIsFeas, 
				       modelNvars, cp->getFeasTol());
      switch (retcomp) {
      case -1: printf("No solution found\n"); break;
      case 0: printf("Best solution found by Cbc. Value: %10.4f. Tolerance: %10g\n", modCbcSolVal, modCbcSolMaxViol); break;
      case 1: //printf("Best solution found by Couenne  Value: %10.4f  Tolerance: %10g\n", modCouenneSolVal, modCouenneSolMaxViol); break;
      default: break; // never happens
      }

      if(rs->getHasSol()) {
	if(cbcLb > rs->getVal()) { // Best sol found by Couenne and not
	                           // transmitted to Cbc
	  cbcLb = rs->getVal();
	}
	printObj = rs->getVal();
	//rs->printSol(fSol);
      }
      delete[] modCbcSol;
      delete[] modCouenneSol;
    }
    //fclose(fSol);
#endif /* FM_TRACE_OPTSOL */

#ifdef FM_FRES
    if(cp != NULL) {
      FILE *f_res = NULL;
      f_res = fopen("fres.xxx", "r");
      if(f_res == NULL) {
	f_res = fopen("fres.xxx", "w");
	fprintf(f_res, "END_OF_HEADER\n");
      }
      else {
	fclose(f_res);
	f_res = fopen("fres.xxx", "a");
      }
      char *pbName, shortName[256];  
      
      pbName = strdup(cp -> problemName ().c_str ());
      char *f_name_pos = strrchr(pbName, '/');
      if(f_name_pos != NULL) {
	strcpy(shortName, &(f_name_pos[1]));
      }
      else {
	strcpy(shortName, pbName);
      }
      

      fprintf(f_res, "%20s ", shortName);
      if((cbcLb > 1e20) || (cbcLb < -1e20)) {
        fprintf(f_res, "%10.4g", cbcLb);
      }
      else {
        fprintf(f_res, "%10.4f", cbcLb);
      }
      if(foundSol) {
	fprintf(f_res, " %10.4f", printObj);
      }
      else {
	fprintf(f_res, "         *");
      }
      fprintf(f_res, " %10d %10.4f\n", bb.numNodes (),
	      CoinCpuTime () - time_start);
      fclose(f_res);
    }
#endif

    // retrieve test value to check
    double global_opt;
    couenne.options () -> GetNumericValue ("couenne_check", global_opt, "couenne.");

    double 
      ub = bb. model (). getObjValue (),
      lb = bb. model (). getBestPossibleObjValue ();

    if (cp -> getRecordBestSol () &&
	(ub > cp -> getRecordBestSol () -> getVal ())) ub = cp -> getRecordBestSol () -> getVal ();

    if ((fabs (lb) > COUENNE_INFINITY / 1e4) ||
	(lb > ub))                                     lb = ub;

    char 
      *gapstr = new char [80],
      *lbstr  = new char [80],
      *ubstr  = new char [80];

    sprintf (lbstr,  "%10g",     lb);
    sprintf (ubstr,  "%10g",     ub);
    if(ub > COUENNE_INFINITY/1e4) {
      sprintf (gapstr, "--");
    }
    else {
      sprintf (gapstr, "%.2f%%", fabs (100. * (ub - lb) / (1. + fabs (lb))));
    }

    jnlst -> Printf (J_ERROR, J_COUENNE, "\n\
Linearization cuts added at root node:   %8d\n\
Linearization cuts added in total:       %8d  (separation time: %gs)\n\
Total solving time:                      %8gs (%gs in branch-and-bound)\n\
Lower bound:                           %s\n\
Upper bound:                           %s  (gap: %s)\n\
Branch-and-bound nodes:                  %8d\n",
		     nr, nt, st, 
		     CoinCpuTime () - time_start,
		     cg ? (CoinCpuTime () - CoinMax (time_start, cg -> rootTime ())) : CoinCpuTime () - time_start,
		     (lb < -9e12) ||
		     (lb > COUENNE_INFINITY/1e4) ? "      -inf" : lbstr,
		     (ub > COUENNE_INFINITY/1e4) ? "       inf" : ubstr,
		     (ub > COUENNE_INFINITY/1e4) ? "--" : gapstr,
		     bb.numNodes ());

    if (fabs (ub - bb. model (). getObjValue ()) > COUENNE_EPS * fabs (ub))
      jnlst -> Printf (J_ERROR, J_COUENNE, 
		       "Warning: upper bounds differ between Couenne and Cbc. Saving Couenne's (more reliable).\n");

    delete [] lbstr;
    delete [] ubstr;
    delete [] gapstr;

    if (global_opt < COUENNE_INFINITY) { // some value found in couenne.opt

      double opt = bb.model (). getBestPossibleObjValue ();

      printf ("Global Optimum Test on %-40s %s\n", 
	      cp ? cp -> problemName ().c_str () : "unknown", 
	      (fabs (opt - global_opt) / 
	       (1. + CoinMax (fabs (opt), fabs (global_opt))) < PRINTED_PRECISION) ? 
	      (const char *) "OK" : (const char *) "FAILED");
	      //opt, global_opt,
	      //fabs (opt - global_opt));

    } else // good old statistics

    if (couenne.displayStats ()) { // print statistics

      if (cg && !cp) printf ("Warning, could not get pointer to problem\n");
      else
	printf ("Stats: %-15s %4d [var] %4d [int] %4d [con] %4d [aux] "
		"%6d [root] %8d [tot] %6g [sep] %8g [time] %8g [bb] "
		"%20e [lower] %20e [upper] %7d [nodes]\n",// %s %s\n",
		cp ? cp -> problemName (). c_str () : "unknown",
		(cp) ? cp -> nOrigVars     () : -1, 
		(cp) ? cp -> nOrigIntVars  () : -1, 
		(cp) ? cp -> nOrigCons     () : -1,
		(cp) ? (cp -> nVars     () - 
			cp -> nOrigVars ()): -1,
		nr, nt, st, 
		CoinCpuTime () - time_start,
		cg ? (CoinCpuTime () - cg -> rootTime ()) : CoinCpuTime (),
		lb, //bb.model (). getBestPossibleObjValue (),
		ub, //bb.model (). getObjValue (),
		//bb.bestBound (),
		//bb.bestObj (),
		bb.numNodes ());
		//bb.iterationCount ());
		//status.c_str (), message.c_str ());
    }

//    nlp_and_solver -> writeAmplSolFile (message, bb.bestSolution (), NULL);
  }
  catch(Bonmin::TNLPSolver::UnsolvedError *E) {
     E->writeDiffFiles();
     E->printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
  }
  catch (Bonmin::OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch (CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch (Ipopt::OPTION_INVALID &E)
  {
   std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
  }
  catch (int generic_error) {
    if (generic_error == infeasible)
      printf ("problem infeasible\n");
  }

//  catch(...) {
//    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
//    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
//    throw;
//  }

  delete [] pbName;
  return 0;
}
