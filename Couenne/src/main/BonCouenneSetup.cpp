// $Id$
//
// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007

// Ampl includes

#include "CouenneConfig.h"

#include "OsiClpSolverInterface.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_GRB
#include "OsiGrbSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SPX
#include "OsiSpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_XPR
#include "OsiXprSolverInterface.hpp"
#endif

// MILP cuts
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"
#include "CglLandP.hpp"
#include "CglRedSplit.hpp"

#include "BonCouenneSetup.hpp"
#include "CouenneFeasPump.hpp"
#include "CouenneIterativeRounding.hpp"
#include "BonCouenneInterface.hpp"
#include "BonInitHeuristic.hpp"
#include "BonNlpHeuristic.hpp"

#include "Heuristics/BonFixAndSolveHeuristic.hpp"
#include "Heuristics/BonDummyPump.hpp"
#include "Heuristics/BonPumpForMinlp.hpp"
#include "Heuristics/BonHeuristicRINS.hpp"
#include "Heuristics/BonHeuristicLocalBranching.hpp"
#include "Heuristics/BonHeuristicFPump.hpp"
#include "Heuristics/BonHeuristicDiveFractional.hpp"
#include "Heuristics/BonHeuristicDiveVectorLength.hpp"
#include "Heuristics/BonHeuristicDiveMIPFractional.hpp"
#include "Heuristics/BonHeuristicDiveMIPVectorLength.hpp"
#include "Heuristics/BonMilpRounding.hpp"

#include "BonGuessHeuristic.hpp"
#include "CbcCompareActual.hpp"

#include "CouenneObject.hpp"
#include "CouenneVarObject.hpp"
#include "CouenneVTObject.hpp"
#include "CouenneOrbitObj.hpp"
#include "CouenneChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneSolverInterface.hpp"
#include "CouenneFixPoint.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneDisjCuts.hpp"
//#include "CouenneCrossConv.hpp"
#include "CouenneTwoImplied.hpp"

#include "BonCouenneInfo.hpp"
#include "BonCbcNode.hpp"
#include "BonCbc.hpp"

// ASL includes need to come behind OsiClp and Bonmin, because it defines "filename",
// which is used as variablename in Clp
// (similar bad as windows.h, which defines "small")
#ifdef COIN_HAS_ASL
#include "asl.h"
#include "getstub.h"
#endif

using namespace Ipopt;
using namespace Couenne;
  
CouenneSetup::~CouenneSetup(){
  if (couenneProb_ && couenneProb_is_own_)
    delete couenneProb_;

#ifdef COIN_HAS_ASL
  // free (aslfg_ -> asl); // triggers segfault
#endif
}

bool CouenneSetup::InitializeCouenne (char ** argv,
				      CouenneProblem *couenneProb,
				      Ipopt::SmartPtr<Bonmin::TMINLP> tminlp,
				      CouenneInterface *ci,
				      Bonmin::Bab *bb) {

  std::string s;

  if (couenneProb) {
    //TODO create a copy of user problem, since we modify it? 
    couenneProb_ = couenneProb;
    couenneProb_is_own_ = false;
  }

  /* Get the basic options. */
  readOptionsFile();
 
  // Suppress iteration output from nonlinear solver
  options () -> SetStringValue ("sb", "yes", false, true);

  // in check mode, avoid pop-up error message (there are quite a few messages)
  options_ -> GetStringValue ("test_mode", s, "couenne.");
  if (s == "yes")
    WindowsErrorPopupBlocker();

  /** Change default value for failure behavior so that code doesn't crash 
      when Ipopt does not solve a sub-problem.*/

  options_ -> SetStringValue ("nlp_failure_behavior", "fathom", "couenne.");

  gatherParametersValues (options_);

  if (!ci) {

    ci = new CouenneInterface;

    if (!couenneProb_ && argv) {
#ifdef COIN_HAS_ASL
      /* Read the model in various places. */
      ci -> readAmplNlFile (argv, roptions (), options (), journalist ());
      aslfg_ = new SmartAsl;
      aslfg_ -> asl = readASLfg (argv);
#else
      std::cerr << 
	"Couenne was compiled without AMPL Solver Library. Cannot initialize from AMPL NL File." 
		<< std::endl;
      return false;
#endif
    } else {
      assert (couenneProb_ != NULL);
      assert (IsValid (tminlp)); //TODO would be great to setup own TMINLP based on CouenneProblem formulation
      ci -> initialize (roptions_, options_, journalist_, tminlp);
    }
  }

  nonlinearSolver_ = ci;

  /** Set the output level of the journalist for all Couenne
      categories.  We probably want to make that a bit more flexible
      later. */

  int i;

  /// trying to avoid repetitions here...

  // FIXME: doesn't work if suppress_all_output is true (gives
  // segfault on options(), but checking options()!=NULL won't work as
  // options() is a SmartPtr

#define addJournalist(optname,jlevel) {				\
    options    () -> GetIntegerValue ((optname), i, "couenne."); \
    journalist () -> GetJournal      ("console") -> SetPrintLevel ((jlevel), (EJournalLevel) i); \
}

  addJournalist ("output_level",                J_COUENNE);
  addJournalist ("boundtightening_print_level", J_BOUNDTIGHTENING);
  addJournalist ("branching_print_level",       J_BRANCHING);
  addJournalist ("convexifying_print_level",    J_CONVEXIFYING);
  addJournalist ("problem_print_level",         J_PROBLEM);
  addJournalist ("nlpheur_print_level",         J_NLPHEURISTIC);
  addJournalist ("disjcuts_print_level",        J_DISJCUTS);
  addJournalist ("reformulate_print_level",     J_REFORMULATE);

  /* Initialize Couenne cut generator.*/
  //int ivalue, num_points;
  //options()->GetEnumValue("convexification_type", ivalue,"couenne.");
  //options()->GetIntegerValue("convexification_points",num_points,"couenne.");

  if (!couenneProb_)
    couenneProb_ = new CouenneProblem (aslfg_ -> asl, this, journalist ());

  CouenneCutGenerator * couenneCg = 
    new CouenneCutGenerator (ci, this, couenneProb_, NULL);

  options_ -> GetStringValue ("lp_solver", s, "couenne.");

  if (s == "clp") {

    CouenneSolverInterface <OsiClpSolverInterface> *CSI = new CouenneSolverInterface <OsiClpSolverInterface>;
    continuousSolver_ = CSI;
    CSI -> setCutGenPtr (couenneCg);

  } else if (s == "cplex") {

#ifdef COIN_HAS_CPX
    CouenneSolverInterface <OsiCpxSolverInterface> *CSI = new CouenneSolverInterface <OsiCpxSolverInterface>;
    continuousSolver_ = CSI;
    CSI -> setCutGenPtr (couenneCg);
#else
    journalist()->Printf(J_ERROR, J_INITIALIZATION, "Couenne was compiled without CPLEX interface. Please reconfigure, recompile, and try again.\n");
    return false;
#endif
  } else if (s == "xpress-mp") {

#ifdef COIN_HAS_XPR
    CouenneSolverInterface <OsiXprSolverInterface> *CSI = new CouenneSolverInterface <OsiXprSolverInterface>;
    continuousSolver_ = CSI;
    CSI -> setCutGenPtr (couenneCg);
#else
    journalist()->Printf(J_ERROR, J_INITIALIZATION, "Couenne was compiled without Xpress-MP interface. Please reconfigure, recompile, and try again.\n");
    return false;
#endif
  } else if (s == "gurobi") {

#ifdef COIN_HAS_GRB
    CouenneSolverInterface <OsiGrbSolverInterface> *CSI = new CouenneSolverInterface <OsiGrbSolverInterface>;
    continuousSolver_ = CSI;
    CSI -> setCutGenPtr (couenneCg);
#else
    journalist()->Printf(J_ERROR, J_INITIALIZATION, "Couenne was compiled without GUROBI interface. Please reconfigure, recompile, and try again.\n");
    return false;
#endif
  } else if (s == "soplex") {

#ifdef COIN_HAS_SPX
    CouenneSolverInterface <OsiSpxSolverInterface> *CSI = new CouenneSolverInterface <OsiSpxSolverInterface>;
    continuousSolver_ = CSI;
    CSI -> setCutGenPtr (couenneCg);
#else
    journalist()->Printf(J_ERROR, J_INITIALIZATION, "Couenne was compiled without Soplex. Please reconfigure, recompile, and try again.\n");
    return false;
#endif
  } else {
    journalist ()-> Printf (J_ERROR, J_INITIALIZATION, "The LP solver you specified hasn't been added to Couenne yet.\n");
    return false;
  }

  continuousSolver_ -> passInMessageHandler(ci -> messageHandler());

  couenneProb_ -> setBase (this);

  assert (couenneProb_);

  couenneProb_ -> reformulate (couenneCg);

  Bonmin::BabInfo * extraStuff = new CouenneInfo (0);

  // as per instructions by John Forrest, to get changed bounds
  extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);

  continuousSolver_ -> setAuxiliaryInfo (extraStuff);
  delete extraStuff;
    
  extraStuff = dynamic_cast <Bonmin::BabInfo *> (continuousSolver_ -> getAuxiliaryInfo ());

  // Setup log level of LP solver
  int lpLogLevel;
  options () -> GetIntegerValue ("lp_log_level", lpLogLevel, "couenne.");
  continuousSolver_ -> messageHandler () -> setLogLevel (lpLogLevel);

  //////////////////////////////////////////////////////////////

  couenneCg -> Problem () -> setMaxCpuTime (getDoubleParameter (BabSetupBase::MaxTime));

  ci -> extractLinearRelaxation (*continuousSolver_, *couenneCg);

  // In case there are no discrete variables, we have already a
  // heuristic solution for which create a initialization heuristic
  if (!(extraStuff -> infeasibleNode ()) &&
      ci -> isProvenOptimal () && 
      ci -> haveNlpSolution ()) {

    /// setup initial heuristic (in principle it should only run once...)
    InitHeuristic* initHeuristic = new InitHeuristic 
      (ci -> getObjValue (), ci -> getColSolution (), *couenneProb_);
    HeuristicMethod h;
    h.id = "Couenne Rounding NLP"; // same name as the real rounding one
    h.heuristic = initHeuristic;
    heuristics_.push_back(h);
  }

  if (extraStuff -> infeasibleNode ()){
    journalist() -> Printf(J_SUMMARY, J_PROBLEM, "Initial linear relaxation constructed by Couenne is infeasible, exiting...\n");
    return false;
  }

  //continuousSolver_ -> findIntegersAndSOS (false);
  //addSos (); // only adds embedded SOS objects

  // Add Couenne SOS ///////////////////////////////////////////////////////////////

  int 
    nSOS  = 0,
    nVars = couenneProb_ -> nVars ();

  OsiObject ** objects = NULL;

  options () -> GetStringValue ("enable_sos", s, "couenne.");

  if (s == "yes") {

    // allocate sufficient space for both nonlinear variables and SOS's
    objects = new OsiObject* [couenneProb_ -> nCons () + nVars];

    nSOS = couenneProb_ -> findSOS (&(bb -> model()), dynamic_cast <OsiSolverInterface *> (nonlinearSolver ()), objects);

    nonlinearSolver () -> addObjects (nSOS, objects);

    //printf ("==================== found %d SOS\n", nSOS);
    //nonlinearSolver () -> addObjects (nSOS, objects);
    //continuousSolver () -> addObjects (nSOS, objects);

    //printf ("found %d SOS!\n", nSOS);

    /*for (int i=0; i<nSOS; i++)
      delete objects [i];
      delete [] objects;*/

    if (!nSOS) {
      delete [] objects;
      objects = NULL;
    } 
  }

  //model -> assignSolver (continuousSolver_, true);
  //continuousSolver_ = model -> solver();

  // Add Couenne objects for branching /////////////////////////////////////////////

  options () -> GetStringValue ("display_stats", s, "couenne.");
  displayStats_ = (s == "yes");

  options () -> GetStringValue ("branching_object", s, "couenne.");

  enum CouenneObject::branch_obj objType = CouenneObject::VAR_OBJ;

  if      (s ==   "vt_obj") objType = CouenneObject::VT_OBJ;
  else if (s ==  "var_obj") objType = CouenneObject::VAR_OBJ;
  else if (s == "expr_obj") objType = CouenneObject::EXPR_OBJ;
  else {
    printf ("CouenneSetup: Unknown branching object type\n");
    exit (-1);
  }

  int 
    nobj = nSOS; // if no SOS then objects is empty

  if (!objects)
    objects = new OsiObject* [nVars];

  int 
    contObjPriority, 
    intObjPriority;

  options () -> GetIntegerValue ("cont_var_priority", contObjPriority, "couenne.");
  options () -> GetIntegerValue ( "int_var_priority",  intObjPriority, "couenne.");

  int varSelection;
  if (!options_ -> GetEnumValue ("variable_selection", varSelection, "couenne.")) {
    // change the default for Couenne
    varSelection = Bonmin::BabSetupBase::OSI_SIMPLE;
  }

  for (int i = 0; i < nVars; i++) { // for each variable

    exprVar *var = couenneProb_ -> Var (i);

    // we only want enabled variables
    if (var -> Multiplicity () <= 0) 
      continue;

    switch (objType) {

    case CouenneObject::EXPR_OBJ:

      // if this variable is associated with a nonlinear function
      if (var -> isInteger () || 
	  ((var -> Type  () == AUX) && 
	   (var -> Image () -> Linearity () > LINEAR))) {

	/*if ((var -> Image () -> code () == COU_EXPRMUL) &&
	  (var -> Image () -> ArgList () [0] -> Index () >= 0) &&
	  (var -> Image () -> ArgList () [1] -> Index () >= 0) &&
	  (fabs (var -> lb ()) < COUENNE_EPS) &&
	  (fabs (var -> ub ()) < COUENNE_EPS))

	  // it's a complementarity constraint object!
	  objects    [nobj] = new CouenneComplObject (couenneProb_, var, this, journalist ());
	  else*/
	objects [nobj] = new CouenneObject (couenneCg, couenneProb_, var, this, journalist ());

	objects [nobj++] -> setPriority (var -> isInteger () ? intObjPriority : contObjPriority);
	//objects [nobj++] -> setPriority (contObjPriority + var -> rank ());
      }

      break;

    case CouenneObject::VAR_OBJ:

      // branching objects on variables
      if // comment three lines below for linear variables too
	(var -> isInteger () || 
	 (couenneProb_ -> Dependence () [var -> Index ()] . size () > 0)) {  // has indep
	//|| ((var -> Type () == AUX) &&                                  // or, aux 
	//    (var -> Image () -> Linearity () > LINEAR))) {              // of nonlinear

	objects [nobj] = new CouenneVarObject (couenneCg, couenneProb_, var, this, journalist (), varSelection);
	objects [nobj++] -> setPriority (var -> isInteger () ? intObjPriority : contObjPriority);
	//objects [nobj++] -> setPriority (contObjPriority + var -> rank ());
      }

      break;

    default:
    case CouenneObject::VT_OBJ:

      // branching objects on variables
      if // comment three lines below for linear variables too
	(var -> isInteger () || 
	 (couenneProb_ -> Dependence () [var -> Index ()] . size () > 0)) { // has indep
	//|| ((var -> Type () == AUX) &&                      // or, aux 
	//(var -> Image () -> Linearity () > LINEAR))) { // of nonlinear

	objects [nobj] = new CouenneVTObject (couenneCg, couenneProb_, var, this, journalist ());
	objects [nobj++] -> setPriority (var -> isInteger () ? intObjPriority : contObjPriority);
	//objects [nobj++] -> setPriority (contObjPriority + var -> rank ());
      }

      break;
    }
  }

  // // Experimental: orbital branching //////////////////////////////////////////////

  // options () -> GetStringValue ("orbital_branching", s, "couenne.");

  // if (s == "yes") {

  //   objects [nobj] = new CouenneOrbitObj (couenneCg, couenneProb_, NULL, this, journalist ());
  //   objects [nobj++] -> setPriority (contObjPriority);
  // }

  // Add objects /////////////////////////////////

  continuousSolver_ -> addObjects (nobj, objects);

  for (int i = 0 ; i < nobj ; i++)
    delete objects [i];

  delete [] objects;

  int freq;

  // Setup Fix Point bound tightener /////////////////////////////////////////////

  options () -> GetIntegerValue ("fixpoint_bt", freq, "couenne.");

  if (freq != 0) {

    CuttingMethod cg;
    cg.frequency = freq;
    cg.cgl = new CouenneFixPoint (couenneProb_, options ());
    cg.id = "Couenne fixed point FBBT";
    cutGenerators (). push_back (cg);
  }

  // Setup Convexifier generators ////////////////////////////////////////////////

  options () -> GetIntegerValue ("convexification_cuts", freq, "couenne.");

  if (freq != 0) {

    CuttingMethod cg;
    cg.frequency = freq;
    cg.cgl = couenneCg;
    cg.id = "Couenne convexifier cuts";
    cutGenerators().push_back (cg);

    // set cut gen pointer
    //dynamic_cast <CouenneSolverInterface <OsiClpSolverInterface> *> 
    //(continuousSolver_)

    // this is done on an explicitly declared CSI pointer, however
    // CSI == continuousSolver_
  }

  // add other cut generators -- test for integer variables first
  if (couenneCg -> Problem () -> nIntVars () > 0)
    addMilpCutGenerators ();

  CouennePtr_ = couenneCg;

  // Add two-inequalities based bound tightening ///////////////////////////////////////////////////////

  options () -> GetIntegerValue ("two_implied_bt", freq, "couenne.");

  if (freq != 0) {

    CouenneTwoImplied * couenne2I = 
      new CouenneTwoImplied (couenneProb_,
			     journalist (),
			     options    ());
    CuttingMethod cg;
    cg.frequency = freq;
    cg.cgl = couenne2I;
    cg.id = "Couenne two-implied cuts";
    cutGenerators (). push_back(cg);
  }

  // check branch variable selection for disjunctive cuts

  // Setup heuristic to solve MINLP problems. /////////////////////////////////

  std::string doHeuristic;

  options () -> GetStringValue ("local_optimization_heuristic", doHeuristic, "couenne.");

  if (doHeuristic == "yes") {

    int numSolve;
    options()->GetIntegerValue("log_num_local_optimization_per_level",numSolve,"couenne.");
    NlpSolveHeuristic * nlpHeuristic = new NlpSolveHeuristic;
    nlpHeuristic->setNlp(*ci,false);
    nlpHeuristic->setCouenneProblem(couenneProb_);
    nlpHeuristic->setMaxNlpInf(maxNlpInf_0);
    nlpHeuristic->setNumberSolvePerLevel(numSolve);
    HeuristicMethod h;
    h.id = "Couenne Rounding NLP";
    h.heuristic = nlpHeuristic;
    heuristics_.push_back(h);
  }

  options () -> GetStringValue ("iterative_rounding_heuristic", doHeuristic, "couenne.");
  
  if (doHeuristic == "yes") {
    CouenneIterativeRounding * nlpHeuristic = new CouenneIterativeRounding(nonlinearSolver_, ci, couenneProb_, options());
    HeuristicMethod h;
    h.id = "Couenne Iterative Rounding";
    h.heuristic = nlpHeuristic;
    heuristics_.push_back(h);
  }

  options () -> GetStringValue ("feas_pump_heuristic", doHeuristic, "couenne.");

  if (doHeuristic == "yes") {

    int numSolve;
    options () -> GetIntegerValue ("feas_pump_level", numSolve, "couenne.");

    CouenneFeasPump *nlpHeuristic = new CouenneFeasPump (couenneProb_, couenneCg, options ());

    nlpHeuristic -> setNumberSolvePerLevel (numSolve);

    HeuristicMethod h;

    h.id = "Couenne Feasibility Pump";
    h.heuristic = nlpHeuristic;
    heuristics_. push_back (h);
  }


  if (0) { // inactive as of yet -- segfaults 

    Ipopt::Index doHeuristicDiveFractional = false;
    options()->GetEnumValue("heuristic_dive_fractional",doHeuristicDiveFractional,prefix_.c_str());
    if(doHeuristicDiveFractional){
      Bonmin::HeuristicDiveFractional* dive_fractional = new Bonmin::HeuristicDiveFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_fractional;
      h.id = "DiveFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveVectorLength = false;
    options()->GetEnumValue("heuristic_dive_vectorLength",doHeuristicDiveVectorLength,prefix_.c_str());
    if(doHeuristicDiveVectorLength){
      Bonmin::HeuristicDiveVectorLength* dive_vectorLength = new Bonmin::HeuristicDiveVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_vectorLength;
      h.id = "DiveVectorLength";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPFractional = false;
    if(!options()->GetEnumValue("heuristic_dive_MIP_fractional",doHeuristicDiveMIPFractional,prefix_.c_str())){
      doHeuristicDiveMIPFractional = true;
      std::string o_name = prefix_ + "heuristic_dive_MIP_fractional";
      options_->SetStringValue(o_name.c_str(), "yes",true,true);
    }
    if(doHeuristicDiveMIPFractional){
      Bonmin::HeuristicDiveMIPFractional* dive_MIP_fractional = new Bonmin::HeuristicDiveMIPFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_fractional;
      h.id = "DiveMIPFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPVectorLength = false;
    options()->GetEnumValue("heuristic_dive_MIP_vectorLength",doHeuristicDiveMIPVectorLength,prefix_.c_str());
    if(doHeuristicDiveMIPVectorLength){
      Bonmin::HeuristicDiveMIPVectorLength* dive_MIP_vectorLength = new Bonmin::HeuristicDiveMIPVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_vectorLength;
      h.id = "DiveMIPVectorLength";
      heuristics_.push_back(h);
    }
    Ipopt::Index doHeuristicFPump = false;
    if(!nonlinearSolver_->model()->hasGeneralInteger() && !options()->GetEnumValue("heuristic_feasibility_pump",doHeuristicFPump,prefix_.c_str())){
      doHeuristicFPump = true;
      std::string o_name = prefix_ + "heuristic_feasibility_pump";
      options_->SetStringValue(o_name.c_str(), "yes",true,true);
    }
    if(doHeuristicFPump){
      Bonmin::HeuristicFPump* feasibility_pump = new Bonmin::HeuristicFPump(this);
      HeuristicMethod h;
      h.heuristic = feasibility_pump;
      h.id = "FPump";
      heuristics_.push_back(h);
    }

    Ipopt::Index doFixAndSolve = false;
    options()->GetEnumValue("fix_and_solve_heuristic",doFixAndSolve,prefix_.c_str());
    if(doFixAndSolve){
      Bonmin::FixAndSolveHeuristic* fix_and_solve = new Bonmin::FixAndSolveHeuristic(this);
      HeuristicMethod h;
      h.heuristic = fix_and_solve;
      h.id = "Fix and Solve";
      heuristics_.push_back(h);
    }

    Ipopt::Index doDummyPump = false;
    options()->GetEnumValue("dummy_pump_heuristic",doDummyPump,prefix_.c_str());
    if(doDummyPump){
      Bonmin::DummyPump* fix_and_solve = new Bonmin::DummyPump(this);
      HeuristicMethod h;
      h.heuristic = fix_and_solve;
      h.id = "Dummy pump";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicRINS = false;
    options()->GetEnumValue("heuristic_RINS",doHeuristicRINS,prefix_.c_str());
    if(doHeuristicRINS){
      Bonmin::HeuristicRINS* rins = new Bonmin::HeuristicRINS(this);
      HeuristicMethod h;
      h.heuristic = rins;
      h.id = "RINS";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicLocalBranching = false;
    options()->GetEnumValue("heuristic_local_branching",doHeuristicLocalBranching,prefix_.c_str());
    if(doHeuristicLocalBranching){
      Bonmin::HeuristicLocalBranching* local_branching = new Bonmin::HeuristicLocalBranching(this);
      HeuristicMethod h;
      h.heuristic = local_branching;
      h.id = "LocalBranching";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicPumpForMinlp = false;
    options()->GetEnumValue("pump_for_minlp",doHeuristicPumpForMinlp,prefix_.c_str());
    if(doHeuristicPumpForMinlp){
      Bonmin::PumpForMinlp * pump = new Bonmin::PumpForMinlp(this);
      HeuristicMethod h;
      h.heuristic = pump;
      h.id = "Pump for MINLP";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicMilpRounding = false;
    options()->GetEnumValue("MILP_rounding_heuristic",doHeuristicMilpRounding,prefix_.c_str());
    if(doHeuristicMilpRounding){
      Bonmin::MilpRounding * round = new Bonmin::MilpRounding(this);
      HeuristicMethod h;
      h.heuristic = round;
      h.id = "MILP Rounding";
      heuristics_.push_back(h);
    }
  }

  // Add Branching rules ///////////////////////////////////////////////////////

  switch (varSelection) {

  case OSI_STRONG: { // strong branching
    CouenneChooseStrong * chooseVariable = new CouenneChooseStrong
      (*this, couenneProb_, journalist ());
    chooseVariable->setTrustStrongForSolution(false);
    chooseVariable->setTrustStrongForBound(false);
    chooseVariable->setOnlyPseudoWhenTrusted(true);
    branchingMethod_ = chooseVariable;
    break;
  }

  case OSI_SIMPLE: // default choice
    branchingMethod_ = new CouenneChooseVariable 
      (continuousSolver_, couenneProb_, journalist ());
    break;

  default:
    std::cerr << "Unknown variable_selection for Couenne\n" << std::endl;
    throw;
    break;
  }

  // Node comparison method ///////////////////////////////////////////////////////////////////////////

  int ival;
  if (!options_->GetEnumValue("node_comparison", ival, "bonmin.")) {
    // change default for Couenne
    nodeComparisonMethod_ = bestBound;
  }
  else {
    nodeComparisonMethod_ = NodeComparison(ival);
  }

  if (intParam_[NumCutPasses] < 2)
    intParam_[NumCutPasses] = 2;

  // Tell Cbc not to check again if a solution returned from
  // heuristic is indeed feasible
  intParam_ [BabSetupBase::SpecialOption] = 16 | 4;

  // Add disjunctive cuts ///////////////////////////////////////////////////////

  options () -> GetIntegerValue ("minlp_disj_cuts", freq, "couenne.");

  if (freq != 0) {

    CouenneDisjCuts * couenneDisj = 
      new CouenneDisjCuts (ci, this, 
			   couenneCg, 
			   branchingMethod_, 
			   varSelection == OSI_STRONG, // if true, use strong branching candidates
			   journalist (),
			   options ());

    CuttingMethod cg;
    cg.frequency = freq;
    cg.cgl = couenneDisj;
    cg.id = "Couenne disjunctive cuts";
    cutGenerators (). push_back(cg);
  }

  // Add cross-aux redundant cuts ///////////////////////////////////////////////////////

  // options () -> GetIntegerValue ("crossconv_cuts", freq, "couenne.");

  // if (freq != 0) {

  //   CouenneCrossConv * couenneCross = 
  //     new CouenneCrossConv (couenneProb,
  // 			    journalist (),
  // 			    options ());

  //   CuttingMethod cg;
  //   cg.frequency = freq;
  //   cg.cgl = couenneCross;
  //   cg.id = "Couenne cross-aux cuts";
  //   cutGenerators (). push_back(cg);
  // }

  // Add sdp cuts ///////////////////////////////////////////////////////

  // options () -> GetIntegerValue ("sdp_cuts", freq, "couenne.");

  // if (freq != 0) {

  //   CouenneDisjCuts * couenneDisj = 
  //     new CouenneDisjCuts (ci, this, 
  // 			   couenneCg, 
  // 			   branchingMethod_, 
  // 			   varSelection == OSI_STRONG, // if true, use strong branching candidates
  // 			   journalist (),
  // 			   options ());

  //   CuttingMethod cg;
  //   cg.frequency = freq;
  //   cg.cgl = couenneDisj;
  //   cg.id = "Couenne disjunctive cuts";
  //   cutGenerators (). push_back(cg);
  // }

  return true;
}
 
void CouenneSetup::registerOptions ()
{registerAllOptions (roptions ());}


void CouenneSetup::registerAllOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> SetRegisteringCategory ("Couenne options", Bonmin::RegisteredOptions::CouenneCategory);

  BabSetupBase        ::registerAllOptions (roptions);
  Bonmin::BonCbcFullNodeInfo  ::registerOptions (roptions);

  /// Heuristics
  Bonmin::LocalSolverBasedHeuristic    ::registerOptions (roptions);
  Bonmin::FixAndSolveHeuristic         ::registerOptions (roptions);
  Bonmin::DummyPump                    ::registerOptions (roptions);
  Bonmin::MilpRounding                 ::registerOptions (roptions);
  Bonmin::PumpForMinlp                 ::registerOptions (roptions);
  Bonmin::HeuristicRINS                ::registerOptions (roptions);
  Bonmin::HeuristicLocalBranching      ::registerOptions (roptions);
  Bonmin::HeuristicFPump               ::registerOptions (roptions);
  Bonmin::HeuristicDiveFractional      ::registerOptions (roptions);
  Bonmin::HeuristicDiveVectorLength    ::registerOptions (roptions);
  Bonmin::HeuristicDiveMIPFractional   ::registerOptions (roptions);
  Bonmin::HeuristicDiveMIPVectorLength ::registerOptions (roptions);  

  roptions -> AddStringOption3 ("milp_solver",
				"Choose the subsolver to solve MILP sub-problems in OA decompositions.",
				"Cbc_D",
				"Cbc_D","Coin Branch and Cut with its default",
				"Cbc_Par", "Coin Branch and Cut with passed parameters",
				"Cplex","Ilog Cplex",
				" To use Cplex, a valid license is required and you should have compiled OsiCpx in COIN-OR  (see Osi documentation).");

  roptions -> setOptionExtraInfo ("milp_solver",64);

  roptions -> AddStringOption2 ("milp_strategy",
				"Choose a strategy for MILPs.",
				"find_good_sol",
				"find_good_sol","Stop sub milps when a solution improving the incumbent is found",
				"solve_to_optimality", "Solve MILPs to optimality",
				"");

  roptions -> AddStringOption6 ("algorithm",
				"Choice of the algorithm.",
				"B-BB",
				"B-BB","simple branch-and-bound algorithm,",
				"B-OA","OA Decomposition algorithm,",
				"B-QG","Quesada and Grossmann branch-and-cut algorithm,",
				"B-Hyb","hybrid outer approximation based branch-and-cut,",
				"B-Ecp","ecp cuts based branch-and-cut a la FilMINT.",
				"B-iFP","Iterated Feasibility Pump for MINLP.",
				"This will preset some of the options of bonmin depending on the algorithm choice."
				);

  CouenneProblem          ::registerOptions (roptions);
  CouenneCutGenerator     ::registerOptions (roptions);
  CouenneChooseStrong     ::registerOptions (roptions);
  CouenneChooseVariable   ::registerOptions (roptions);
  CouenneFixPoint         ::registerOptions (roptions);
  CouenneDisjCuts         ::registerOptions (roptions);
  //  CouenneCrossConv        ::registerOptions (roptions);
  CouenneTwoImplied       ::registerOptions (roptions);
  NlpSolveHeuristic       ::registerOptions (roptions);
  CouenneFeasPump         ::registerOptions (roptions);
  CouenneIterativeRounding::registerOptions (roptions);

  /// TODO: move later!
  roptions -> AddStringOption2
    ("local_branching_heuristic",
     "Apply local branching heuristic",
     "no",
     "no","",
     "yes","",
     "A local-branching heuristic based is used to find feasible solutions.");


  roptions -> AddNumberOption  ("couenne_check",
				"known value of a global optimum (for debug purposes only)",
				COIN_DBL_MAX,
				"Default value is +infinity.");

  roptions -> AddStringOption2 ("display_stats",
				"display statistics at the end of the run",
				"no",
				"yes", "",
				"no", "");

  roptions -> AddStringOption2 ("test_mode",
				"set to true if this is Couenne unit test",
				"no",
				"yes", "",
				"no", "");

  roptions -> AddStringOption4 ("lp_solver",
				"Linear Programming solver for the linearization",
				"clp",
				"clp",    "Use the Coin-OR Open Source solver CLP",
				"cplex",  "Use the commercial solver Cplex (license is needed)",
				"gurobi", "Use the commercial solver Gurobi (license is needed)",
				"soplex", "Use the freely available Soplex");

#define addLevOption(optname,comment) roptions -> AddBoundedIntegerOption (optname, comment, -2, J_LAST_LEVEL-1, J_NONE, "")

  addLevOption ("output_level",                "Output level");
  addLevOption ("branching_print_level",       "Output level for braching code in Couenne");
  addLevOption ("boundtightening_print_level", "Output level for bound tightening code in Couenne");
  addLevOption ("convexifying_print_level",    "Output level for convexifying code in Couenne");
  addLevOption ("problem_print_level",         "Output level for problem manipulation code in Couenne");
  addLevOption ("nlpheur_print_level",         "Output level for NLP heuristic in Couenne");
  addLevOption ("disjcuts_print_level",        "Output level for disjunctive cuts in Couenne");
  addLevOption ("reformulate_print_level",     "Output level for reformulating problems in Couenne");

  roptions -> AddNumberOption
    ("feas_tolerance",
     "Tolerance for constraints/auxiliary variables",
     feas_tolerance_default,
     "Default value is 1e-5.");

  roptions -> AddStringOption2 
    ("feasibility_bt",
     "Feasibility-based (cheap) bound tightening (FBBT)",
     "yes",
     "no","",
     "yes","",
     "A pre-processing technique to reduce the bounding box, before the generation of linearization cuts. "
     "This is a quick and effective way to reduce the solution set, and it is highly recommended to keep it active."
    );

  // copied from BonminSetup::registerMilpCutGenerators(), in
  // BonBonminSetup.cpp

  struct cutOption_ {

    const char *cgname;
    int         defaultFreq;

  } cutOption [] = {
    {(const char *) "Gomory_cuts",             0},
    {(const char *) "probing_cuts",            0},
    {(const char *) "cover_cuts",              0},
    {(const char *) "mir_cuts",                0},
    {(const char *) "2mir_cuts",               0},
    {(const char *) "flow_covers_cuts",        0},
    {(const char *) "lift_and_project_cuts",   0},
    {(const char *) "reduce_split_cuts",       0},
    {(const char *) "clique_cuts",             0},
    {NULL, 0}};

  for (int i=0; cutOption [i].cgname; i++) {

    char descr [150];

    sprintf (descr, "Frequency k (in terms of nodes) for generating %s cuts in branch-and-cut.",
	     cutOption [i].cgname);

    roptions -> AddLowerBoundedIntegerOption 
      (cutOption [i].cgname,
       descr,
       -100, cutOption [i].defaultFreq,
       "If k > 0, cuts are generated every k nodes, "
       "if -99 < k < 0 cuts are generated every -k nodes but "
       "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
       "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");

    roptions->setOptionExtraInfo (cutOption [i].cgname, 5);
  }
}



/** Add milp cut generators according to options.*/
void CouenneSetup::addMilpCutGenerators () {

  enum extraInfo_ {CUTINFO_NONE, CUTINFO_MIG, CUTINFO_PROBING, CUTINFO_CLIQUE};

  // extra data structure to avoid repeated code below

  struct cutInfo {

    const char      *optname;
    CglCutGenerator *cglptr;
    const char      *cglId;
    enum extraInfo_  extraInfo;

  } cutList [] = {
    {(const char*)"Gomory_cuts",new CglGomory,      (const char*)"Mixed Integer Gomory",CUTINFO_MIG},
    {(const char*)"probing_cuts",new CglProbing,        (const char*) "Probing", CUTINFO_PROBING},
    {(const char*)"mir_cuts",new CglMixedIntegerRounding2, (const char*) "Mixed Integer Rounding", 
     CUTINFO_NONE},
    {(const char*)"2mir_cuts",    new CglTwomir,         (const char*) "2-MIR",    CUTINFO_NONE},
    {(const char*)"cover_cuts",   new CglKnapsackCover,  (const char*) "Cover",    CUTINFO_NONE},
    {(const char*)"clique_cuts",  new CglClique,         (const char*) "Clique",   CUTINFO_CLIQUE},
    {(const char*)"lift_and_project_cuts",new CglLandP,(const char*)"Lift and Project",CUTINFO_NONE},
    {(const char*)"reduce_split_cuts",new CglRedSplit,(const char*) "Reduce and Split",CUTINFO_NONE},
    {(const char*)"flow_covers_cuts",new CglFlowCover,(const char*) "Flow cover cuts", CUTINFO_NONE},
    {NULL, NULL, NULL, CUTINFO_NONE}};

  int freq;

  for (int i=0; cutList [i]. optname; i++) {

    options_ -> GetIntegerValue (std::string (cutList [i]. optname), freq, "couenne.");

    if (!freq) {
      delete cutList [i].cglptr;
      continue;
    }

    CuttingMethod cg;
    cg.frequency = freq;
    cg.cgl       = cutList [i].cglptr;
    cg.id        = std::string (cutList [i]. cglId);
    cutGenerators_.push_back (cg);

    // options for particular cases
    switch (cutList [i].extraInfo) {

    case CUTINFO_MIG: {
      CglGomory *gc = dynamic_cast <CglGomory *> (cutList [i].cglptr);

      if (!gc) break;

      gc -> setLimitAtRoot(512);
      gc -> setLimit(50);
    }
      break;

    case CUTINFO_PROBING: {
      CglProbing *pc = dynamic_cast <CglProbing *> (cutList [i].cglptr);

      if (!pc) break;

      pc->setUsingObjective(1);
      pc->setMaxPass(3);
      pc->setMaxPassRoot(3);
      // Number of unsatisfied variables to look at
      pc->setMaxProbe(10);
      pc->setMaxProbeRoot(50);
      // How far to follow the consequences
      pc->setMaxLook(10);
      pc->setMaxLookRoot(50);
      pc->setMaxLookRoot(10);
      // Only look at rows with fewer than this number of elements
      pc->setMaxElements(200);
      pc->setRowCuts(3);
    }
      break;

    case CUTINFO_CLIQUE: {
      CglClique *clique = dynamic_cast <CglClique *> (cutList [i].cglptr);

      if (!clique) break;

      clique -> setStarCliqueReport(false);
      clique -> setRowCliqueReport(false);
      clique -> setMinViolation(0.1);
    }
      break;

      //case CUTINFO_NONE:
    default:
      break;
    }
  }

  double givenAllowFGap2 = 0.0;
  options_->GetNumericValue(std::string("allowable_fraction_gap"), 
			    givenAllowFGap2, "bonmin.");
  double upval = 1e50;

#ifdef FM_UP_BND
  printf("CutOff value:\n");
  scanf("%lf", &upval);
#else /* not FM_UP_BND */
  options_->GetNumericValue(std::string("art_cutoff"), upval, "bonmin.");
#endif /* FM_UP_BND */

  if(upval < 1e50) {
    double newCO = (1-givenAllowFGap2) * upval;
    couenneProb_->setCutOff(newCO);
    printf("CutOff set to %f\n", newCO);

#ifdef FM_TRACE_OPTSOL
    if(couenneProb_->getRecordBestSol()->getHasSol()) {
      if(newCO < couenneProb_->getRecordBestSol()->getVal()) {
	couenneProb_->getRecordBestSol()->setVal(newCO);
      }
    }
#endif
  }
}
