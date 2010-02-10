// $Id$
//
// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007

#include "BonCouenneSetup.hpp"
#include "BonInitHeuristic.hpp"
#include "BonNlpHeuristic.hpp"
#include "BonCouenneInterface.hpp"

#include "BonGuessHeuristic.hpp"
#include "CbcCompareActual.hpp"

#include "CouenneObject.hpp"
//#include "CouenneComplObject.hpp"
#include "CouenneVarObject.hpp"
#include "CouenneVTObject.hpp"
#include "CouenneOrbitObj.hpp"
#include "CouenneChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneSolverInterface.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneDisjCuts.hpp"

#include "BonCouenneInfo.hpp"
#include "BonCbcNode.hpp"

#include "OsiClpSolverInterface.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
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

// Ampl includes
#ifdef COIN_HAS_ASL
#include "asl.h"
#include "getstub.h"
#endif


namespace Bonmin{
  
  SmartAsl::~SmartAsl(){
#ifdef COIN_HAS_ASL
    //Code from Ipopt::AmplTNLP to free asl
    if(asl != NULL){
        if (X0) {
          delete [] X0;
          X0 = NULL;
        }
        if (havex0) {
          delete [] havex0;
          havex0 = NULL;
        }
        if (pi0) {
          delete [] pi0;
          pi0 = NULL;
        }
        if (havepi0) {
          delete [] havepi0;
          havepi0 = NULL;
        }
        ASL* asl_to_free = (ASL*)asl;
        ASL_free(&asl_to_free);
        asl = NULL;
    }
    ASL_free(&asl);
#endif
  }
  
  CouenneSetup::~CouenneSetup(){
    if (couenneProb_ && couenneProb_is_own_)
      delete couenneProb_;
  }

  bool CouenneSetup::InitializeCouenne (char ** argv, 
					CouenneProblem *couenneProb,
          Ipopt::SmartPtr<Bonmin::TMINLP> tminlp,
					Bonmin::CouenneInterface *ci) {
    std::string s;

    if (couenneProb) {
    	//TODO create a copy of user problem, since we modify it? 
    	couenneProb_ = couenneProb;
    	couenneProb_is_own_ = false;
    }

    /* Get the basic options. */
    readOptionsFile();
 
    // in check mode, avoid pop-up error message (there are quite a few messages)
    options_ -> GetStringValue ("test_mode", s, "couenne.");
    if (s == "yes")
      WindowsErrorPopupBlocker();

    /** Change default value for failure behavior so that code doesn't crash 
	when Ipopt does not solve a sub-problem.*/

    options_ -> SetStringValue ("nlp_failure_behavior", "fathom", "bonmin.");

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
      	assert(couenneProb_ != NULL);
      	assert(IsValid(tminlp)); //TODO would be great to setup own TMINLP based on CouenneProblem formulation
      	ci -> initialize(roptions_, options_, journalist_, tminlp);
      }
    }

    nonlinearSolver_ = ci;

    /** Set the output level for the journalist for all Couenne
     categories.  We probably want to make that a bit more flexible
     later. */
    int i;

    options()->GetIntegerValue("boundtightening_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_BOUNDTIGHTENING, (EJournalLevel) i);

    options()->GetIntegerValue("branching_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_BRANCHING, (EJournalLevel) i);

    options()->GetIntegerValue("convexifying_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_CONVEXIFYING, (EJournalLevel) i);

    options()->GetIntegerValue("problem_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_PROBLEM, (EJournalLevel) i);

    options()->GetIntegerValue("nlpheur_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_NLPHEURISTIC, (EJournalLevel) i);

    options()->GetIntegerValue("disjcuts_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_DISJCUTS, (EJournalLevel) i);

    options()->GetIntegerValue("reformulate_print_level", i, "bonmin.");
    journalist()->GetJournal("console")-> SetPrintLevel(J_REFORMULATE, (EJournalLevel) i);

    /* Initialize Couenne cut generator.*/
    //int ivalue, num_points;
    //options()->GetEnumValue("convexification_type", ivalue,"bonmin.");
    //options()->GetIntegerValue("convexification_points",num_points,"bonmin.");

    if (!couenneProb_)
      couenneProb_ = new CouenneProblem (aslfg_ -> asl, this, journalist ());

    CouenneCutGenerator * couenneCg = 
      new CouenneCutGenerator (ci, this, couenneProb_, NULL);

    options_ -> GetStringValue ("lp_solver", s, "couenne.");

    if (s == "clp") {

      CouenneSolverInterface <OsiClpSolverInterface> *CSI 
	= new CouenneSolverInterface <OsiClpSolverInterface>;

      continuousSolver_ = CSI;
      CSI -> setCutGenPtr (couenneCg);

    } else if (s == "cplex") {

#ifdef COIN_HAS_CPX
      CouenneSolverInterface <OsiCpxSolverInterface> *CSI 
	= new CouenneSolverInterface <OsiCpxSolverInterface>;

      continuousSolver_ = CSI;
      CSI -> setCutGenPtr (couenneCg);
#else
      journalist()->Printf(J_ERROR, J_INITIALIZATION, "Couenne was compiled without CPLEX interface. Please reconfigure, recompile, and try again.\n");
      return false;
#endif
    }
    continuousSolver_ -> passInMessageHandler(ci -> messageHandler());

    couenneProb_ -> setBase (this);

    assert (couenneProb_);

    couenneProb_ -> reformulate (couenneCg);

    Bonmin::BabInfo * extraStuff = new Bonmin::CouenneInfo(0);

    // as per instructions by John Forrest, to get changed bounds
    extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);

    continuousSolver_ -> setAuxiliaryInfo (extraStuff);
    delete extraStuff;
    
    extraStuff = dynamic_cast <Bonmin::BabInfo *> (continuousSolver_ -> getAuxiliaryInfo ());

    /* Setup log level*/
    int lpLogLevel;
    options()->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);

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
      h.id = "Init Rounding NLP";
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

      nSOS = couenneProb_ -> findSOS (nonlinearSolver (), objects);

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

    int contObjPriority = 2000; // default object priority -- it is 1000 for integers and 10 for SOS

    options () -> GetIntegerValue ("cont_var_priority", contObjPriority, "bonmin.");

    for (int i = 0; i < nVars; i++) { // for each variable

      exprVar *var = couenneProb_ -> Var (i);

      // we only want enabled variables
      if (var -> Multiplicity () <= 0) 
	continue;

      switch (objType) {

      case CouenneObject::EXPR_OBJ:

	// if this variable is associated with a nonlinear function
	if (var -> isInteger () || 
	    (var -> Type  () == AUX) && 
	    (var -> Image () -> Linearity () > LINEAR)) {

	  /*if ((var -> Image () -> code () == COU_EXPRMUL) &&
	      (var -> Image () -> ArgList () [0] -> Index () >= 0) &&
	      (var -> Image () -> ArgList () [1] -> Index () >= 0) &&
	      (fabs (var -> lb ()) < COUENNE_EPS) &&
	      (fabs (var -> ub ()) < COUENNE_EPS))

	    // it's a complementarity constraint object!
	    objects    [nobj] = new CouenneComplObject (couenneProb_, var, this, journalist ());
	    else*/
	  objects [nobj] = new CouenneObject (couenneCg, couenneProb_, var, this, journalist ());

	  objects [nobj++] -> setPriority (contObjPriority);
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

	  objects [nobj] = new CouenneVarObject (couenneCg, couenneProb_, var, this, journalist ());
	  objects [nobj++] -> setPriority (contObjPriority);
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
	  objects [nobj++] -> setPriority (contObjPriority);
	  //objects [nobj++] -> setPriority (contObjPriority + var -> rank ());
	}

	break;
      }
    }

    // Experimental: orbital branching //////////////////////////////////////////////
    options () -> GetStringValue ("orbital_branching", s, "couenne.");

    if (s == "yes") {

      objects [nobj] = new CouenneOrbitObj (couenneCg, couenneProb_, NULL, this, journalist ());
      objects [nobj++] -> setPriority (contObjPriority);
    }

    // Add objects /////////////////////////////////

    continuousSolver_ -> addObjects (nobj, objects);

    for (int i = 0 ; i < nobj ; i++)
      delete objects [i];

    delete [] objects;

    // Setup Convexifier generators ////////////////////////////////////////////////

    int freq;

    options()->GetIntegerValue("convexification_cuts",freq,"couenne.");

    if (freq != 0) {

      CuttingMethod cg;
      cg.frequency = freq;
      cg.cgl = couenneCg;
      cg.id = "Couenne convexifier cuts";
      cutGenerators().push_back(cg);

      // set cut gen pointer
      //dynamic_cast <CouenneSolverInterface <OsiClpSolverInterface> *> 
      //(continuousSolver_)

      // this is done on an explicitly declared CSI pointer, however
      // CSI == continuousSolver_
    }

    // disjunctive cuts generator added AFTER 

    // add other cut generators -- test for integer variables first
    if (couenneCg -> Problem () -> nIntVars () > 0)
      addMilpCutGenerators ();

    CouennePtr_ = couenneCg;

    // Setup heuristic to solve nlp problems. /////////////////////////////////

    int doNlpHeurisitic = 0;
    options()->GetEnumValue("local_optimization_heuristic", doNlpHeurisitic, "couenne.");
    if(doNlpHeurisitic)
    {
      int numSolve;
      options()->GetIntegerValue("log_num_local_optimization_per_level",numSolve,"couenne.");
      NlpSolveHeuristic * nlpHeuristic = new NlpSolveHeuristic;
      nlpHeuristic->setNlp(*ci,false);
      nlpHeuristic->setCouenneProblem(couenneProb_);
      //nlpHeuristic->setMaxNlpInf(1e-4);
      nlpHeuristic->setMaxNlpInf(maxNlpInf_0);
      nlpHeuristic->setNumberSolvePerLevel(numSolve);
      HeuristicMethod h;
      h.id = "Couenne Rounding NLP";
      h.heuristic = nlpHeuristic;
      heuristics_.push_back(h);
    }

    // Add Branching rules ///////////////////////////////////////////////////////

    int varSelection;
    if (!options_->GetEnumValue("variable_selection",varSelection,"bonmin.")) {
      // change the default for Couenne
      varSelection = OSI_SIMPLE;
    }

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

    int ival;
    if (!options_->GetEnumValue("node_comparison",ival,"bonmin.")) {
      // change default for Couenne
      nodeComparisonMethod_ = bestBound;
    }
    else {
      nodeComparisonMethod_ = NodeComparison(ival);
    }

    if(intParam_[NumCutPasses] < 2)
    intParam_[NumCutPasses] = 2;

    // Tell Cbc not to check again if a solution returned from
    // heuristic is indeed feasible
    intParam_[BabSetupBase::SpecialOption] = 16 | 4;

    return true;
  }

  bool CouenneSetup::InitializeCouenne (char ** argv,
          CouenneProblem *couenneProb,
          Bonmin::CouenneInterface *ci) {
    return InitializeCouenne(argv, couenneProb, ci ? ci->model() : NULL, ci);
  }
 
  void CouenneSetup::registerOptions(){
    registerAllOptions(roptions());
  }


  void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    BonCbcFullNodeInfo::registerOptions(roptions);
    CouenneCutGenerator::registerOptions (roptions);
    CouenneDisjCuts::registerOptions (roptions);

    roptions -> AddNumberOption
      ("couenne_check",
       "known value of a global optimum",
       COIN_DBL_MAX,
       "Default value is +infinity.");

    roptions -> AddStringOption2 (
      "display_stats",
      "display statistics at the end of the run",
      "no",
      "yes", "",
      "no", "");

    roptions -> AddStringOption2 (
      "test_mode",
      "set to true if this is Couenne unit test",
      "no",
      "yes", "",
      "no", "");

    roptions -> AddStringOption3 (
      "lp_solver",
      "Linear Programming solver for the linearization",
      "clp",
      "clp", "Use the Coin-OR Open Source solver CLP",
      "cplex", "Use the commercial solver Cplex (license is needed)",
      "soplex", "Use the freely available Soplex (not available yet)");

    roptions->AddBoundedIntegerOption(
      "branching_print_level", "Output level for braching code in Couenne",
      -2, J_LAST_LEVEL-1, J_NONE, "");

    roptions->AddBoundedIntegerOption(
      "boundtightening_print_level", "Output level for bound tightening code in Couenne",
      -2, J_LAST_LEVEL-1, J_NONE, "");

    roptions->AddBoundedIntegerOption(
      "convexifying_print_level", "Output level for convexifying code in Couenne",
      -2, J_LAST_LEVEL-1, J_NONE, "");

    roptions->AddBoundedIntegerOption(
      "problem_print_level", "Output level for problem manipulation code in Couenne",
      -2, J_LAST_LEVEL-1, J_ERROR, "");

    roptions->AddBoundedIntegerOption(
      "nlpheur_print_level", "Output level for NLP heuristic in Couenne",
      -2, J_LAST_LEVEL-1, J_NONE, "");

    roptions->AddBoundedIntegerOption(
    "disjcuts_print_level", "Output level for disjunctive cuts in Couenne",
    -2, J_LAST_LEVEL-1, J_NONE, "");

    roptions->AddBoundedIntegerOption(
    "reformulate_print_level", "Output level for reformulating problems in Couenne",
    -2, J_LAST_LEVEL-1, J_NONE, "");


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

      options_ -> GetIntegerValue (std::string (cutList [i]. optname), freq, "bonmin.");

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
  }
}
