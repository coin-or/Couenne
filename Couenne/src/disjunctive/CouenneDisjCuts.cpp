/* $Id$
 *
 * Name:    CouenneDisjCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for the disjunctive cuts
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"


/// constructor
CouenneDisjCuts::CouenneDisjCuts (Bonmin::OsiTMINLPInterface *minlp,
				  Bonmin::BabSetupBase *base,
				  CouenneCutGenerator *cg,
				  OsiChooseVariable *bcv,
				  bool is_strong,
				  JnlstPtr journalist,
				  const Ipopt::SmartPtr<Ipopt::OptionsList> options):
  couenneCG_          (cg),
  nrootcuts_          (-1), // to indicate first iteration not done yet
  ntotalcuts_         (0),
  septime_            (0.),
  objValue_           (-COIN_DBL_MAX),
  minlp_              (minlp),
  branchingMethod_    (bcv),
  isBranchingStrong_  (is_strong),
  jnlst_              (journalist),
  activeRows_         (false),
  activeCols_         (false),
  addPreviousCut_     (false),
  cpuTime_            (-1.) {

  options -> GetNumericValue ("time_limit", cpuTime_,  "couenne.");

  options -> GetNumericValue ("disj_init_perc",   initDisjPercentage_,  "couenne.");
  options -> GetIntegerValue ("disj_init_number", initDisjNumber_,      "couenne.");
  options -> GetIntegerValue ("disj_depth_level", depthLevelling_,      "couenne.");
  options -> GetIntegerValue ("disj_depth_stop",  depthStopSeparate_,   "couenne.");

  std::string s;
  options -> GetStringValue ("disj_active_rows", s, "couenne."); activeRows_     = (s == "yes");
  options -> GetStringValue ("disj_active_cols", s, "couenne."); activeCols_     = (s == "yes");
  options -> GetStringValue ("disj_cumulative",  s, "couenne."); addPreviousCut_ = (s == "yes");
}


/// copy constructorProduceOutput
CouenneDisjCuts::CouenneDisjCuts (const CouenneDisjCuts &src):
  couenneCG_          (src.couenneCG_),
  nrootcuts_          (src.nrootcuts_),
  ntotalcuts_         (src.ntotalcuts_),
  septime_            (src.septime_),
  objValue_           (src.objValue_),
  minlp_              (src.minlp_),
  branchingMethod_    (src.branchingMethod_),
  isBranchingStrong_  (src.isBranchingStrong_),
  jnlst_              (src.jnlst_),
  initDisjPercentage_ (src.initDisjPercentage_),
  initDisjNumber_     (src.initDisjNumber_),
  depthLevelling_     (src.depthLevelling_),
  depthStopSeparate_  (src.depthStopSeparate_),
  activeRows_         (src.activeRows_),
  activeCols_         (src.activeCols_),
  addPreviousCut_     (src.addPreviousCut_),
  cpuTime_            (src.cpuTime_) {}



/// destructor
CouenneDisjCuts::~CouenneDisjCuts ()
{if (septime_ > 1e-9) jnlst_ -> Printf (J_ERROR, J_DISJCUTS, "Disjunctive cuts: total time %g\n", septime_);}


/// Add list of options to be read from file
void CouenneDisjCuts::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("minlp_disj_cuts",
     "The frequency (in terms of nodes) at which Couenne disjunctive cuts are generated.",
     -99, 0,
     "A frequency of 0 (default) means these cuts are never generated. "
     "Any positive number n instructs Couenne to generate them at every n nodes of the B&B tree. "
     "A negative number -n means that generation should be attempted at the root node, and if successful it can be repeated at every n nodes, otherwise it is stopped altogether."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("disj_init_number",
     "Maximum number of disjunction to consider at each iteration.",
     -1, 10, "-1 means no limit.");

  roptions -> AddBoundedNumberOption
    ("disj_init_perc",
     "The maximum fraction of all disjunctions currently violated by the problem to consider for generating disjunctions.",
     0., false,
     1., false,
     0.5, "");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_level",
     "Depth of the B&B tree when to start decreasing the number of objects that generate disjunctions.",
     -1, 5, "This has a similar behavior as log_num_obbt_per_level. "
     "A value of -1 means that generation can be done at all nodes.");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_stop",
     "Depth of the B&B tree where separation of disjunctive cuts is stopped.",
     -1, 20, "A value of -1 means that generation can be done at all nodes");

  roptions -> AddStringOption2
    ("disj_active_rows",
     "Only include violated linear inequalities in the CGLP.",
     "no", 
     "yes", "",
     "no", "",
     "This reduces the size of the CGLP, but may produce less efficient cuts.");

  roptions -> AddStringOption2
    ("disj_active_cols",
     "Only include violated variable bounds in the Cut Generating LP (CGLP).",
     "no", 
     "yes", "",
     "no", "",
     "This reduces the size of the CGLP, but may produce less efficient cuts."
     );

  roptions -> AddStringOption2
    ("disj_cumulative",
     "Add previous disjunctive cut to current CGLP.",
     "no", 
     "yes", "",
     "no", "",
     "When generating disjunctive cuts on a set of disjunctions 1, 2, ..., k, introduce the cut relative to the previous disjunction i-1 in the CGLP used for disjunction i. "
     "Notice that, although this makes the cut generated more efficient, it increases the rank of the disjunctive cut generated."
    );
}
