/* $Id$
 *
 * Name:    CouenneProblemConstructors.cpp
 * Author:  Pietro Belotti
 * Purpose: Constructors and destructors of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2009-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "BonBabSetupBase.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprQuad.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprIVar.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprOpp.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneGlobalCutOff.hpp"
#include "CouenneDepGraph.hpp"
#include "CouenneLQelems.hpp"

#include "CouenneObject.hpp"

#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

#define MAX_FBBT_ITER 3

/// constructor
CouenneProblem::CouenneProblem (struct ASL *asl,
				Bonmin::BabSetupBase *base,
				JnlstPtr jnlst):
  problemName_ (""),
  auxSet_    (NULL), 
  curnvars_  (-1),
  nIntVars_  (0),
  optimum_   (NULL),
  bestObj_   (COIN_DBL_MAX),
  quadIndex_ (NULL),
  commuted_  (NULL),
  numbering_ (NULL),
  ndefined_  (0),
  graph_     (NULL),
  nOrigVars_ (0),
  nOrigIntVars_ (0),
  pcutoff_   (new GlobalCutOff (COIN_DBL_MAX)),
  created_pcutoff_ (true),
  doFBBT_    (true),
  doRCBT_    (true),
  doOBBT_    (true),
  doABT_     (true),
  logObbtLev_(0),
  logAbtLev_ (0),
  jnlst_     (jnlst),
  opt_window_ (COIN_DBL_MAX),
  useQuadratic_ (false),
  feas_tolerance_ (feas_tolerance_default),
  integerRank_ (NULL),
  maxCpuTime_  (COIN_DBL_MAX),
  bonBase_     (base),
#ifdef COIN_HAS_ASL
  asl_         (asl),
#endif
  unusedOriginalsIndices_ (NULL),
  nUnusedOriginals_ (-1),
  multilinSep_ (CouenneProblem::MulSepNone),
  max_fbbt_iter_ (MAX_FBBT_ITER),
  orbitalBranching_ (false) {

  double now = CoinCpuTime ();

  if (asl) {
#if COIN_HAS_ASL
    // read problem from AMPL structure
    readnl (asl);
#else
    jnlst_ -> Printf (Ipopt::J_ERROR, J_PROBLEM, "Couenne was compiled without the ASL library. Cannot process ASL structure.\n");
    throw -1;
#endif

    if ((now = (CoinCpuTime () - now)) > 10.)
      jnlst_ -> Printf (Ipopt::J_WARNING, J_PROBLEM,
			"Couenne: reading time %.3fs\n", now);
  }

  // create expression set for binary search
  auxSet_ = new std::set <exprAux *, compExpr>;

  if (base)
    initOptions (base -> options());

  recBSol = new struct Couenne::CouenneRecordBestSol();
  lastPrioSort_ = 1000000;

  minDepthPrint_ = -1;
  minNodePrint_ = -1;
  doPrint_ = false;
}


/// copy constructor

CouenneProblem::CouenneProblem (const CouenneProblem &p):
  problemName_  (p.problemName_),
  domain_       (p.domain_),
  curnvars_     (-1),
  nIntVars_     (p.nIntVars_),
  optimum_      (NULL),
  bestObj_      (p.bestObj_),
  commuted_     (NULL),
  numbering_    (NULL),
  ndefined_     (p.ndefined_),
  graph_        (NULL),
  nOrigVars_    (p.nOrigVars_),
  nOrigCons_    (p.nOrigCons_),
  nOrigIntVars_ (p.nOrigIntVars_),
  pcutoff_      (p.pcutoff_),
  created_pcutoff_ (false),
  doFBBT_       (p. doFBBT_),
  doRCBT_       (p. doRCBT_),
  doOBBT_       (p. doOBBT_),
  doABT_        (p. doABT_),
  logObbtLev_   (p. logObbtLev_),
  logAbtLev_    (p. logAbtLev_),
  jnlst_        (p.jnlst_),
  opt_window_   (p.opt_window_),    // needed only in standardize (), unnecessary to update it
  useQuadratic_ (p.useQuadratic_),  // ditto
  feas_tolerance_ (p.feas_tolerance_),
  dependence_   (p.dependence_),
  objects_      (p.objects_), // NO! have to copy all of them 
  integerRank_  (NULL),
  numberInRank_ (p.numberInRank_),
  maxCpuTime_   (p.maxCpuTime_),
  bonBase_      (p.bonBase_),
#ifdef COIN_HAS_ASL
  asl_          (p.asl_),
#endif
  unusedOriginalsIndices_ (NULL),
  nUnusedOriginals_ (p.nUnusedOriginals_),
  multilinSep_  (p.multilinSep_),
  max_fbbt_iter_  (p.max_fbbt_iter_),
  orbitalBranching_  (p.orbitalBranching_) {

  for (int i=0; i < p.nVars (); i++)
    variables_ . push_back (NULL);

  for (int i=0; i < p.nVars (); i++) {
    int ind = p.numbering_ [i];
    variables_ [ind] = p.Var (ind) -> clone (&domain_);
  }

  for (std::vector <CouenneObject *>::iterator i = objects_.begin ();
       i != objects_.end (); ++i)
    (*i) = (*i) -> clone ();

  if (p.numbering_)
    numbering_ = CoinCopyOfArray (p.numbering_, nVars ());

  // clone objectives and constraints (there's a leak around here)
  for (int i=0; i < p.nObjs (); i++) objectives_  . push_back (p.Obj (i) -> clone (&domain_));
  for (int i=0; i < p.nCons (); i++) constraints_ . push_back (p.Con (i) -> clone (&domain_));

  if (p.optimum_) 
    optimum_ = CoinCopyOfArray (p.optimum_, nVars ());
    
  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();

  // copy integer rank (used in getIntegerCandidate)
  if (p.integerRank_) {
    integerRank_ = new int [nVars ()];
    CoinCopyN (p.integerRank_, nVars (), integerRank_);
  }

  // copy unusedOriginals
  if (nUnusedOriginals_ > 0) {
    unusedOriginalsIndices_ = (int *) malloc (nUnusedOriginals_ * sizeof (int));
    CoinCopyN (p.unusedOriginalsIndices_, nUnusedOriginals_, unusedOriginalsIndices_);
  }

  recBSol = new CouenneRecordBestSol(*(p.recBSol));
  lastPrioSort_ = p.lastPrioSort_;

  minDepthPrint_ = p.minDepthPrint_;
  minNodePrint_ = p.minNodePrint_;
  doPrint_ = p.doPrint_;
}


/// Destructor

CouenneProblem::~CouenneProblem () {

  // delete optimal solution (if any)
  if (optimum_)
    free (optimum_);

  // delete objectives
  for (std::vector <CouenneObjective *>::iterator i  = objectives_ . begin ();
       i != objectives_ . end (); ++i)
    delete (*i);

  // delete constraints
  for (std::vector <CouenneConstraint *>::iterator i = constraints_ . begin ();
       i != constraints_ . end (); ++i)
    delete (*i);

  // deletion of variables should be done in reverse order w.r.t. the
  // dependence. If numbering_ is available, use it.
  if (numbering_) for (int i=nVars (); i--;) delete variables_ [numbering_ [i]];
  else            for (int i=nVars (); i--;) delete variables_ [i];

  // delete extra structures
  if (graph_)     delete    graph_;
  if (commuted_)  delete [] commuted_;
  if (numbering_) delete [] numbering_;

  if (created_pcutoff_) delete pcutoff_;

  if (integerRank_) delete [] integerRank_;

  if (unusedOriginalsIndices_)
    free (unusedOriginalsIndices_);

  for (std::vector <CouenneObject *>::iterator i = objects_.begin ();
       i != objects_.end (); ++i)
    delete (*i);

  delete recBSol;
}


/// initializes parameters like doOBBT
void CouenneProblem::initOptions(Ipopt::SmartPtr<Ipopt::OptionsList> options) {

  assert(IsValid(options));

  std::string s;

  options -> GetStringValue ("use_quadratic",   s, "couenne."); useQuadratic_ = (s == "yes");
  options -> GetStringValue ("feasibility_bt",  s, "couenne."); doFBBT_ = (s == "yes");
  options -> GetStringValue ("redcost_bt",      s, "couenne."); doRCBT_ = (s == "yes");
  options -> GetStringValue ("optimality_bt",   s, "couenne."); doOBBT_ = (s == "yes");
  options -> GetStringValue ("aggressive_fbbt", s, "couenne."); doABT_  = (s == "yes");

  options -> GetIntegerValue ("log_num_obbt_per_level", logObbtLev_, "couenne.");
  options -> GetIntegerValue ("log_num_abt_per_level",  logAbtLev_,  "couenne.");

  options -> GetIntegerValue ("max_fbbt_iter",  max_fbbt_iter_,  "couenne.");

  options -> GetNumericValue ("feas_tolerance",  feas_tolerance_, "couenne.");
  options -> GetNumericValue ("opt_window",      opt_window_,     "couenne.");

  options -> GetStringValue ("multilinear_separation", s, "couenne.");
  multilinSep_ = (s == "none"   ? CouenneProblem::MulSepNone   :
		  s == "simple" ? CouenneProblem::MulSepSimple :
                 		  CouenneProblem::MulSepTight);

  options -> GetStringValue ("orbital_branching", s, "couenne."); orbitalBranching_ = (s == "yes");

  options -> GetStringValue ("quadrilinear_decomp",  s, "couenne."); 

  if      (s == "rAI")     trilinDecompType_ = rAI;
  else if (s == "tri+bi")  trilinDecompType_ = tri_bi;
  else if (s == "bi+tri")  trilinDecompType_ = bi_tri;
  else if (s == "hier-bi") trilinDecompType_ = treeDecomp;
}
