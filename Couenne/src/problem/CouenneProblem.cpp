/*
 * Name:    CouenneProblem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CouenneTypes.hpp"

#include "expression.hpp"
#include "exprConst.hpp"
#include "exprQuad.hpp"
#include "exprClone.hpp"
#include "exprIVar.hpp"
#include "exprAux.hpp"
#include "exprOpp.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "depGraph.hpp"
#include "lqelems.hpp"

/// constructor
CouenneProblem::CouenneProblem (const struct ASL *asl,
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
  doOBBT_    (false),
  doABT_     (true),
  logObbtLev_(0),
  logAbtLev_ (0),
  jnlst_     (jnlst),
  opt_window_ (COIN_DBL_MAX),
  useQuadratic_ (false),
  feas_tolerance_ (feas_tolerance_default),
  integerRank_ (NULL),
  maxCpuTime_  (COIN_DBL_MAX),
  bonBase_     (base)
{

  double now = CoinCpuTime ();

  if (asl) {

    // read problem from AMPL structure
    readnl (asl);

    if ((now = (CoinCpuTime () - now)) > 10.)
      jnlst_ -> Printf (Ipopt::J_WARNING, J_PROBLEM,
			"Couenne: reading time %.3fs\n", now);
  }

  // create expression set for binary search
  auxSet_ = new std::set <exprAux *, compExpr>;

  if (base) {
    std::string s;
    base -> options() -> GetStringValue ("use_quadratic", s, "couenne."); 
    useQuadratic_ = (s == "yes");
  }

  if (base) {

    std::string s;

    base -> options() -> GetStringValue ("feasibility_bt",  s, "couenne."); doFBBT_ = (s == "yes");
    base -> options() -> GetStringValue ("optimality_bt",   s, "couenne."); doOBBT_ = (s == "yes");
    base -> options() -> GetStringValue ("aggressive_fbbt", s, "couenne."); doABT_  = (s == "yes");

    base -> options() -> GetIntegerValue ("log_num_obbt_per_level", logObbtLev_, "couenne.");
    base -> options() -> GetIntegerValue ("log_num_abt_per_level",  logAbtLev_,  "couenne.");

    base -> options() -> GetNumericValue ("feas_tolerance",  feas_tolerance_, "couenne.");
    base -> options() -> GetNumericValue ("opt_window",      opt_window_,     "couenne.");
  }
}


/// preprocess problem in order to extract linear relaxations etc.
void CouenneProblem::reformulate () {

  double now = CoinCpuTime ();

  if (domain_.current () == NULL) {

    // create room for problem's variables and bounds, if no domain exists
    CouNumber 
      *x  = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber)),
      *lb = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber)),
      *ub = (CouNumber *) malloc ((nVars()) * sizeof (CouNumber));

    for (int i = nVars(); i--;) {
      x  [i] =  0.;
      lb [i] = -COUENNE_INFINITY;
      ub [i] =  COUENNE_INFINITY;
    }

    domain_.push (nVars (), x, lb, ub);
  }

  // link initial variables to problem's domain
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)
    (*i) -> linkDomain (&domain_);

  if (jnlst_ -> ProduceOutput(Ipopt::J_SUMMARY, J_PROBLEM))
    print (std::cout);

  // save -- for statistic purposes -- number of original
  // constraints. Some of them will be deleted as definition of
  // auxiliary variables.
  nOrigCons_    = constraints_. size ();
  nOrigIntVars_ = nIntVars ();

  jnlst_->Printf (Ipopt::J_SUMMARY, J_PROBLEM,
		  "Problem size before reformulation: %d variables (%d integer), %d constraints.\n",
		  nOrigVars_, nOrigIntVars_, nOrigCons_);

  // reformulation
  if (!standardize ()) { // problem is infeasible if standardize returns false

    jnlst_->Printf(Ipopt::J_WARNING, J_PROBLEM,
		   "Couenne: problem infeasible after reformulation\n");
    // fake infeasible bounds for Couenne to bail out
    for (int i = nVars (); i--;)
      Ub (i) = - (Lb (i) = 1.);
  }

  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();

  // give a value to all auxiliary variables. Do it now to be able to
  // recognize complementarity constraints in fillDependence()
  initAuxs ();

  // fill dependence_ structure
  fillDependence (bonBase_);

  // quadratic handling
  fillQuadIndices ();

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_->Printf(Ipopt::J_WARNING, J_PROBLEM,
    "Couenne: reformulation time %.3fs\n", now);

  // give a value to all auxiliary variables
  initAuxs ();

  int nActualVars = nIntVars_ = 0;

  // check how many integer variables we have now (including aux)
  for (int i=0; i<nVars(); i++)
    if (variables_ [i] -> Multiplicity () > 0) {

      nActualVars++;
      if (variables_ [i] -> isDefinedInteger ())
	nIntVars_++;
    }

  jnlst_->Printf(Ipopt::J_SUMMARY, J_PROBLEM,
		  "Problem size after  reformulation: %d variables (%d integer), %d constraints.\n",
		  nActualVars, nIntVars_, nCons());

  // check if optimal solution is available (for debug purposes)
  readOptimum ();

  if (bonBase_) {

    CouNumber 
      art_cutoff =  COIN_DBL_MAX,
      art_lower  = -COIN_DBL_MAX;

    bonBase_ -> options() -> GetNumericValue ("art_cutoff", art_cutoff, "couenne.");
    bonBase_ -> options() -> GetNumericValue ("art_lower",  art_lower,  "couenne.");

    if (art_cutoff <  1.e50) setCutOff (art_cutoff);
    if (art_lower  > -1.e50) {
      int indobj = objectives_ [0] -> Body () -> Index ();
      if (indobj >= 0)
	domain_.lb (indobj) = art_lower;
    }
  }

  if (jnlst_->ProduceOutput(Ipopt::J_DETAILED, J_PROBLEM)) {
    // We should route that also through the journalist
    print (std::cout);
  }

  //writeAMPL ("extended-aw.mod", true);
  //writeAMPL ("original.mod", false);
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
  doOBBT_       (p. doOBBT_),
  doABT_        (p. doABT_),
  logObbtLev_   (p. logObbtLev_),
  logAbtLev_    (p. logAbtLev_),
  jnlst_        (p.jnlst_),
  opt_window_   (p.opt_window_),    // needed only in standardize (), unnecessary to update it
  useQuadratic_ (p.useQuadratic_),  // ditto
  feas_tolerance_ (p.feas_tolerance_),
  dependence_   (p.dependence_),
  objects_      (p.objects_),
  integerRank_  (NULL),
  numberInRank_ (p.numberInRank_),
  maxCpuTime_   (p.maxCpuTime_),
  bonBase_      (p.bonBase_) {

  for (int i=0; i < p.nVars (); i++)
    variables_ . push_back (NULL);

  for (int i=0; i < p.nVars (); i++) {
    int ind = p.numbering_ [i];
    variables_ [ind] = p.Var (ind) -> clone (&domain_);
  }

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

  // delete variables
  //for (std::vector <exprVar *>::iterator i = variables_ . begin ();
  //i != variables_ . end (); ++i)
  //delete (*i);

  for (int i=nVars (); i--;) { // delete in inverse order
    int ind = numbering_ [i];
    delete variables_ [ind];
  }

  // delete extra structures
  if (graph_)     delete    graph_;
  if (commuted_)  delete [] commuted_;
  if (numbering_) delete [] numbering_;

  if (created_pcutoff_) delete pcutoff_;

  if (integerRank_) delete [] integerRank_;
}


/// methods to add objective function. 

void CouenneProblem::addObjective (expression *newobj, const std::string &sense = "min") {
  objectives_ . push_back
    (new CouenneObjective ((sense == "min") ? newobj : new exprOpp (newobj)));//, MINIMIZE));
}


/// methods to add nonlinear constraints:

/// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs = NULL) {

  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

/// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (COUENNE_INFINITY)));
}

/// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (-COUENNE_INFINITY), rhs));
}

/// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb=NULL, expression *ub=NULL) {
  if (!lb) lb = new exprConst (0.);
  if (!ub) ub = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}


/// add variable to the problem -- check whether it is integer (isDiscrete)

expression *CouenneProblem::addVariable (bool isDiscrete, Domain *d) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size (), d)) :
    (new exprVar  (variables_ . size (), d));

  variables_ . push_back (var);

  if (isDiscrete) 
    nIntVars_++;

  return var;
}


/// add auxiliary variable and associate it with pointer to expression
/// given as argument
exprAux *CouenneProblem::addAuxiliary (expression *symbolic) {

  // check if image is already in the expression database auxSet_
  std::set <exprAux *, compExpr>::iterator i;

  int var_ind = variables_ . size ();
  domain_. current () -> resize (var_ind + 1);

  symbolic -> getBounds (domain_. lb (var_ind), 
			 domain_. ub (var_ind));

  // create new aux associated with that expression
  exprAux *w = new exprAux (symbolic,
			    var_ind,
			    1 + symbolic -> rank (), 
			    exprAux::Unset, 
			    &domain_);
  //symbolic -> isInteger () ? exprAux::Integer : exprAux::Continuous);

  //  w -> linkDomain (&domain_);

  // seek expression in the set
  if ((i = auxSet_ -> find (w)) == auxSet_ -> end ()) {

    // no such expression found in the set, create entry therein
    variables_ . push_back (w);
    auxSet_ -> insert (w); // insert into repetition checking structure
    graph_  -> insert (w); // insert into acyclic structure

  } else {  // otherwise, just return the entry's pointer

    delete w;
    w = *i;
    (*i) -> increaseMult ();
  }

  return w;
}


/// translates pair (indices, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexL, 
				    CouNumber *coeff,
				    std::vector <std::pair <exprVar *, CouNumber> > &lcoeff) {
  // TODO: sort

  for (int i=0; indexL [i] >= 0; i++)
    lcoeff.push_back (std::pair <exprVar *, CouNumber> (Var (indexL [i]), coeff [i]));
}


/// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexI,
				    int *indexJ,
				    CouNumber *coeff,
				    std::vector <quadElem> &qcoeff) {
  // TODO: sort

  for (int i=0; indexI [i] >= 0; i++)
    qcoeff.push_back (quadElem (Var (indexI [i]), Var (indexJ [i]), coeff [i]));
}


/// fill in the integerRank_ array
void CouenneProblem::fillIntegerRank () const {

  if (integerRank_)
    return;

  int nvars = nVars ();

  integerRank_ = new int [nvars];

  // 0: fractional
  // 1: integer
  // k: integer,    depending on at least one integer with associated value k-1, or
  // k: fractional, depending on at least one integer with associated value k

  for (int ii = 0; ii < nvars; ii++) {

    int index = numbering_ [ii];

    if (Var (index) -> Multiplicity () <= 0) {
      integerRank_ [index] = 0;
      continue;
    }

    bool isInt = Var (index) -> isDefinedInteger ();

    integerRank_ [index] = (isInt) ? 1 : 0;

    if (Var (index) -> Type () == AUX) {

      std::set <int> deplist;

      if (Var (index) -> Image () -> DepList (deplist, STOP_AT_AUX) != 0) // depends on something
	for (std::set <int>::iterator i = deplist.begin (); i != deplist.end (); ++i) {

	  int token = integerRank_ [*i];
	  if (isInt) token++;

	  if (token > integerRank_ [index]) // there's a free integer below us
	    integerRank_ [index] = token;
	}
    }
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "Free (original) integers\n");
  for (int i=0; i<nOrigVars_; i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, integerRank_ [i]);

  // fill in numberInRank_
  for (int i=0; i<nOrigVars_; i++)
    if ((variables_ [i] -> isDefinedInteger ()) &&
	(variables_ [i] -> Multiplicity () > 0)) {

    int rank = integerRank_ [i];

    if (numberInRank_.size () <= (unsigned int) rank)
      for (int j=numberInRank_.size (); j <= rank; j++)
	numberInRank_ .push_back (0);

    numberInRank_ [rank] ++;
  }

  jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "numInteger [neglect non-originals]\n");
  for (unsigned int i=0; i<numberInRank_.size(); i++)
    jnlst_->Printf (Ipopt::J_VECTOR, J_PROBLEM, "%d: %d\n", i, numberInRank_ [i]);
}
