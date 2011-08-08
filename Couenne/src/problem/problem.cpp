/* $Id$
 *
 * Name:    problem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <vector>

#include "BonRegisteredOptions.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "BonBabSetupBase.hpp"

#include "CouenneTypes.hpp"

#include "CouenneExpression.hpp"
#include "CouenneExprConst.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneProblem.hpp"
#include "CouenneGlobalCutOff.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneLQelems.hpp"

using namespace Couenne;

#define MAX_FBBT_ITER 3

/// save CouenneBase
void CouenneProblem::setBase (Bonmin::BabSetupBase *base) {
  bonBase_ = base;
  jnlst_   = base -> journalist ();
}

// tricky... smaller values cut the optimum in OS unitTest
const CouNumber SafeCutoff = 1e-4; 

// absolute difference
const CouNumber SafeDelta = 1e-2; 

/// initialize auxiliary variables from original variables in the
/// nonlinear problem

void CouenneProblem::initAuxs () const {

  domain_.current () -> resize (nVars ());

  // initially, auxiliary variables are unbounded, their bounds only
  // depending on their function

  int nvars = nVars ();

  for (int i=0; i < nvars; i++) {

    int indvar = variables_ [i] -> Index ();

    if (((variables_ [i] -> Type () == AUX) &&     // this is an auxiliary
	 (indvar >= nOrigVars_))            ||     // and not an original, originally
	(variables_ [i] -> Multiplicity () == 0))  // or a useless one
      //int index = variables_ [i] -> Index ();
      Lb (indvar) = - (Ub (indvar) = COIN_DBL_MAX);
  }

  // first initialize with values from constraints

  //Jnlst()->Printf(Ipopt::J_VECTOR, J_PROBLEM, "Initial bounds for aux (initAuxs):\n");

  for (std::vector <CouenneConstraint *>::const_iterator con = constraints_.begin ();
       con != constraints_.end (); ++con) {

    CouNumber
      lb = (*((*con) -> Lb ())) (),
      ub = (*((*con) -> Ub ())) ();

    int index = (*con) -> Body () -> Index ();

    assert (index >= 0);

    if ((Lb (index) = CoinMax (Lb (index), lb)) <= -COUENNE_INFINITY) Lb (index) = -COIN_DBL_MAX;
    if ((Ub (index) = CoinMin (Ub (index), ub)) >=  COUENNE_INFINITY) Ub (index) =  COIN_DBL_MAX;
  }

  // only one loop is sufficient here, since auxiliary variable are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  Jnlst () -> Printf (Ipopt::J_MOREMATRIX, J_PROBLEM, "InitAux -- assigning bounds\n");

  for (int j=0, i=nVars (); i--; j++) {

    int ord = numbering_ [j];

    // ignore these variables!
    if (variables_ [ord] -> Multiplicity () == 0) {
      Lb (ord) = - (Ub (ord) = COIN_DBL_MAX);
      X (ord) = 0.;
      continue;
    }

    exprVar *var = variables_ [ord];

    // and handle only those with nonzero multiplicity
    if (var -> Type () == AUX) {

      Jnlst () -> Printf (Ipopt::J_MOREMATRIX, J_PROBLEM, 
			  "w_%04d [%10g,%10g] ", ord, Lb (ord), Ub (ord));

      CouNumber l, u;

      var -> Image () -> getBounds (l, u);

      /*printf ("printing bounds: [%g %g]\n", Lb (ord), Ub (ord));
      var -> Lb () -> print (); printf ("\n");
      var -> Ub () -> print (); printf ("\n");*/

      Jnlst () -> Printf (Ipopt::J_MOREMATRIX, J_PROBLEM, 
			  " ( --> w_%04d [%10g,%10g] ) vs [%10g %10g]", 
			  ord, l, u, Lb (ord), Ub (ord));

      // set bounds 
      if ((var -> sign () != expression::AUX_LEQ) && ((Lb (ord) = CoinMax (Lb (ord), l)) <= -COUENNE_INFINITY)) Lb (ord) = -COIN_DBL_MAX;
      if ((var -> sign () != expression::AUX_GEQ) && ((Ub (ord) = CoinMin (Ub (ord), u)) >=  COUENNE_INFINITY)) Ub (ord) =  COIN_DBL_MAX;

      //if ((lb_ [ord] = (*(aux -> Lb ())) ()) <= -COUENNE_INFINITY) lb_ [ord] = -DBL_MAX;
      //if ((ub_ [ord] = (*(aux -> Ub ())) ()) >=  COUENNE_INFINITY) ub_ [ord] =  DBL_MAX;

      bool integer = var -> isDefinedInteger ();

      if (integer) {
	if (var -> sign () != expression::AUX_GEQ) Lb (ord) = ceil  (Lb (ord) - COUENNE_EPS);
	if (var -> sign () != expression::AUX_LEQ) Ub (ord) = floor (Ub (ord) + COUENNE_EPS);
      }

      X (ord) = CoinMax (Lb (ord), CoinMin (Ub (ord), (*(var -> Image ())) ()));

      Jnlst () -> Printf (Ipopt::J_MOREMATRIX, J_PROBLEM, 
			  " --> [%10g,%10g] (%g)\n", Lb (ord), Ub (ord), X (ord));
    }
  }

  restoreUnusedOriginals (); // no argument as the local x vector will be used
}


/// get auxiliary variables from original variables in the nonlinear
/// problem
void CouenneProblem::getAuxs (CouNumber * x) const {

  // set point at x, don't copy
  domain_.push (nVars (), x, domain_.lb (), domain_.ub (), false);

  // if there are common expressions, set them (they usually don't get
  // set by the AMPL interface as they are not really variables)
  if (ndefined_ > 0)
    for (int i = 0; i < nVars (); ++i) {
      int ii = numbering_ [i];
      if (ii >= nOrigVars_ - ndefined_ && 
	  ii <  nOrigVars_)
	X (ii) = (*(commonexprs_ [ii - nOrigVars_ + ndefined_])) ();
    }

  // set auxiliary w to f(x). This procedure is exact even though the
  // auxiliary variables have an incomplete image, i.e. they have been
  // decomposed previously, since they are updated with increasing
  // index.

  for (int j=0, i=nVars (); i--; j++) {

    int index = numbering_ [j];
    exprVar *var = variables_ [index];

    if (var -> Multiplicity () > 0) {

      CouNumber l, u;

      if (var -> Type () == AUX)
	var -> Image () -> getBounds (l,u);
      else {
	l = Lb (index);
	u = Ub (index);
      }

      if (var -> Type () == AUX) {

	CouNumber &x = X (index);

	bool isInt = var -> isDefinedInteger ();

	/*printf ("checking "); var -> print ();
	printf (" = %g = %g. Sign: %d, int: %d, [%g,%g]", 
		X (index), (*(var -> Image ())) (),
		var -> sign (), isInt, l, u);*/

	if ((var -> sign () == expression::AUX_EQ) &&
	    ((index >= nOrigVars_) ||
	     (index  < nOrigVars_ - ndefined_)))
	  x = (*(var -> Image ())) ();  // addresses of x[] and X() are equal
    
	x = 
	  CoinMax ((var -> sign () != expression::AUX_LEQ) ? (isInt ? ceil  (l - COUENNE_EPS) : l) : -COIN_DBL_MAX, 
	  CoinMin ((var -> sign () != expression::AUX_GEQ) ? (isInt ? floor (u + COUENNE_EPS) : u) :  COIN_DBL_MAX, x));

	// heuristic feasibility heuristic: if semiaux, round value to nearest integer if variable is integer

	if (isInt) {
	  if (var -> sign () == expression::AUX_GEQ) x = ceil  (x - COUENNE_EPS);
	  if (var -> sign () == expression::AUX_LEQ) x = floor (x + COUENNE_EPS);
	}

	//printf (" -> %g\n", X (index));
      }
    } else X (index) = 0.;
  }

  restoreUnusedOriginals ();

  domain_.pop ();
}


/// fill obj vector with coefficient of the (linearized) obj function
/// (depends on sense of optimization -- invert if sense()==MAXIMIZE)

void CouenneProblem::fillObjCoeff (double *&obj) {

  // linearized objective can be an exprAux, an exprSub, an exprGroup,
  // or an exprSum. In the last two cases, the components are
  // variables or constants

  expression *body = objectives_ [0] -> Body ();
  //int sense = objectives_ [0] -> Sense ();

  switch (body -> code ()) {

  case COU_EXPRVAR:   //
    obj [body -> Index ()] = 1; //(sense == MINIMIZE) ? 1 : -1;
    break;

  case COU_EXPRSUB: { // 

    expression **arglist = body -> ArgList ();

    obj [arglist [0] -> Index ()] =  1; //(sense == MINIMIZE) ?  1 : -1;
    obj [arglist [1] -> Index ()] = -1; //(sense == MINIMIZE) ? -1 :  1;

  } break;

  case COU_EXPRGROUP: { // 

    exprGroup *eg    = dynamic_cast <exprGroup *> (body -> isaCopy () ? 
						   body -> Copy () :
						   body);

    const exprGroup::lincoeff &lcoe = eg -> lcoeff ();

    //    if (sense == MINIMIZE) while (*index >= 0) obj [*index++] =  *coeff++;
    //    else                   while (*index >= 0) obj [*index++] = -*coeff++;      

    for (int n = (int) lcoe.size (), i=0; n--; i++)
      //exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el)
      obj [lcoe [i]. first -> Index ()] = lcoe [i]. second;
    //(sense == MINIMIZE) ? 
    //(lcoe [i]. second) : 
    //-(lcoe [i]. second);

  } // no break, as exprGroup is derived from exprSum

  case COU_EXPRSUM: { // 

    expression **arglist = body -> ArgList ();

    for (int i = body -> nArgs (); i--;)
      switch ((arglist [i]) -> code ()) {

      case COU_EXPRCONST: 
	break;

      case COU_EXPRVAR: 
	obj [arglist [i] -> Index ()] = 1; //(sense == MINIMIZE) ? 1 : -1;
	break;

      case COU_EXPRMUL: {

	expression **mulArgList = arglist [i] -> ArgList ();
	int index = mulArgList [0] -> Index ();

	if (index >= 0) obj [index]                      = mulArgList [1] -> Value ();
	else            obj [mulArgList [1] -> Index ()] = mulArgList [0] -> Value ();
      }	break;

      default: 
	Jnlst()->Printf(Ipopt::J_ERROR, J_PROBLEM,
			"Couenne: invalid element of sum\nAborting\n");
	exit (-1);
      }
  } break;

  case COU_EXPRCONST: break; // a constant objective

  default:
    Jnlst()->Printf(Ipopt::J_WARNING, J_PROBLEM,
		    "Couenne: warning, objective function not recognized\n");
    break;
  }
}


/// set cutoff from NLP solution
void CouenneProblem::setCutOff (CouNumber cutoff, const double *s) const {

  int indobj = objectives_ [0] -> Body () -> Index ();

  // AW: Should we use the value of the objective variable computed by 
  //     Couenne here?
  if ((indobj >= 0) && (cutoff < pcutoff_ -> getCutOff () - COUENNE_EPS)) {

    //if (fabs (cutoff - pcutoff_ -> getCutOff ()) > (1 + fabs (cutoff)) * 2 * SafeCutoff) // avoid too many printouts
    Jnlst () -> Printf (Ipopt::J_WARNING, J_PROBLEM, "Couenne: new MINLP solution, value %.10e\n", cutoff);
			//pcutoff_ -> getCutOff ());

    if (Var (indobj) -> isInteger ())
      pcutoff_    -> setCutOff (this, floor (cutoff + COUENNE_EPS), s);
    else pcutoff_ -> setCutOff (this, cutoff, s);
  }
}

/// Reset cutoff to a given value
void CouenneProblem::resetCutOff (CouNumber value) const {

  int indobj = objectives_ [0] -> Body () -> Index ();

  if ((indobj >= 0)) {

    if (Var (indobj) -> isInteger ())
      pcutoff_    -> setCutOff (this, floor (value + COUENNE_EPS), NULL);
    else pcutoff_ -> setCutOff (this, value, NULL);
  }
}


/// Tell problem that auxiliary related to obj has a cutoff, to be
/// used in bound tightening
void CouenneProblem::installCutOff () const {

  // all problem are assumed to be minimization
  double cutoff = pcutoff_ -> getCutOff();

  if (cutoff > COUENNE_INFINITY) 
    return;

  int indobj = objectives_ [0] -> Body () -> Index ();

  assert (indobj >= 0);

  //Jnlst () -> Printf (Ipopt::J_WARNING, J_PROBLEM,
  //"installing cutoff %.10e vs current ub %.10e\n",
  //cutoff, Ub (indobj));

  cutoff = (Var (indobj) -> isInteger ()) ?
    floor (cutoff + COUENNE_EPS) :
    (cutoff + CoinMin (SafeDelta, SafeCutoff * (1. + fabs (cutoff))));  // tolerance needed to retain feasibility

  if (cutoff < Ub (indobj))
    Ub (indobj) = cutoff;
}


// clear all spurious variables pointers not referring to the variables_ vector
void CouenneProblem::realign () {

  // link variables to problem's domain
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i) {

    (*i) -> linkDomain (&domain_);
    (*i) -> realign (this);
    if ((*i) -> Type () == AUX)
      (*i) -> Image () -> realign (this);
  }

  // link variables to problem's domain
  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) 
    (*i) -> Body () -> realign (this);


  // link variables to problem's domain
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i)
    (*i) -> Body () -> realign (this);
}


/// Add list of options to be read from file
void CouenneProblem::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddNumberOption
    ("art_cutoff",
     "Artificial cutoff",
     COIN_DBL_MAX,
     "Default value is infinity.");

  roptions -> AddNumberOption
    ("opt_window",
     "Window around known optimum",
     COIN_DBL_MAX,
     "Default value is infinity.");

  roptions -> AddStringOption2
    ("use_semiaux",
     "Use semiauxiliaries, i.e. auxiliaries defined as w >= f(x) rather than w := f(x))",
     "yes",
     "no",   "Only use auxiliaries assigned with \"=\" ",
     "yes",  "Use auxiliaries defined by w <= f(x), w >= f(x), and w = f(x)"
    );

  roptions -> AddStringOption2
    ("use_auxcons",
     "Use constraints-defined auxiliaries, i.e. auxiliaries w = f(x) defined by original constraints f(x) - w = 0",
     "yes",
     "no",   "",
     "yes",  ""
    );

  roptions -> AddStringOption2 
    ("redcost_bt",
     "Reduced cost bound tightening",
     "yes",
     "no","",
     "yes","",
     "This bound reduction technique uses the reduced costs of the LP in order to infer better variable bounds.");

  roptions -> AddStringOption2 
    ("use_quadratic",
     "Use quadratic expressions and related exprQuad class",
     "no",
     "no","Use an auxiliary for each bilinear term",
     "yes","Create only one auxiliary for a quadratic expression",
     "If enabled, then quadratic forms are not reformulated and therefore decomposed as a sum of auxiliary variables, each associated with a bilinear term, but rather taken as a whole expression. "
     "Envelopes for these expressions are generated through alpha-convexification."
    );

  roptions -> AddStringOption2 
    ("optimality_bt",
     "Optimality-based (expensive) bound tightening (OBBT)",
     "yes",
     "no","",
     "yes","",
     "This is another bound reduction technique aiming at reducing the solution set by looking at the initial LP relaxation. "
     "This technique is computationally expensive, and should be used only when necessary."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_obbt_per_level",
     "Specify the frequency (in terms of nodes) for optimality-based bound tightening.",
     -1,1,
     "\
If -1, apply at every node (expensive!). \
If 0, apply at root node only. \
If k>=0, apply with probability 2^(k - level), level being the current depth of the B&B tree.");

  roptions -> AddLowerBoundedIntegerOption
    ("max_fbbt_iter",
     "Number of FBBT iterations before stopping even with tightened bounds.",
     -1, MAX_FBBT_ITER,
     "Set to -1 to impose no upper limit");

  roptions -> AddStringOption2 
    ("aggressive_fbbt",
     "Aggressive feasibility-based bound tightening (to use with NLP points)",
     "yes",
     "no","",
     "yes","",
     "Aggressive FBBT is a version of probing that also allows to reduce the solution set, although it is not as quick as FBBT. "
     "It can be applied up to a certain depth of the B&B tree -- see ``log_num_abt_per_level''. "
     "In general, this option is useful but can be switched off if a problem is too large and seems not to benefit from it."
    );

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_abt_per_level",
     "Specify the frequency (in terms of nodes) for aggressive bound tightening.",
     -1,2,
     "\
If -1, apply at every node (expensive!). \
If 0, apply at root node only. \
If k>=0, apply with probability 2^(k - level), level being the current depth of the B&B tree.");

  roptions -> AddNumberOption
    ("art_lower",
     "Artificial lower bound",
     -COIN_DBL_MAX,
     "Default value is -COIN_DBL_MAX.");

  roptions -> AddStringOption3
    ("branching_object",
     "type of branching object for variable selection",
     "var_obj",
     "vt_obj",   "use Violation Transfer from Tawarmalani and Sahinidis",
     "var_obj",  "use one object for each variable",
     "expr_obj", "use one object for each nonlinear expression");

  roptions -> AddStringOption2 
    ("delete_redundant",
     "Eliminate redundant variables, which appear in the problem as x_k = x_h",
     "yes",
     "no","Keep redundant variables, making the problem a bit larger",
     "yes","Eliminate redundant variables (the problem will be equivalent, only smaller)");

  roptions -> AddStringOption4
    ("quadrilinear_decomp",
     "type of decomposition for quadrilinear terms (see work by Cafieri, Lee, Liberti)",
     "rAI",
     "rAI",     "Recursive decomposition in bilinear terms (as in Ryoo and Sahinidis): x5 = ((x1 x2) x3) x4)",
     "tri+bi",  "Trilinear and bilinear term: x5 = (x1 (x2 x3 x4))",
     "bi+tri",  "Bilinear, THEN trilinear term: x5 = ((x1 x2) x3 x4))",
     "hier-bi", "Hierarchical decomposition: x5 = ((x1 x2) (x3 x4))");
}
