/*
 *
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

#include "CoinHelperFunctions.hpp"

using namespace Ipopt;
using namespace Couenne;

static const double init_estimate = COUENNE_EPS;

/// Empty constructor
CouenneObject::CouenneObject ():

  OsiObject (),
  cutGen_          (NULL),
  problem_         (NULL),
  reference_       (NULL),
  strategy_        (MID_INTERVAL),
  jnlst_           (NULL),
  alpha_           (default_alpha),
  lp_clamp_        (default_clamp),
  feas_tolerance_  (feas_tolerance_default),
  doFBBT_          (true),
  doConvCuts_      (true),
  downEstimate_    (init_estimate),
  upEstimate_      (init_estimate),
  pseudoMultType_  (INFEASIBILITY) {}


/// Constructor with information for branching point selection strategy
CouenneObject::CouenneObject (CouenneCutGenerator *cutgen,
			      CouenneProblem *p,
			      exprVar *ref,
			      Bonmin::BabSetupBase *base,
			      JnlstPtr jnlst):
  OsiObject (),

  cutGen_         (cutgen),
  problem_        (p),
  reference_      (ref),
  strategy_       (MID_INTERVAL),
  jnlst_          (jnlst),
  alpha_          (default_alpha),
  lp_clamp_       (default_clamp),
  feas_tolerance_ (feas_tolerance_default),
  doFBBT_         (true),
  doConvCuts_     (true),
  downEstimate_   (init_estimate),
  upEstimate_     (init_estimate),
  pseudoMultType_ (INFEASIBILITY) {

  // read options
  setParameters (base);

  // debug output ////////////////////////////////////////////

  if (reference_ &&
      (reference_ -> Type () == AUX) &&
      jnlst_ -> ProduceOutput (J_SUMMARY, J_BRANCHING)) {

    printf ("created Expression Object: "); reference_ -> print ();
    if (reference_ -> Image ()) {
      printf (" := ");
      reference_ -> Image () -> print ();
    }

    printf (" with %s strategy [clamp=%g, alpha=%g]\n",
	    (strategy_ == LP_CLAMPED)   ? "lp-clamped" :
	    (strategy_ == LP_CENTRAL)   ? "lp-central" :
	    (strategy_ == BALANCED)     ? "balanced"   :
	    (strategy_ == MIN_AREA)     ? "min-area"   :
	    (strategy_ == MID_INTERVAL) ? "mid-point"  :
	    (strategy_ == NO_BRANCH)    ? "no-branching (null infeasibility)" :
	                                  "no strategy",
	    lp_clamp_, alpha_);
  }
}


/// Constructor with lesser information, used for infeasibility only
CouenneObject::CouenneObject (exprVar *ref,
			      Bonmin::BabSetupBase *base,
			      JnlstPtr jnlst):

  OsiObject (),
  cutGen_         (NULL),
  problem_        (NULL),
  reference_      (ref),
  strategy_       (MID_INTERVAL),
  jnlst_          (jnlst),
  alpha_          (default_alpha),
  lp_clamp_       (default_clamp),
  feas_tolerance_ (feas_tolerance_default),
  doFBBT_         (true),
  doConvCuts_     (true),
  downEstimate_   (init_estimate),
  upEstimate_     (init_estimate),
  pseudoMultType_ (INFEASIBILITY) {

  // read options
  setParameters (base);
}


/// Copy constructor
CouenneObject::CouenneObject (const CouenneObject &src):
  OsiObject       (src),
  cutGen_         (src.cutGen_),
  problem_        (src.problem_),
  reference_      (src.reference_),
  strategy_       (src.strategy_),
  jnlst_          (src.jnlst_),
  alpha_          (src.alpha_),
  lp_clamp_       (src.lp_clamp_),
  feas_tolerance_ (src.feas_tolerance_),
  doFBBT_         (src.doFBBT_),
  doConvCuts_     (src.doConvCuts_),
  downEstimate_   (src.downEstimate_),
  upEstimate_     (src.upEstimate_),
  pseudoMultType_ (src.pseudoMultType_) {}


/// integer infeasibility: min {value - floor(value), ceil(value) -
/// value}. Handle special cases so that Cbc and Couenne agree. Also,
/// check out-of-bounds values.

double CouenneObject::intInfeasibility (double value, double lb, double ub) const {

  if      (value < lb) value = lb;
  else if (value > ub) value = ub;

  return CoinMin (value - floor (value + COUENNE_EPS),
		          ceil  (value - COUENNE_EPS) - value);
}


/// apply the branching rule
OsiBranchingObject *CouenneObject::createBranch (OsiSolverInterface *si,
						 const OsiBranchingInformation *info,
						 int way) const {

  // a nonlinear constraint w = f(x) is violated. The infeasibility
  // should be given by something more elaborate than |w-f(x)|, that
  // is, it is the minimum, among the two branching nodes, of the
  // distance from the current optimum (w,x) and the optimum obtained
  // after convexifying the two subproblems. We call selectBranch for
  // the purpose, and save the output parameter into the branching
  // point that should be used later in createBranch.

  if (jnlst_ -> ProduceOutput (J_ITERSUMMARY, J_BRANCHING)) {
    printf ("CouObj::createBranch on ");
    reference_ -> print (); printf ("\n");
  }

  // copy current point into Couenne
  problem_ -> domain () -> push
    (problem_ -> nVars (),
     info -> solution_,
     info -> lower_,
     info -> upper_,
     false);

  CouNumber
    *brPts  = NULL,         // branching point(s)
    *brDist = NULL;         // distances from current LP point to each
			    // new convexification (usually two)
  expression *brVar = NULL; // branching variable
  int whichWay = 0;

  CouNumber improv;

  if (reference_ -> Image ())  // not used out of debug
    improv = reference_ -> Image () ->
      selectBranch (this, info,                      // input parameters
		    brVar, brPts, brDist, whichWay); // result: who, where, distance, and direction
  else {

    // reference_ -> Image () == NULL implies this is an integer object

    brVar = reference_;
    brPts  = (double *) realloc (brPts,      sizeof (double));
    brDist = (double *) realloc (brDist, 2 * sizeof (double));

    double point = info -> solution_ [reference_ -> Index ()];

    *brPts = point;
    improv = 0.;

    if (point > floor (point)) {improv =                  brDist [0] = point - floor (point);}
    if (point < ceil  (point)) {improv = CoinMin (improv, brDist [1] = ceil (point) - point);}

    point -= floor (point);
    whichWay = (point < 0.45) ? TWO_LEFT : (point > 0.55) ? TWO_RIGHT : TWO_RAND;
  }

  assert (brVar); // MUST have a branching variable

  /*  if        (pseudoMultType_ == INTERVAL) {
    int index = brVar -> Index ();
    downEstimate_ = *brPts - info -> lower_ [index];
    upEstimate_   =          info -> upper_ [index] - *brPts;
  } else if (pseudoMultType_ == REV_INTERVAL) {
    int index = brVar -> Index ();
    downEstimate_ = *brPts - info -> upper_ [index];          // notice upper&lower inverted
    upEstimate_   =          info -> lower_ [index] - *brPts;
  } else
  */

  if (pseudoMultType_ == PROJECTDIST) {
    downEstimate_ = brDist [0];
    upEstimate_   = brDist [1];
  } else setEstimates (info, NULL, brPts);

  /// Debug output
  if (jnlst_ -> ProduceOutput (J_MOREMATRIX, J_BRANCHING)) {
    printf ("brpts for "); reference_ -> print ();
    if (reference_ -> Image ()) {printf (" := "); reference_ -> Image () -> print ();}
    printf (" is on "); brVar -> print ();
    printf (" @ %.12g [%.12g,%.12g]\n", *brPts,
	    problem_ -> Lb (brVar -> Index ()),
	    problem_ -> Ub (brVar -> Index ()));

    if (brVar) {

      if (improv <= COUENNE_EPS) {
	printf ("### warning, infeas = %g for ", improv);
	reference_ -> print ();
	if (reference_ -> Image ()) {printf (":="); reference_ -> Image () -> print ();}
	printf ("\n");
      }

      int index = brVar -> Index ();
      if (info -> lower_ [index] >=
	  info -> upper_ [index] - COUENNE_EPS) {
	printf ("### warning, tiny bounding box [%g,%g] for x_%d\n",
		info -> lower_ [index],
		info -> upper_ [index], index);
      }
    }
  }

  // create branching object ///////////////////////////////////////
  OsiBranchingObject *brObj = new CouenneBranchingObject
    (si, this, jnlst_, cutGen_, problem_, brVar, way, *brPts, doFBBT_, doConvCuts_);

  problem_ -> domain () -> pop (); // Couenne discards current point

  if (brPts)  free (brPts);
  if (brDist) free (brDist);

  return brObj;
}


/// computes a not-too-bad point where to branch, in the "middle" of an interval
CouNumber CouenneObject::midInterval (CouNumber x, CouNumber l, CouNumber u,
				      const OsiBranchingInformation *info) const {

  CouNumber curAlpha = alpha_;

#if 1
  if (info) {

    // adaptive scheme: make LP point count more when gap approaches zero
    int objInd = problem_ -> Obj (0) -> Body () -> Index ();

    double
      lb = objInd >= 0 ? info -> lower_ [objInd] : problem_ -> Obj (0) -> Body () -> Value (),
      ub = problem_ -> getCutOff (),
      currentGap =
      (ub >  COUENNE_INFINITY    / 10 ||
       lb < -Couenne_large_bound / 10) ? 1.e3 :
      fabs (ub - lb) / (1.e-3 + CoinMin (fabs (ub), fabs (lb)));

#if 1
    if (currentGap < 1e-3) {

      currentGap *= 1e3;

      assert ((currentGap >= 0.) &&
	      (currentGap <= 1.));

      curAlpha = currentGap * alpha_ + (1 - currentGap);
    }
#else
    // make curAlpha closer to 1 by adding remaining (1-alpha_)
    // weighted inversely proportional to gap
    curAlpha += (1 - alpha_) / (1. + currentGap);
#endif

    // printf ("using %g rather than %g. cutoff %g, lb %g, gap %g\n",
    // 	    curAlpha, alpha_,
    // 	    problem_ -> getCutOff (),
    // 	    lb, currentGap);
  }
#endif

  if (u < l + COUENNE_EPS)
    return (0.5 * (l + u));

  if      (x<l) x = l;
  else if (x>u) x = u;

  if   (l < -large_bound)
    if (u >  COUENNE_EPS) return 0.;                                        // ]-inf,+inf[
    else                  return CoinMax ((l+u)/2, (AGGR_MUL * (-1. + u))); // ]-inf,u]
  else
    if (u >  large_bound)                                                   // [l,+inf[
      if (l < - COUENNE_EPS) return 0.;
      else                   return CoinMin ((l+u)/2, (AGGR_MUL * (1. + l)));
    else {                                                                  // [l,u]

      CouNumber point = curAlpha * x + (1. - curAlpha) * (l + u) / 2.;

      if      ((point-l) / (u-l) < closeToBounds) point = l + (u-l) * closeToBounds;
      else if ((u-point) / (u-l) < closeToBounds) point = u + (l-u) * closeToBounds;

      return point;
    }
}


/// pick branching point based on current strategy
CouNumber CouenneObject::getBrPoint (funtriplet *ft, CouNumber x0, CouNumber l, CouNumber u,
				     const OsiBranchingInformation *info) const {

  if ((l < -COUENNE_EPS) &&
      (u >  COUENNE_EPS) &&
      (-l/u >= THRES_ZERO_SYMM) &&
      (-u/l >= THRES_ZERO_SYMM))
    return 0.; // this happens when [l,u] is quite symmetrical w.r.t. 0

  CouNumber width = lp_clamp_ * (u-l);

  switch (strategy_) {

  case CouenneObject::MIN_AREA:     return maxHeight   (ft, l, u);
  case CouenneObject::BALANCED:     return minMaxDelta (ft, l, u);
  case CouenneObject::LP_CLAMPED:   return CoinMax (l + width, CoinMin (x0, u - width));
  case CouenneObject::LP_CENTRAL:   return ((x0 < l + width) || (x0 > u - width)) ? (l+u)/2 : x0;
  case CouenneObject::MID_INTERVAL: return midInterval (x0, l, u, info);

  default:
    printf ("Couenne: unknown branching point selection strategy\n");
    exit (-1);
  }
}


/// non-linear infeasibility -- not called by independent's CouenneVarObject
double CouenneObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  if (strategy_ == NO_BRANCH) {
    upEstimate_ = downEstimate_ = init_estimate;
    return 0;
  }

  problem_ -> domain () -> push
    (problem_ -> nVars (),
     info -> solution_,
     info -> lower_,
     info -> upper_,
     false);

  double retval = checkInfeasibility (info);

  problem_ -> domain () -> pop ();

  bool isInt = reference_ -> isInteger ();

  int refInd = reference_ -> Index ();
  double point = info -> solution_ [refInd];

  if (pseudoMultType_ == INFEASIBILITY) {

    if (isInt) {
      CouNumber intInfeas = intInfeasibility (point, info -> lower_ [refInd], info -> upper_ [refInd]);

      if (retval < intInfeas) {
	if (downEstimate_ <       point  - floor (point)) downEstimate_ =       point  - floor (point);
	if (upEstimate_   < ceil (point) -        point)  upEstimate_   = ceil (point) -        point;
	retval = intInfeas;
      }
    }
    else upEstimate_ = downEstimate_ = retval;
  }
  else setEstimates (info, &retval, NULL);

  return (isInt ?
    CoinMax (retval, intInfeasibility (point, info -> lower_ [refInd], info -> upper_ [refInd])) :
    retval);
}


/// non-linear infeasibility -- no need for the domain's push
/// instruction as this is called from
/// CouenneVarObject::infeasibility()
double CouenneObject::checkInfeasibility (const OsiBranchingInformation *info) const {

  int refInd = reference_ -> Index ();

  if (reference_ -> Type () == VAR)
    return (reference_ -> isInteger ()) ?
      intInfeasibility (info -> solution_ [refInd],
			info -> lower_    [refInd],
			info -> upper_    [refInd]) : 0.;

  double
    vval = info -> solution_ [reference_ -> Index ()],
    fval = (*(reference_ -> Image ())) (),
    denom  = CoinMax (1., reference_ -> Image () -> gradientNorm (info -> solution_));

  // check if fval is a number (happens with e.g. w13 = w12/w5 and w5=0, see test/harker.nl)
  if (CoinIsnan (fval)) {
    fval = vval + 1.;
    denom = 1.;
  }

  if (fabs (fval) > COUENNE_INFINITY)
    fval = COUENNE_INFINITY;

  double
    retval =
    ((reference_ -> sign () == expression::AUX_GEQ) && (vval >= fval)) ? 0. :
    ((reference_ -> sign () == expression::AUX_LEQ) && (vval <= fval)) ? 0. : fabs (vval - fval),

    ratio = (CoinMax (1., fabs (vval)) /
	     CoinMax (1., fabs (fval)));

  //jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "checkinf --> v=%e f=%e den=%e ret=%e ratio=%e\n", vval, fval, denom, retval, ratio);

  if ((ratio < 2)  &&
      (ratio > .5) &&
      ((retval /= denom) < CoinMin (COUENNE_EPS, feas_tolerance_)))
    retval = 0.;

  if (//(retval > 0.) &&
      (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING))) {

    printf ("  infeas %g: ", retval);
    reference_             -> print ();
    if (reference_ -> Image ()) {printf (" := "); reference_ -> Image () -> print ();}
    printf ("\n");
  }

  // Need to allow infeasibility for a variable without making the
  // whole problem infeasible with an infeasibility = 1e+50. Check
  // BonChooseVariable.cpp:382

  if (retval > 1.e40)
    retval = 1.e20;

  return (reference_ -> isInteger ()) ?
    CoinMax (retval, intInfeasibility (info -> solution_ [refInd],
				       info -> lower_    [refInd],
				       info -> upper_    [refInd])) :
    retval;
}


/// set parameters by reading command line or parameter file
void CouenneObject::setParameters (Bonmin::BabSetupBase *base) {

  if (!base) return;

  std::string s;

  base -> options () -> GetStringValue ("pseudocost_mult", s, "couenne.");

  if      (s == "interval_lp")     pseudoMultType_ = INTERVAL_LP;
  else if (s == "interval_lp_rev") pseudoMultType_ = INTERVAL_LP_REV;
  else if (s == "interval_br")     pseudoMultType_ = INTERVAL_BR;
  else if (s == "interval_br_rev") pseudoMultType_ = INTERVAL_BR_REV;
  else if (s == "infeasibility")   pseudoMultType_ = INFEASIBILITY;
  else if (s == "projectDist")     pseudoMultType_ = PROJECTDIST;

  base -> options() -> GetStringValue ("branch_fbbt",      s, "couenne."); doFBBT_     = (s=="yes");
  base -> options() -> GetStringValue ("branch_conv_cuts", s, "couenne."); doConvCuts_ = (s=="yes");

  base -> options() -> GetNumericValue ("feas_tolerance", feas_tolerance_, "couenne.");

  std::string brtype;
  base -> options () -> GetStringValue ("branch_pt_select", brtype, "couenne.");

  if      (brtype == "balanced")    strategy_ = BALANCED;
  else if (brtype == "lp-clamped")  strategy_ = LP_CLAMPED;
  else if (brtype == "lp-central")  strategy_ = LP_CENTRAL;
  else if (brtype == "min-area")    strategy_ = MIN_AREA;
  else if (brtype == "no-branch")   strategy_ = NO_BRANCH;
  else if (brtype == "mid-point") {
    strategy_ = MID_INTERVAL;
    base -> options () -> GetNumericValue ("branch_midpoint_alpha", alpha_, "couenne.");
  }

  if (strategy_ == LP_CLAMPED ||
      strategy_ == LP_CENTRAL)
    base -> options () -> GetNumericValue ("branch_lp_clamp", lp_clamp_, "couenne.");

  // accept options for branching rules specific to each operator
  if (reference_ && reference_ -> Type () == AUX) {

    std::string br_operator = "";

    switch (reference_ -> Image () -> code ()) {

    case COU_EXPRPOW: {

      // begin with default value in case specific exponent are not given
      base -> options () -> GetStringValue ("branch_pt_select_pow", brtype, "couenne.");

      CouNumber expon = reference_ -> Image () -> ArgList () [1] -> Value ();

      if      (fabs (expon - 2.) < COUENNE_EPS) br_operator = "sqr";
      else if (fabs (expon - 3.) < COUENNE_EPS) br_operator = "cube";
      else if (expon             < 0.)          br_operator = "negpow";
      else                                      br_operator = "pow";
    } break;

    case COU_EXPRMUL:
      br_operator = (reference_ -> Image () -> ArgList () [0] -> Index () !=
		     reference_ -> Image () -> ArgList () [1] -> Index ()) ?
	"prod" : "sqr";
      break;
    case COU_EXPRINV: br_operator = "negpow"; break;
    case COU_EXPRDIV: br_operator = "div"; break;
    case COU_EXPRLOG: br_operator = "log"; break;
    case COU_EXPREXP: br_operator = "exp"; break;
    case COU_EXPRSIN:
    case COU_EXPRCOS: br_operator = "trig"; break;
    default: break;
    }

    if (br_operator != "") {
      // read option
      char select [40], sel_clamp [40];
      double lp_clamp_fun = default_clamp;
      sprintf (select,    "branch_pt_select_%s", br_operator.c_str ());
      sprintf (sel_clamp, "branch_lp_clamp_%s",  br_operator.c_str ());
      base -> options () -> GetStringValue (select, brtype, "couenne.");
      base -> options () -> GetNumericValue (sel_clamp, lp_clamp_fun, "couenne.");

      if (lp_clamp_fun != default_clamp)
	lp_clamp_ = lp_clamp_fun;

      if      (brtype == "balanced")    strategy_ = BALANCED;
      else if (brtype == "lp-clamped")  strategy_ = LP_CLAMPED;
      else if (brtype == "lp-central")  strategy_ = LP_CENTRAL;
      else if (brtype == "min-area")    strategy_ = MIN_AREA;
      else if (brtype == "no-branch")   strategy_ = NO_BRANCH;
      else if (brtype == "mid-point") {
	strategy_ = MID_INTERVAL;
	double alpha_fun = default_alpha;
	base -> options () -> GetNumericValue ("branch_midpoint_alpha", alpha_fun, "couenne.");
	if (alpha_fun != default_alpha)
	  alpha_ = alpha_fun;
      }
    }
  }
}


/// set up/down estimates based on branching information
void CouenneObject::setEstimates (const OsiBranchingInformation *info,
				  CouNumber *infeasibility,
				  CouNumber *brpoint) const {

  int index = reference_ -> Index ();

  //bool isInteger = reference_ -> isInteger ();

  CouNumber
    *up   = &upEstimate_,
    *down = &downEstimate_,
     point = 0.,
     lower = info -> lower_ [index],
     upper = info -> upper_ [index];

  ////////////////////////////////////////////////////////////
  //
  // these rules invert the interval, so we just flip the pointers

  if ((pseudoMultType_ == INTERVAL_LP_REV) ||
      (pseudoMultType_ == INTERVAL_BR_REV)) {

    up   = &downEstimate_;
    down = &upEstimate_;
  }

  ///////////////////////////////////////////////////////////
  if (info &&
      ((pseudoMultType_ == INTERVAL_LP) ||
       (pseudoMultType_ == INTERVAL_LP_REV)))

    point = info -> solution_ [index];

  else if (brpoint &&
	   ((pseudoMultType_ == INTERVAL_BR) ||
	    (pseudoMultType_ == INTERVAL_BR_REV)))

    point = *brpoint;

  //printf ("point = %g [%s]\n", point, reference_ -> isInteger () ? "int" : "frac");

  // now move it away from the bounds, unless this is an integer variable

  //if (!isInteger)
  point = midInterval (point, lower, upper, info);

  //printf ("point = %g\n", point);

  if ((lower > -COUENNE_INFINITY) &&
      (upper <  COUENNE_INFINITY)) {

    CouNumber delta = closeToBounds * (upper - lower);

    if      (point < lower + delta)
      point        = lower + delta;
    else if (point > upper - delta)
      point        = upper - delta;
  }

  ///////////////////////////////////////////////////////////
  switch (pseudoMultType_) {

  case INFEASIBILITY:

    if (infeasibility)
      upEstimate_ = downEstimate_ = *infeasibility;

    break;

  case INTERVAL_LP:
  case INTERVAL_LP_REV:
  case INTERVAL_BR:
  case INTERVAL_BR_REV:
    assert (info);
    *up   = CoinMin (max_pseudocost, COUENNE_EPS + fabs (upper - point));
    *down = CoinMin (max_pseudocost, COUENNE_EPS + fabs (        point - lower));
    break;

  case PROJECTDIST: // taken care of in selectBranch procedure
    break;

  default:
    printf ("Couenne: invalid estimate setting procedure\n");
    exit (-1);
  }

  /*if (downEstimate_ <= 0.0 || upEstimate_ <= 0.0)
     printf ("%g [%g,%g] ---> [%g,%g]\n",
 	    point,
 	    lower,
 	    upper,
 	    downEstimate_, upEstimate_);*/

  /*if (reference_ -> isInteger ()) {

    CouNumber
      fracDn =       point  - floor (point),
      fracUp = ceil (point) -        point;

    if (downEstimate_ < fracDn) downEstimate_ = fracDn;
    if (upEstimate_   < fracUp)   upEstimate_ = fracUp;
    }*/

  assert (downEstimate_ > 0. &&
          upEstimate_   > 0.);
}
