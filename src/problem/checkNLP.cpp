/*
 *
 * Name:    checkNLP.cpp
 * Author:  Pietro Belotti
 *          Francois Margot
 * Purpose: check NLP feasibility of incumbent integer solution
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

// check if solution is MINLP feasible
bool CouenneProblem::checkNLP (const double *solution, double &obj, bool recompute) const {

  if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {

    printf ("Checking solution: %.12e (", obj);

    for (int i=0; i<nOrigVars_ - ndefined_; i++)
      printf ("%g ", solution [i]);
    printf (")\n");
  }

  // pre-check on original variables --- this is done after every LP,
  // and should be efficient
  for (int i=0; i < nOrigVars_ - ndefined_; i++) {

    if (variables_ [i] -> Multiplicity () <= 0) 
      continue;

    CouNumber val = solution [i];

    // check (original and auxiliary) variables' integrality

    exprVar *v = variables_ [i];

    if ((v -> Type ()      == VAR) &&
	(v -> isDefinedInteger ()) &&
	(v -> Multiplicity () > 0) &&
	(fabs (val - COUENNE_round (val)) > feas_tolerance_)) {

      Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM,
		      "checkNLP: integrality %d violated: %.6f [%g,%g]\n", 
		      i, val, domain_.lb (i), domain_.ub (i));

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "Done: (0)\n");

      return false;
    }
  }

  const int infeasible = 1;
  const int wrong_obj  = 2;

  // copy solution, evaluate the corresponding aux, and then replace
  // the original variables again for checking
  CouNumber *sol = new CouNumber [nVars ()];
  CoinZeroN (sol     + nOrigVars_ - ndefined_, nVars() - (nOrigVars_ - ndefined_));
  CoinCopyN (solution, nOrigVars_ - ndefined_, sol);
  getAuxs (sol);
  CoinCopyN (solution, nOrigVars_ - ndefined_, sol);

  // install NL solution candidate in evaluation structure
  domain_.push (nVars (), sol, domain_.lb (), domain_.ub (), false);

  if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {
    printf ("  checkNLP: %d vars -------------------\n    ", domain_.current () -> Dimension ());
    for (int i=0; i<domain_.current () -> Dimension (); i++) {
      if (i && !(i%5)) printf ("\n    ");
      printf ("%4d %16g [%16e %16e]  ", i, domain_.x (i), domain_.lb (i), domain_.ub (i));
    }
  }

  expression *objBody = Obj (0) -> Body ();

  // BUG: if Ipopt solution violates bounds of original variables and
  // objective depends on originals, we may have a "computed object"
  // out of bounds

  //CouNumber realobj = (*(objBody -> Image () ? objBody -> Image () : objBody)) ();
  CouNumber realobj = obj;

  if (objBody) 
    realobj = 
      (objBody -> Index () >= 0) ?
      sol [objBody -> Index ()] : 
      (*(objBody -> Image () ? objBody -> Image () : objBody)) ();

  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "  Objective: %.12e %.12e %.12e\n", 
		      realobj, objBody -> Index () >= 0 ? sol [objBody -> Index ()] : objBody -> Value (), 
		      (*(objBody -> Image () ? objBody -> Image () : objBody)) ());

  bool retval = true;

  try {

    // check if objective corresponds
    
    if (fabs (realobj - obj) / (1. + fabs (realobj)) > feas_tolerance_) {

      Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM,
		      "  checkNLP, false objective: computed %g != %g xQ (diff. %g)\n", 
		      realobj, obj, realobj - obj);

      if (!recompute)
	throw wrong_obj;
    }

    if (recompute)
      obj = realobj;

    Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "  recomputed: %.12e\n", obj);

    for (int i=0; i < nOrigVars_ - ndefined_; i++) {

      if (variables_ [i] -> Multiplicity () <= 0) 
	continue;

      CouNumber val = domain_.x (i);

      // check bounds

      if ((val > domain_.ub (i) + feas_tolerance_) ||
	  (val < domain_.lb (i) - feas_tolerance_)) {

	Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM,
			"  checkNLP: variable %d out of bounds: %.6f [%g,%g] (diff %g)\n", 
			i, val, domain_.lb (i), domain_.ub (i),
			CoinMax (fabs (val - domain_.lb (i)), 
				 fabs (val - domain_.ub (i))));
	throw infeasible;
      }

      // check (original and auxiliary) variables' integrality

      if (variables_ [i] -> isDefinedInteger () &&
	  (fabs (val - COUENNE_round (val)) > feas_tolerance_)) {

	Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM,
			"  checkNLP: integrality %d violated: %.6f [%g,%g]\n", 
			i, val, domain_.lb (i), domain_.ub (i));

	throw infeasible;
      }
    }

    // check ALL auxs

    for (int i=0; i < nVars (); i++) {

      exprVar *v = variables_ [i];

      if ((v -> Type         () != AUX) ||
	  (v -> Multiplicity () <= 0))
	continue;

      if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {
	printf ("    "); v -> print (); 
	CouNumber
	  val = (*(v)) (), 
	  img = (*(v -> Image ())) (), 
	  diff = fabs (val - img);
	printf (": val = %15g; img = %-15g ", val, img);
	if (diff > 1e-9)
	  printf ("[diff %12e] ", diff);
	//for (int j=0; j<nVars (); j++) printf ("%.12e ", (*(variables_ [j])) ());
	v -> Image () -> print (); 
	printf ("\n");
      }
      
      // check if auxiliary has zero infeasibility

      // same as in CouenneObject::checkInfeasibility -- main difference is use of gradientNorm()

      double 
	vval = (*v) (),
	fval = (*(v -> Image ())) (),
	denom  = CoinMax (1., v -> Image () -> gradientNorm (X ()));

      // check if fval is a number (happens with e.g. w13 = w12/w5 and w5=0, see test/harker.nl)
      if (CoinIsnan (fval)) {
	fval = vval + 1.;
	denom = 1.;
      }

      if (fabs (fval) > COUENNE_INFINITY)
	fval = COUENNE_INFINITY;

      double
	delta = 
	((v -> sign () == expression::AUX_GEQ) && (vval >= fval)) ? 0. : 
	((v -> sign () == expression::AUX_LEQ) && (vval <= fval)) ? 0. : fabs (vval - fval),

	ratio = (CoinMax (1., fabs (vval)) / 
		 CoinMax (1., fabs (fval)));

      // printf ("checkinf --> v=%e f=%e den=%e ret=%d ratio=%e delta=%e, delta/denom=%e, thres=%e [", 
      // 	      vval, fval, denom, retval, ratio, delta, delta/denom, CoinMin (COUENNE_EPS, feas_tolerance_));
      // v -> print ();
      // printf (" %c= ", v -> sign () == expression::AUX_LEQ ? '<' : 
      // 	               v -> sign () == expression::AUX_GEQ ? '>' : ':');
      // v -> Image () -> print ();
      // printf ("]\n");

      if ((delta > 0.) &&
	  (((ratio > 2.)  ||  // check delta > 0 to take into account semi-auxs
	    (ratio <  .5)) ||
	   ((delta /= denom) > CoinMin (COUENNE_EPS, feas_tolerance_)))) {

	Jnlst () -> Printf (Ipopt::J_MOREVECTOR, J_PROBLEM,
			    "  checkNLP: auxiliary %d violates tolerance %g by %g/%g = %g\n", 
			    i, CoinMin (COUENNE_EPS, feas_tolerance_), delta*denom, denom, delta);

	throw infeasible;
      }
    }

    // check constraints

    for (int i=0; i < nCons (); i++) {

      CouenneConstraint *c = Con (i);

      CouNumber
	body = (*(c -> Body ())) (),
	lhs  = (*(c -> Lb   ())) (),
	rhs  = (*(c -> Ub   ())) ();

      if (((rhs <  COUENNE_INFINITY) && (body > rhs + feas_tolerance_ * (1. + CoinMax (fabs (body), fabs (rhs))))) || 
	  ((lhs > -COUENNE_INFINITY) && (body < lhs - feas_tolerance_ * (1. + CoinMax (fabs (body), fabs (lhs)))))) {

	if (Jnlst () -> ProduceOutput (Ipopt::J_MOREVECTOR, J_PROBLEM)) {

	  printf ("  checkNLP: constraint %d violated (lhs=%+e body=%+e rhs=%+e, violation %g): ",
		  i, lhs, body, rhs, CoinMax (lhs-body, body-rhs));

	  c -> print ();
	}

	throw infeasible;
      }
    }
  }

  catch (int exception) {

    switch (exception) {

    case wrong_obj:
      retval = false;
      break;

    case infeasible:
    default:
      retval = false;
      break;
    }
  }

  delete [] sol;
  domain_.pop ();

  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "Done: %d\n", retval);

  return retval;
}

/************************************************************************/
// Recompute objective value for sol
double CouenneProblem::checkObj(const CouNumber *sol, const double &precision) 
  const {

  expression *objBody = Obj(0)->Body();

  // BUG: if Ipopt couSol violates bounds of original variables and
  // objective depends on originals, we may have a "computed object"
  // out of bounds

  //CouNumber realobj = (*(objBody -> Image () ? objBody -> Image () : objBody)) ();
  CouNumber realObj = 0;

  if (objBody) {
    realObj = 
      (objBody ->Index() >= 0) ?
      sol[objBody->Index()] : 
      (*(objBody->Image() ? objBody->Image() : objBody)) ();
    
    Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, 
			"%.12e %.12e %.12e ------------------------------\n", 
			realObj, objBody -> Index () >= 0 ? sol[objBody -> Index ()] : 0., 
			(*(objBody -> Image () ? objBody -> Image () : objBody)) ());
  } 
  else {
    printf("### ERROR: CouenneProblem::checkObj(): no objective body\n");
    exit(1);
  }
  return realObj;
}

/************************************************************************/
// check integrality of original vars in sol; return true if all
// original integer vars are within precision of an integer value
bool CouenneProblem::checkInt(const CouNumber *sol,
			      const int from, const int upto, 
			      const std::vector<int> listInt,
			      const bool origVarOnly,  
			      const bool stopAtFirstViol,  
			      const double precision, double &maxViol) const {

  bool isFeas = true;

  for(unsigned int i=0; i<listInt.size(); i++) {

    int ind = listInt[i];

    if((ind < from) || (variables_ [ind] -> Multiplicity () <= 0))
      continue;

    if(ind >= upto)
      break;

    CouNumber val = sol[ind];
    exprVar *v = variables_ [ind];

    if ((!origVarOnly) || (v -> Type () == VAR)) {

      double viol = fabs (val - COUENNE_round (val));
      if (viol > maxViol) maxViol = viol;
      if (viol > precision) {

	Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM,
			"checkInt(): integrality %d violated: %.6f [%g,%g]: integer distance %e > %e (by %e)\n", 
			i, val, domain_.lb (i), domain_.ub (i), 
			fabs (val - COUENNE_round (val)),  feas_tolerance_, 
			fabs (val - COUENNE_round (val)) - feas_tolerance_);

	Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkInt(): integrality %d violated: %.6f [%g,%g]\n", ind, val, domain_.lb (ind), domain_.ub (ind));

	isFeas = false;

	if (stopAtFirstViol)
	  break;
      }
    }
  }
  return(isFeas);
}

/************************************************************************/
// Check bounds; returns true iff feasible for given precision
bool CouenneProblem::checkBounds(const CouNumber *sol,
				 const bool stopAtFirstViol,  
				 const double precision, double &maxViol) const {

  bool isFeas = true;
  for(int i=0; i<nOrigVars_ - ndefined_; i++) {
    
    if (variables_[i]-> Multiplicity () <= 0) 
      continue;
    
    CouNumber val = domain_.x (i);
    double viol = 0;
    double violUb = val - domain_.ub (i);
    double violLb = domain_.lb (i) - val;

    if (viol < violUb) viol = violUb;
    if (viol < violLb) viol = violLb; 
    
    maxViol = (maxViol > viol ? maxViol : viol);
    
    if (viol > precision) {
      
      Jnlst()->Printf(Ipopt::J_MOREVECTOR, J_PROBLEM, "checkBounds(): variable %d out of bounds: %.6f [%g,%g] (diff %g)\n", 
		      i, val, domain_.lb (i), domain_.ub (i),
		      CoinMax (fabs (val - domain_.lb (i)), 
			       fabs (val - domain_.ub (i))));

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkBounds: variable %d out of bounds: %.6f [%g,%g] (diff %g)\n", 
			  i, val, domain_.lb (i), domain_.ub (i),
			  CoinMax (fabs (val - domain_.lb (i)), 
				   fabs (val - domain_.ub (i))));
      
      isFeas = false;
      if(stopAtFirstViol)
	break;
    }
  }
  return isFeas;
}

/************************************************************************/
// returns true iff value of all auxiliaries are within bounds
bool CouenneProblem::checkAux(const CouNumber *sol,
			      const bool stopAtFirstViol,
			      const double precision, double &maxViol) const {

  bool isFeas = true;
  for (int i=0; i<nVars(); i++) {

    exprVar *v = variables_ [i];

    if ((v -> Type         () != AUX) || 
	(v -> Multiplicity () <= 0)) 
      continue;

    if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {

      printf ("before check\n");

      double 
	vdb = (*(variables_ [i])) (), 
	fdb = (*(variables_ [i] -> Image ())) ();

      double
	del = 
	((v -> sign () == expression::AUX_GEQ) && (vdb >= fdb)) ? 0. : 
	((v -> sign () == expression::AUX_LEQ) && (vdb <= fdb)) ? 0. : 
	fabs (vdb - fdb);

      printf ("[%g,%g]\n", vdb, fdb);

      printf ("checking aux -- %+.12e = %+.12e [%+.12e] ", vdb, fdb, del);
      variables_ [i] -> print ();
      if(v->sign()== expression::AUX_GEQ) printf(" >= ");
      if(v->sign()== expression::AUX_LEQ) printf(" <= ");
      if(v->sign()== expression::AUX_EQ) printf(" := ");
      variables_ [i] -> Image () -> print (); printf ("\n");
    }
    
    double 
      vval = (*v) (),
      fval = (*(v -> Image ())) (),
      denom  = CoinMax (1., v -> Image () -> gradientNorm (X ()));
    
    // check if fval is a number (happens with e.g. w13 = w12/w5 and w5=0, see test/harker.nl)
    if (CoinIsnan (fval)) {
      fval = vval + 1.;
      denom = 1.;
    }
    
    if (fabs (fval) > COUENNE_INFINITY)
      fval = COUENNE_INFINITY;

    double
      delta = 
      ((v -> sign () == expression::AUX_GEQ) && (vval >= fval)) ? 0. : 
      ((v -> sign () == expression::AUX_LEQ) && (vval <= fval)) ? 0. : fabs (vval - fval),
      
      ratio = (CoinMax (1., fabs (vval)) / 
	       CoinMax (1., fabs (fval)));
    
    //printf ("checkinf --> v=%e f=%e den=%e ret=%e ratio=%e\n", vval, fval, denom, retval, ratio);
    
    double deldenom = delta/denom;
    if (maxViol <= deldenom) maxViol =  deldenom;

    if ((delta > 0.) &&
	(((ratio > 2.)  ||  // check delta > 0 to take into account semi-auxs
	  (ratio <  .5)) ||
	 ((delta /= denom) > CoinMin (COUENNE_EPS, feas_tolerance_)))) {

      Jnlst () -> Printf (Ipopt::J_MOREVECTOR, J_PROBLEM, "checkAux(): auxiliary %d violates tolerance %g by %g (deldenom: %g ratio %g)\n", i, feas_tolerance_, delta, deldenom, ratio);
      
      isFeas = false;

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkAux(): auxiliary %d violates tolerance %g by %g (deldenom: %g  ratio %g  COUENNE_EPS: %g)\n", 
			  i, feas_tolerance_, delta, deldenom, ratio, COUENNE_EPS);

      if (stopAtFirstViol)
	break;
    }
  }

  return isFeas;
}


/************************************************************************/
bool CouenneProblem::checkCons(const CouNumber *sol,
			       const bool stopAtFirstViol,
			       const double precision, double &maxViol) const {

  bool isFeas = true;
  for (int i=0; i<nCons(); i++) {
    
    CouenneConstraint *c = Con(i);
    
    CouNumber
      body = (*(c -> Body ())) (),
      lhs  = (*(c -> Lb   ())) (),
      rhs  = (*(c -> Ub   ())) ();
    
    double denomUb = 1 + CoinMax (fabs (body), fabs (rhs));
    double denomLb = 1 + CoinMax (fabs (body), fabs (lhs));
    double violUb = 0, violRelUb = 0, violAbsUb = 0;
    if(rhs < COUENNE_INFINITY) {
      violAbsUb = body - rhs;
      violRelUb = violAbsUb / denomUb;
      violUb = violAbsUb - precision * denomUb;

#ifdef FM_USE_ABS_VIOL_CONS
      if (maxViol <= violAbsUb) maxViol = violAbsUb;
#else
      if (maxViol <= violRelUb) maxViol = violRelUb;
#endif

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "violAbsUb: %12.10f  violRelUb: %12.10f  violUb: %12.10f maxViol: %12.10f\n", violAbsUb, violRelUb, violUb, maxViol);
    }

    double violLb = 0, violRelLb = 0, violAbsLb = 0;
    if(lhs > -COUENNE_INFINITY) {
      violAbsLb = - body + lhs;
      violRelLb = violAbsLb / denomLb;
      violLb = violAbsLb - precision * denomLb;

#ifdef FM_USE_ABS_VIOL_CONS
      if (maxViol <= violAbsLb) maxViol = violAbsLb;
#else
      if (maxViol <= violRelLb) maxViol = violRelLb;
#endif

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "violAbsLb: %12.10f  violRelLb: %12.10f  violLb: %12.10f maxViol: %12.10f\n", violAbsLb, violRelLb, violLb, maxViol);
    }

#ifdef FM_USE_ABS_VIOL_CONS
    if((violAbsUb > precision) || (violAbsLb > precision)) {
      if (Jnlst()->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
	
	printf ("CouenneProblem::checkCons(): constraint %d violated (lhs=%+e body=%+e rhs=%+e, absolute violation: %g)\n", i, lhs, body, rhs, CoinMax(violAbsUb, violAbsLb));
	c -> print ();
      }

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkCons(): constraint %d violated (lhs=%+e body=%+e rhs=%+e, absolute violation: %g)\n", i, lhs, body, rhs, CoinMax (violAbsUb, violAbsLb));

      isFeas = false;
      if(stopAtFirstViol) {
	break;
      }
    }
#else /* not FM_USE_ABS_VIOL_CONS */
    if((violUb > 0) || (violLb > 0)) {
      if (Jnlst()->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
	
	printf ("CouenneProblem::checkCons(): constraint %d violated (lhs=%+e body=%+e rhs=%+e, relative violation: %g)\n", i, lhs, body, rhs, CoinMax(violRelUb, violRelLb));
	c -> print ();
      }

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkCons(): constraint %d violated (lhs=%+e body=%+e rhs=%+e, relative violation: %g)\n", i, lhs, body, rhs, CoinMax (violRelUb, violRelLb));

      isFeas = false;
      if(stopAtFirstViol) {
	break;
      }
    }
#endif /* not FM_USE_ABS_VIOL_CONS */
  }
  return(isFeas);
}

/************************************************************************/
bool CouenneProblem::checkNLP2(const double *solution, 
			       const double obj, const bool careAboutObj,
			       const bool stopAtFirstViol,
			       const bool checkAll,
			       const double precision) const {

  if (careAboutObj && stopAtFirstViol) {
    printf("CouenneProblem::checkNLP2(): ### ERROR: careAboutObj: true and stopAtFirstViol: true are incompatible\n");
    exit(1);
  }

  const std::vector<int> listInt = getRecordBestSol()->getListInt();
  bool isFeas = false;
  double maxViolCouSol = 0;
  double maxViolRecSol = 0;

  if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {

    const bool *initIsInt = getRecordBestSol()->getInitIsInt();
    printf("Integrality:\n");
    for (int i=0; i<nVars(); i++) {

      if (variables_ [i] -> Multiplicity () <= 0) 
	continue;

      if(initIsInt[i]) {
	printf(" %d", i);
      }
    }
    printf("\n");

    printf("VAR:\n");
    for (int i=0; i<nVars(); i++) {

      if (variables_ [i] -> Multiplicity () <= 0) 
	continue;
      exprVar *v = variables_ [i];
      if(	(v -> Type () == VAR)) {
	printf(" %d", i);
      }
    }
    printf("\n");

    printf("AUX:\n");
    for (int i=0; i<nVars(); i++) {

      if (variables_ [i] -> Multiplicity () <= 0) 
	continue;
      exprVar *v = variables_ [i];
      if(	(v -> Type () == AUX)) {
	printf(" %d", i);
      }
    }
    printf("\n");

    printf("mult 0:\n");
    for (int i=0; i<nVars(); i++) {

      if (variables_ [i] -> Multiplicity () <= 0) { 
	printf(" %d", i);
      }
    }
    printf("\n");
  }

  if ((Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM))) {
    printf ("checking solution:\n");
    for (int i=0; i<nOrigVars_ - ndefined_; i++)
      printf ("%.12e ", solution [i]);
    printf ("\nCouenneProblem::checkNLP2(): Start checking recomputed_solution\n");
  }

  // check integrality of integer constrained variables

  int from = 0;
  bool isFeasRec = checkInt(solution, from, nOrigVars_ - ndefined_, listInt,
			    false, stopAtFirstViol,
			    precision, maxViolRecSol);
  bool isFeasCou = isFeasRec;
  maxViolCouSol = maxViolRecSol;

  if (stopAtFirstViol && !isFeasRec && Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) 
    printf("CouenneProblem::checkNLP2(): recomputed_solution is infeasible (some orig vars not integer feasible; violation: %12.10g)\n", maxViolRecSol);

#ifdef CHECK
  if(getRecordBestSol()->getCardInitDom() != nVars()) {
    printf("CouenneProblem::checkNLP2(): ### ERROR: cardInitDom: %d  nVars: %d\n", getRecordBestSol()->getCardInitDom(), nVars());
    exit(1);
  }
  if(getInitDomLb() == NULL) {
    printf("CouenneProblem::checkNLP2(): ### WARNING: initDomLb == NULL\n");
  }
  if(getInitDomUb() == NULL) {
    printf("CouenneProblem::checkNLP2(): ### WARNING: initDomUb == NULL\n");
  }
#endif

  // install NL solution candidate and original bounds in evaluation structure
  // bounds are important so that getAuxs below works properly
  domain_.push(nVars(), solution, getRecordBestSol()->getInitDomLb(), 
	       getRecordBestSol()->getInitDomUb(), false);

  CouNumber *couRecSol = new CouNumber[nVars()];
  CoinCopyN (solution, nOrigVars_ - ndefined_, couRecSol);
  getAuxs(couRecSol);
  //CoinCopyN (solution, nOrigVars_, couRecSol);

  domain_.pop (); // getting rid of current domain now as won't be used again

  // install couRecSol in evaluation structure
  domain_.push(nVars(), couRecSol, 
	       getRecordBestSol()->getInitDomLb(), 
	       getRecordBestSol()->getInitDomUb(), false);

  if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {
    printf ("checkNLP2(): recomputed_solution: %d vars -------------------\n", domain_.current () -> Dimension ());
    double maxDelta = 0;
    for (int i=0; i<nVars (); i++) {
      exprVar *v = variables_ [i];
      if (v -> Multiplicity () <= 0) 
        continue;
      if(i < nOrigVars_ - ndefined_) {
        double soli = solution[i];
        double domi = domain_.x (i);
        double domlbi = domain_.lb (i);
        double domubi = domain_.ub (i);
        printf ("%4d %+e %+e [%+e %+e] %+e\n", i, solution[i], domain_.x (i), domain_.lb (i), domain_.ub (i), solution[i] - domain_.x (i));
        if (maxDelta < fabs(solution[i] - domain_.x(i))) maxDelta = fabs (solution[i] - domain_.x(i));
      }
      else {
        if(v -> Type() == AUX) {
          printf ("%4d  ------     %+e [%+e %+e]\n", i, domain_.x (i), domain_.lb (i), domain_.ub (i));
        }
      }
      printf ("maxDelta: %.12g\n", maxDelta);
    }
  }

  if(checkAll) {
    if(!stopAtFirstViol || isFeasRec) {
      bool isFeasInt = checkInt(couRecSol, 0, nVars(), listInt,
				false, stopAtFirstViol,
				precision, maxViolRecSol);

      if (!isFeasInt) {
	Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): recomputed_solution is infeasible (some aux vars not integer feasible; violation: %12.10g)\n", maxViolRecSol);
	
	isFeasRec = false;
      }
    }
  }

  double objRecSol = checkObj(couRecSol, precision);
  double objCouSol = 0;

  if(checkAll) {
    if(!stopAtFirstViol || isFeasRec) {
      bool isFeasBnd = checkBounds(couRecSol, stopAtFirstViol, 
				   precision, maxViolRecSol);

      if(!isFeasBnd) {
	
	Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): recomputed_solution is infeasible (violated bounds; violation: %12.10g)\n", maxViolRecSol);

	isFeasRec = false;
      }
    }

    if(!stopAtFirstViol || isFeasRec) {
      bool isFeasAux = checkAux(couRecSol, stopAtFirstViol, 
				precision, maxViolRecSol);

      if(!isFeasAux) {
	
	Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): recomputed_solution is infeasible (violated Aux; violation: %12.10g)\n", maxViolRecSol);
	
	isFeasRec = false;
      }
    }
  }

  if(!stopAtFirstViol || isFeasRec) {
    bool isFeasCons = checkCons(couRecSol, stopAtFirstViol, 
				precision, maxViolRecSol);

    if(!isFeasCons) {

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): recomputed_solution is infeasible (violated constraint; violation: %12.10g)\n", maxViolRecSol);

      isFeasRec = false;
    }
  }

  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): end check recomputed_solution (maxViol: %12.10g)\n", maxViolRecSol);

  double objErrorRecSol = objRecSol - obj;
  if(!careAboutObj)
    objErrorRecSol = 0;

  CouNumber *couSol = new CouNumber[nVars()];
  bool useRecSol = false;
  if(isFeasRec && (objErrorRecSol < precision)) {
    useRecSol = true;
    isFeas = true;
  }
  else {

    if(checkAll) { // otherwise duplicates above calculations

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): Start checking solution (maxViol: %g)\n", maxViolCouSol);

      CoinCopyN(solution, nVars(), couSol);
      restoreUnusedOriginals(couSol);
      domain_.push(nVars(), couSol, getRecordBestSol()->getInitDomLb(), 
		   getRecordBestSol()->getInitDomUb(), false);

      if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {
	printf ("checkNLP2(): solution: %d vars -------------------\n", domain_.current () -> Dimension ());	
	double maxDelta = 0;
	for (int i=0; i<domain_.current()->Dimension(); i++) {
	  //printf ("%4d %.12g %.12g [%.12g %.12g]\n", i, solution[i], domain_.x(i), domain_.lb(i), domain_.ub(i));
	  printf ("%4d %+e %+e [%+e %+e] %+e\n",          i, solution[i], domain_.x (i), domain_.lb (i), domain_.ub (i), solution[i] - domain_.x (i));      
	  if (fabs (solution[i] - domain_.x(i)) > maxDelta)
	    maxDelta = fabs(solution[i] - domain_.x(i));
	}
	printf("maxDelta: %.12g\n", maxDelta);
      }

      if(!stopAtFirstViol || isFeasCou) {
	bool isFeasInt = checkInt(couSol, nOrigVars_ - ndefined_, nVars(), listInt,
				  false, stopAtFirstViol,
				  precision, maxViolCouSol);
	if(!isFeasInt) {

	  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): solution is infeasible (some aux vars not integer feasible; violation: %12.10g)\n", maxViolCouSol);
	  
	  isFeasCou = false;
	}
      }
      
      objCouSol = checkObj(couSol, precision);
      
      if(!stopAtFirstViol || isFeasCou) {
	bool isFeasCouBnd = checkBounds(couSol, stopAtFirstViol, 
					precision, maxViolCouSol);
	if(!isFeasCouBnd) {
	  
	  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): solution is infeasible (some bounds are violated; violation: %12.10g)\n", maxViolCouSol);
	  
	  isFeasCou = false;
	}
      }
      
      if(!stopAtFirstViol || isFeasCou) {
	bool isFeasCouAux = checkAux(couSol, stopAtFirstViol, 
				     precision, maxViolCouSol);
	if(!isFeasCouAux) {
	  
	  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): solution is infeasible (violated Aux; violation: %12.10g)\n", maxViolCouSol);
	  
	  isFeasCou = false;
	}
      }
      
      if(!stopAtFirstViol || isFeasCou) {
	bool isFeasCouCons = checkCons(couSol, stopAtFirstViol, 
				       precision, maxViolCouSol);
	if(!isFeasCouCons) {
	  
	  Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): solution is infeasible (violated constraint; violation: %12.10g)\n", maxViolCouSol);
	  
	  isFeasCou = false;
	}
      }
    
      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): end check solution (maxViol: %12.10g)\n", maxViolCouSol);
    
      double objErrorCouSol = objCouSol - obj;
      if(!careAboutObj) {
	objErrorCouSol = 0;
      }
      double delObj = objErrorCouSol - objErrorRecSol;
      double delViol = maxViolCouSol - maxViolRecSol;

      if(isFeasRec) {
	if(isFeasCou) {
	  // careAboutObj must be true
	  if(delObj > 0) {
	    useRecSol = true;
	  }
	  else {
	    useRecSol = false;            
	  }
	  isFeas = true;
	}
	else { /* isFeasRec == true and isFeasCou == false */
	  useRecSol = true;
	  isFeas = true;
	}
      }
      else { /* isFeasRec == false */
	if(isFeasCou) {
	  useRecSol = false;            
	  isFeas = true;
	}
	else { /* isFeasRec == false and isFeasCou == false */
	  isFeas = false;
	  if(careAboutObj) {
	    if(fabs(delViol) < 10 * precision) {
	      useRecSol = (delObj <= 0 ? false : true);
	    }
	    else {
	      if(fabs(delObj) < 10 * precision) {
		useRecSol = (delViol > 0 ? false : true);
	      }
	      else {
		double ratObj = fabs(delObj)/(1+fabs(obj));
		if(ratObj < fabs(delViol)) {
		  useRecSol = (delViol > 0 ? false : true);
		}
		else {
		  useRecSol = (delObj > 0 ? false : true);
		}
	      }
	    }
	  }
	  else {
	    if(delViol < 0) {
	      useRecSol = false;
	    }
	    else {
	      useRecSol = true;
	    }
	  }
	  useRecSol = true;
	}
      }
      domain_.pop (); // pop couSol
    } 
  }
  
  double maxViol = 0;
  
  if(!stopAtFirstViol || isFeas) {
    if(useRecSol) {
      recBSol->setModSol(couRecSol, nVars(), objRecSol, maxViolRecSol);
      maxViol = maxViolRecSol;
      
      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): select recomputed_solution (maxViol: %12.10g)\n", maxViol);      
    }
    else {
      recBSol -> setModSol(couSol, nVars(), objCouSol, maxViolCouSol);
      maxViol = maxViolCouSol;
      
      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "CouenneProblem::checkNLP2(): select solution (maxViol: %12.10g)\n", maxViol);      
    }
  }

  if (Jnlst () -> ProduceOutput (Ipopt::J_ALL, J_PROBLEM)) {
    if(isFeas) {

      // printf ("Solution: [");
      // for (int i = 0; i < nVars (); ++i)
      //   printf ("%g ", couSol [i]);
      // printf ("]\n");

      printf ("checkNLP2(): RETURN: feasible (maxViol: %g)\n", maxViol);
    }
    else {
      printf ("checkNLP2(): RETURN: recomputed_solution and solution are infeasible\n");
      if(!stopAtFirstViol) {
	printf("(maxViol: %g)\n", maxViol);
      }
      else {
	printf("\n");
      }
    }
  }

  delete[] couSol;
  delete[] couRecSol;

  domain_.pop (); // pop bounds
    
  return isFeas;
}

// comprehensive method to call one of the two variants
bool CouenneProblem::checkNLP0 (const double *solution, 
			       double &obj,
			       bool recompute_obj, 
			       const bool careAboutObj,
			       const bool stopAtFirstViol,
			       const bool checkAll,
			       const double precision) const {

  bool retval;

#ifdef FM_CHECKNLP2

  retval = checkNLP2 (solution,
		      obj,
		      careAboutObj,
		      stopAtFirstViol,
		      checkAll,
		      (precision < 0.) ? feas_tolerance_ : precision);

  if (retval)
    obj = getRecordBestSol () -> getModSolVal ();

#else 

  retval = checkNLP1 (solution, obj, recompute_obj);

#endif

  return retval;
}
