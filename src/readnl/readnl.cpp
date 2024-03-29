/*
 *
 * Name:    readnl.cpp
 * Author:  Pietro Belotti
 * Purpose: define a reader for .nl files. Adapted from ampl2ev3 by L. Liberti and S. Galli
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

// include CoinTime.hpp before ASL defines filename for https://github.com/coin-or-tools/ThirdParty-ASL/issues/4
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "asl.h"

// Added to avoid compiler issues with gcc 6.2.1 --- see https://github.com/JuliaOpt/CoinOptServices.jl/issues/27#issuecomment-290960312
#undef strtod

#include "nlp.h"
#include "getstub.h"
#include "opcode.hd"

#include "CouenneProblem.hpp"

#if   defined HAVE_CSTDINT
#include "cstdint"
#elif defined HAVE_STDINT_H
#include "stdint.h"
#endif

#include "CouenneTypes.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprGroup.hpp"

#define OBJ_DE    ((const ASL_fg *) asl) -> I.obj_de_
#define VAR_E     ((const ASL_fg *) asl) -> I.var_e_
#define CON_DE    ((const ASL_fg *) asl) -> I.con_de_
#define OBJ_sense ((const ASL_fg *) asl) -> i.objtype_

#define THRESHOLD_OUTPUT_READNL 10000

//#define DEBUG

using namespace Couenne;

// check if an expression is a null pointer or equals zero
inline bool is_expr_zero (expr* e)
{return ((e==NULL) || ((((Intcast (e->op)) == OPNUM) &&
			  (fabs (((expr_n *)e) -> v)  < COUENNE_EPS)
			  //  && (fabs (e -> dL) < COUENNE_EPS)
			  // *** CHECK THIS! dL is the derivative
			  )));}

void createCommonExpr (CouenneProblem *p, const ASL *asl, int i, int which);

// Reads a MINLP from an AMPL .nl file through the ASL methods
int CouenneProblem::readnl (const ASL *asl) {

  problemName_ = filename;

  // number of defined variables (aka common expressions)
  ndefined_ = como + comc + comb + como1 + comc1;

  // see "hooking your solver to AMPL", by David M. Gay, tables 3, 4, and 5

  // nonlinear in both objectives and constraints
  if (nlvb >= 0) {
    for (int i = 0; i < nlvb - nlvbi; i++) addVariable (false, &domain_);
    for (int i = 0; i < nlvbi;        i++) addVariable (true,  &domain_);
  }

  // nonlinear in either objectives or constraints
  if (nlvo > nlvc) {
    for (int i = 0; i < nlvc - (nlvb + nlvci); i++) addVariable (false, &domain_);
    for (int i = 0; i < nlvci;                 i++) addVariable (true,  &domain_);
    for (int i = 0; i < nlvo - (nlvc + nlvoi); i++) addVariable (false, &domain_);
    for (int i = 0; i < nlvoi;                 i++) addVariable (true,  &domain_);
  } else {
    for (int i = 0; i < nlvo - (nlvb + nlvoi); i++) addVariable (false, &domain_);
    for (int i = 0; i < nlvoi;                 i++) addVariable (true,  &domain_);
    for (int i = 0; i < nlvc - (nlvo + nlvci); i++) addVariable (false, &domain_);
    for (int i = 0; i < nlvci;                 i++) addVariable (true,  &domain_);
  }

  for (int i = 0; i < nwv; i++)                                  addVariable(false, &domain_);//arc
  for (int i = n_var - (CoinMax (nlvc,nlvo) +niv+nbv+nwv); i--;) addVariable(false, &domain_);//other
  for (int i = 0; i < nbv; i++)                                  addVariable(true,  &domain_);//binary
  for (int i = 0; i < niv; i++)                                  addVariable(true,  &domain_);//int.

  // add space for common expressions
  for (int i = ndefined_; i--;)                                  addVariable(false, &domain_);

  double now = CoinCpuTime ();

  if (nVars () > THRESHOLD_OUTPUT_READNL) {
    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "Reading problem: ");
    fflush(stdout);
  }

  // common expressions (or defined variables) ///////////////////////////////////////

#ifdef DEBUG
  printf ("ndefd_  = %d\n", ndefined_);
  printf ("tot var = %d\n", variables_ . size ());
  printf ("c_vars_ = %d\n", ((const ASL_fg *) asl) -> i.c_vars_ );
  printf ("comb_ = %d\n",   ((const ASL_fg *) asl) -> i.comb_  );
  printf ("combc_ = %d\n",  ((const ASL_fg *) asl) -> i.combc_ );
  printf ("comc1_ = %d\n",  ((const ASL_fg *) asl) -> i.comc1_ );
  printf ("comc_ = %d\n",   ((const ASL_fg *) asl) -> i.comc_  );
  printf ("como1_ = %d\n",  ((const ASL_fg *) asl) -> i.como1_ );
  printf ("como_ = %d\n",   ((const ASL_fg *) asl) -> i.como_  );
#endif

  // Each has a linear and a nonlinear part (thanks to Dominique
  // Orban: http://www.gerad.ca/~orban/drampl/def-vars.html)

  for (int i = 0; i < como + comc + comb; i++)
    createCommonExpr (this, asl, i, 0);

  for (int i = 0; i < como1 + comc1; i++)
    createCommonExpr (this, asl, i, 1);

  //commonexprs_ . erase (commonexprs_ . begin (), commonexprs_ . end ());

  // objective functions /////////////////////////////////////////////////////////////

  if (n_obj == 0) {

    // strange, no objective function. Add one equal to zero

    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "Couenne: warning, no objective function found\nAdded fictitious function f(x)=0\n");
    addObjective (new exprConst (0.), "min");
  }

  for (int i = 0; i < n_obj; i++) {

    ////////////////////////////////////////////////
    int nterms = 0;

    // count nonzero terms in linear part

    for (ograd *objgrad = Ograd [i];
	 objgrad;
	 objgrad = objgrad -> next)
      if (fabs (objgrad -> coef) > COUENNE_EPS)
	nterms++;

    expression
      *body,
      *nl = Simplified (nl2e (OBJ_DE [i] . e, asl));
    //   *nls = nl -> simplify ();

    // if (nls) {
    //   delete nl;
    //   nl = nls;
    // }

    if (nterms) { // have linear terms

      int       *indexL = new int       [nterms+1];
      CouNumber *coeff  = new CouNumber [nterms];

      for (ograd *objgrad = Ograd [i]; objgrad; objgrad = objgrad -> next)
	if (fabs (objgrad -> coef) > COUENNE_EPS) {

	  *indexL++ = objgrad -> varno;
	  *coeff++  = objgrad -> coef;
	}

      *indexL = -1;

      indexL -= nterms;
      coeff  -= nterms;

      std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      indcoe2vector (indexL, coeff, lcoeff);

      if (nl -> code () == COU_EXPRSUM) {
	body = exprGroup::genExprGroup (0., lcoeff, nl -> ArgList (), nl -> nArgs ());
	// delete node without deleting children (they are now in body)
	nl -> ArgList (NULL);
	delete nl;
      }
      else {

	expression **nll = new expression * [1];

	*nll = nl;

	// apparently, objconst (i) is included in the obj expression
	body = exprGroup::genExprGroup (0., lcoeff, nll, 1);
	//body = new exprGroup (objconst (i), indexL, coeff, nll, 1);
      }

      delete [] indexL;
      delete [] coeff;

    } else
      // apparently, objconst (i) is included in the obj expression
      body = nl;
      //if (fabs (objconst (i) > COUENNE_EPS))
      //body = new exprSum (nl, new exprConst (objconst (i)));
      //else

    ///////////////////////////////////////////////////

    expression *subst = Simplified (body); //  -> simplify ();

    // if (subst) {
    //   delete body; // VALGRIND
    //   body = subst;
    // }

    // ThirdParty/ASL/solvers/asl.h, line 336: 0 is minimization, 1 is maximization
    addObjective (body, (OBJ_sense [i] == 0) ? "min" : "max");
  }

  // constraints ///////////////////////////////////////////////////////////////////

  int *nterms = new int [n_con];

  // allocate space for argument list of all constraints' summations
  // of linear and nonlinear terms

  // init array with # terms of each constraint
  for (int i = n_con; i--;)
    *nterms++ = 0;
  nterms -= n_con;

  cgrad *congrad;

  // count all linear terms
  if (A_colstarts && A_vals)         // Constraints' linear info is stored in A_vals
    for (int j = A_colstarts [n_var]; j--;) {

      real coeff = A_vals [j];

      if (fabs (coeff) > COUENNE_EPS)
	nterms [A_rownos [j]] ++;
    }
  else {                             // Constraints' linear info is stored in Cgrad
    for (int i = 0; i < n_con; i++)
      for (congrad = Cgrad [i];
	   congrad;
	   congrad = congrad -> next)
	if (fabs (congrad -> coef) > COUENNE_EPS)
	  nterms [i] ++;
  }


  // vectors of the linear part
  CouNumber **coeff  = new CouNumber * [n_con];
  int       **indexL = new int       * [n_con];

  for (int i = n_con; i--;)
    *indexL++ = NULL;

  indexL -= n_con;


  // set linear terms

  if (A_colstarts && A_vals)         // Constraints' linear info is stored in A_vals
    for (int j = 0; j < n_var; j++)
      for (int i = A_colstarts [j], k = A_colstarts [j+1] - i; k--; i++) {

	int rowno = A_rownos [i],
	    nt    = nterms [rowno] --;

	CouNumber **cline = coeff  + rowno;
	int       **iline = indexL + rowno;

	if (*iline==NULL) {
	  *cline = new CouNumber [nt];
	  *iline = new int       [nt+1];
	  (*iline) [nt] = -1;
	}

	(*cline) [--nt] = A_vals [i];
	(*iline)   [nt] = j;

      }
  else {                             // Constraints' linear info is stored in Cgrad
    for (int i=0; i < n_con; i++) {

      int nt = nterms [i];

      CouNumber **cline = coeff + i;
      int       **iline = indexL + i;

      *cline = new CouNumber [nt];
      *iline = new int       [nt+1];
      (*iline) [nt] = -1;

      for (congrad = Cgrad [i]; congrad; congrad = congrad -> next)
	if (fabs (congrad -> coef) > COUENNE_EPS) {
	  (*cline) [--nt] = congrad -> coef;
	  (*iline)   [nt] = congrad -> varno;
	}
    }
  }

  // set constraints' bound and sign and store nonlinear part ///////////////////////////////

  for (int i = 0; i < n_con; i++) {

    enum con_sign sign;
    double lb, ub;

    if (Urhsx) {
      lb = LUrhs [i];
      ub = Urhsx [i];
    } else {
      int j = 2*i;
      lb = LUrhs [j];
      ub = LUrhs [j+1];
    }

    // set constraint sign
    if (lb > -COUENNE_INFINITY /*negInfinity*/)
      if (ub < COUENNE_INFINITY /*Infinity*/) sign = COUENNE_RNG;
      else               sign = COUENNE_GE;
    else                 sign = COUENNE_LE;

    // this is an equality constraint
    if (fabs (lb - ub) < COUENNE_EPS)
      sign = COUENNE_EQ;

    expression
       *body,
      **nll = new expression * [1],
       *nls;

    *nll = Simplified (nl2e (CON_DE [i] . e, asl));

    // nls = (*nll) -> simplify ();

    // if (nls) {
    //   delete *nll;
    //   *nll = nls;
    // }

    if (indexL [i] && (*(indexL [i]) >= 0)) {

      int code = (*nll) -> code ();

      std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      indcoe2vector (indexL [i], coeff [i], lcoeff);

      /*std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      for (int i=0, *ind = indexL; *ind >= 0; *ind++, i++)
      lcoeff.push_back (std::pair <exprVar *, CouNumber> (Var (*ind), coeff [i]));*/

      if ((code == COU_EXPRSUM) ||
	  (code == COU_EXPRGROUP)) {

	body    = exprGroup::genExprGroup (0., lcoeff, (*nll) -> ArgList (), (*nll) -> nArgs ());
	// delete node without deleting children (they are now in body)
	(*nll) -> ArgList (NULL);
	delete *nll;
	delete [] nll;
      }
      else body = exprGroup::genExprGroup (0., lcoeff, nll, 1);
    }
    else {
      body = *nll;
      delete [] nll;
    }

    expression *subst = Simplified (body); // or body->simplify() ?

    //  -> simplify ();
    // if (subst) {
    //   delete body; // VALGRIND
    //   body = subst;
    // }

    // add them (and set lower-upper bound)
    if ((lb < -COUENNE_INFINITY /*negInfinity*/) &&
        (ub > COUENNE_INFINITY /*Infinity*/)) {

      printf ("Free constraint %d ignored\n", i);

    } else

      // add them (and set lower-upper bound)
      switch (sign) {

      case COUENNE_EQ:  addEQConstraint  (body, new exprConst (ub)); break;
      case COUENNE_LE:  addLEConstraint  (body, new exprConst (ub)); break;
      case COUENNE_GE:  addGEConstraint  (body, new exprConst (lb)); break;
      case COUENNE_RNG: addRNGConstraint (body, new exprConst (lb),
                                                new exprConst (ub)); break;
      default: printf ("Could not recognize constraint\n"); return -1;
      }

    delete [] indexL [i];
    delete [] coeff  [i];
  }

  delete [] indexL;
  delete [] coeff;

  // create room for problem's variables and bounds
  CouNumber
    *x  = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber)),
    *lb = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber)),
    *ub = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber));

  for (int i = n_var + ndefined_; i--;) {
    x  [i] =  0.;
    lb [i] = -COUENNE_INFINITY;
    ub [i] =  COUENNE_INFINITY;
  }

  domain_.push (n_var + ndefined_, x, lb, ub);

  free (x); free (lb); free (ub);

  // lower and upper bounds ///////////////////////////////////////////////////////////////

  if (LUv) {

    real *Uvx_copy = Uvx;

    if (!Uvx_copy)
      for (int i=0; i<n_var; i++) {

	int j = 2*i;

	Lb (i) = LUv [j];
	Ub (i) = LUv [j+1];
      }
    else
      for (int i=n_var; i--;) {
	Lb (i) = LUv [i];
	Ub (i) = Uvx_copy [i];
      }

  } else
    for (int i=n_var; i--;) {
      Lb (i) = - COUENNE_INFINITY;
      Ub (i) =   COUENNE_INFINITY;
    }

  // initial x ////////////////////////////////////////////////////////////////////

  for (int i=n_var; i--;)

    if (X0 && havex0 [i]) X (i) = X0 [i];

    else {

      CouNumber x, l = Lb (i), u = Ub (i);

      if      (l < - COUENNE_INFINITY)
	if    (u >   COUENNE_INFINITY)  x = 0.;
	else                            x = u;
      else if (u >   COUENNE_INFINITY)  x = l;
      else                              x = 0.5 * (l+u);

      X (i) = x;
    }

  for (int i=n_var; i<ndefined_; i++) {

    X  (i) =  0.;
    Lb (i) = -COUENNE_INFINITY;
    Ub (i) =  COUENNE_INFINITY;
  }

  delete [] nterms;

  if (nVars () > THRESHOLD_OUTPUT_READNL) {
    jnlst_ -> Printf (Ipopt::J_ERROR, J_COUENNE, "%.1f seconds\n", CoinCpuTime () - now);
    fflush(stdout);
  }

  return 0;
}


//
// Common code for common expressions (aka defined variables)
//

void createCommonExpr (CouenneProblem *p, const ASL *asl, int i, int which) {

  struct cexp  *common  = ((const ASL_fg *) asl) -> I.cexps_  + i;
  struct cexp1 *common1 = ((const ASL_fg *) asl) -> I.cexps1_ + i;

  expression
    *nle = Simplified (p -> nl2e (which ? common1 -> e : common -> e, asl));

  //   *nls = nle -> simplify ();

  // if (nls) {
  //   delete nle;
  //   nle = nls;
  // }

#ifdef DEBUG
  printf ("cexp1 %d [%d]: ", i, p -> Variables () . size ()); nle -> print ();  printf (" ||| ");
#endif

  int nlin = which ? common1 -> nlin : common -> nlin;  // Number of linear terms

  if (nlin > 0) {

    linpart *L = which ? common1 -> L : common -> L;

    std::vector <std::pair <exprVar *, CouNumber> > lcoeff;

    for (int j = 0; j < nlin; j++) {
      //vp = (expr_v *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));

      int indVar = ((uintptr_t) (L [j].v.rp) - (uintptr_t) VAR_E) / sizeof (expr_v);
      CouNumber coeff = L [j]. fac;

      lcoeff.push_back (std::pair <exprVar *, CouNumber> (p -> Var (indVar), coeff));

#ifdef DEBUG
      Printf( " %+g x_%d (%-3d)", L [j]. fac, indVar,
	      (expr_v *) (L [j].v.rp) - VAR_E //((const ASL_fg *) asl) -> I.cexps_
	      //L [j]. v.i
	      );
#endif
    }

    expression **al = new expression * [1];
    *al = nle;

    expression *eg;

    if (lcoeff.size  () == 1 &&
	nle -> Type  () == CONST &&
	nle -> Value () == 0.) {

      CouNumber coeff = lcoeff [0].second;

      if      (coeff ==  1.) eg =                                     new exprClone (lcoeff [0].first);
      else if (coeff == -1.) eg = new exprOpp                        (new exprClone (lcoeff [0].first));
      else                   eg = new exprMul (new exprConst (coeff), new exprClone (lcoeff [0].first));
    } else                   eg = exprGroup::genExprGroup (0, lcoeff, al, 1);

    int indVar = p -> nVars () - p -> nDefVars () + p -> commonExprs () . size ();

    if (eg -> Index () != indVar) {

      expression *body = new exprSub (eg, new exprClone (p -> Var (indVar)));
      p -> addEQConstraint (body, new exprConst (0.));
    }

    p -> commonExprs () . push_back (new exprClone (eg));
  }
  else {

    int indVar = p -> nVars () - p -> nDefVars () + p -> commonExprs () . size ();

    if (nle -> Index () != indVar) {

      expression *body = new exprSub (nle, new exprClone (p -> Var (indVar)));
      p -> addEQConstraint (body, new exprConst (0.));
    }

    p -> commonExprs () . push_back (new exprClone (nle));
  }

#ifdef DEBUG
  printf ("\n");
#endif
  //    addAuxiliary (nl2e (((const ASL_fg *) asl) -> I.cexps1_ [i] . e, asl));
}
