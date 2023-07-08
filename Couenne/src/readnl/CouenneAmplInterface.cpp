// $Id$
//
// (C) Copyright XXX 2009
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pietro Belotti, Lehigh University
// Stefan Vigerske, Humboldt University
//
// Date : 07/18/2009

#include "CouenneConfig.h"
#include "CouenneAmplInterface.hpp"
#include "CoinPragma.hpp"

#include <cstdlib>

#if   defined HAVE_CSTDINT
#include <cstdint>
#elif defined HAVE_STDINT_H
#include <stdint.h>
#endif

#include <string>

#include "BonAmplTMINLP.hpp"
#include "BonCbc.hpp"

#include "CouenneProblem.hpp"
#include "CouenneTypes.hpp"

#include "CouenneExprClone.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneExprAbs.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprInv.hpp"
#include "CouenneExprSin.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprLog.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprCos.hpp"
#include "CouenneExprExp.hpp"

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "opcode.hd"

// get ASL op. code relative to function pointer passed as parameter 
size_t getOperator (efunc *);

#define OBJ_DE    ((const ASL_fg *) asl) -> I.obj_de_
#define VAR_E     ((const ASL_fg *) asl) -> I.var_e_
#define CON_DE    ((const ASL_fg *) asl) -> I.con_de_
#define OBJ_sense ((const ASL_fg *) asl) -> i.objtype_

#include "r_opn.hd" /* for N_OPS */

static fint timing = 0;

static
keyword keywds[] = { // must be alphabetical
   KW(const_cast<char*>("timing"), L_val, &timing, const_cast<char*>("display timings for the run")),
};

static
Option_Info Oinfo = { const_cast<char*>("testampl"), const_cast<char*>("ANALYSIS TEST"),
		      const_cast<char*>("concert_options"), keywds, nkeywds, 0, const_cast<char*>("ANALYSIS TEST") };


// (C++) code starts here ///////////////////////////////////////////////////////////////////////////

using namespace Couenne;

void CouenneAmplInterface::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions) {
	roptions->AddStringOption1("nlfile", "name of an ampl .nl file to get the problem from", "", "*", "name of .nl file");
}

CouenneAmplInterface::~CouenneAmplInterface() {
	delete problem;
	
	if (asl) {
	  delete[] X0;
	  delete[] havex0;
	  delete[] pi0;
	  delete[] havepi0;		
		ASL_free(&asl);
	}
}

// create an AMPL problem by using ASL interface to the .nl file
CouenneProblem* CouenneAmplInterface::getCouenneProblem() {
	if (problem)
		return problem;

  if (!readASLfg())
  	return NULL;

  problem = new CouenneProblem;

  if (!readnl()) {
  	delete problem;
  	problem = NULL;
  	return NULL;
  }
  
  return problem;
}

Ipopt::SmartPtr<Bonmin::TMINLP> CouenneAmplInterface::getTMINLP() {
	if (IsValid(tminlp))
		return tminlp;
	
	if (IsNull(roptions)) {
		jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_INITIALIZATION, "Error: Need registered options to create AmplTMINLP object!\n");
		return NULL;
	}
	
	std::string nlfile;
	options->GetStringValue("nlfile", nlfile, "");
	char** argv = new char*[3];
	argv[0] = const_cast<char*>("dummy");
	argv[1] = strdup(nlfile.c_str());
	argv[2] = NULL;
	tminlp = new Bonmin::AmplTMINLP(GetRawPtr(jnlst), roptions, options, argv);
	
	free(argv[1]);
	delete[] argv;
	
	return tminlp;
}

bool CouenneAmplInterface::writeSolution(Bonmin::Bab& bab) {
	const char* message;

	//TODO setup a nicer message
	if (bab.bestSolution()) {
		message = "Couenne found a solution.\n";
	} else {
		message = "Couenne could not found a solution.\n";
	}

	write_sol(const_cast<char*>(message), const_cast<double*>(bab.bestSolution()), NULL, NULL);
	
	return true;
}

bool CouenneAmplInterface::readASLfg() {
	assert(asl == NULL);
	
  std::string nlfile;
  options->GetStringValue("nlfile", nlfile, "");
  
  if (nlfile == "")
  	return false;

  char** argv = new char*[3];
	argv[0] = const_cast<char*>("dummy");
	argv[1] = strdup(nlfile.c_str());
	argv[2] = NULL;
	
  // Create the ASL structure
  asl = (ASL*) ASL_alloc (ASL_read_fg);

  char* stub = getstub (&argv, &Oinfo);

  // Although very intuitive, we shall explain why the second argument
  // is passed with a minus sign: it is to tell the ASL to retrieve
  // the nonlinear information too.
  FILE* nl = jac0dim (stub, - (fint) strlen (stub));

  // Set options in the asl structure
  want_xpi0 = 1 | 2;  // allocate initial values for primal and dual if available
  obj_no = 0;         // always want to work with the first (and only?) objective

  // allocate space for initial values
  X0      = new real [n_var];
  havex0  = new char [n_var];
  pi0     = new real [n_con];
  havepi0 = new char [n_con];

  // read the rest of the nl file
  fg_read (nl, ASL_return_read_err | ASL_findgroups);
  
  //FIXME freeing argv and argv[1] gives segfault !!!
//  free(argv[1]);
//  delete[] argv;

  return true;
}

// check if an expression is a null pointer or equals zero
//static inline bool is_expr_zero (expr* e) {
//	return ((e==NULL) || ((((Intcast (e->op)) == OPNUM) && 
//			  (fabs (((expr_n *)e) -> v)  < COUENNE_EPS) 
//			  //  && (fabs (e -> dL) < COUENNE_EPS)
//			  // *** CHECK THIS! dL is the derivative
//			  )));
//} 


// Reads a MINLP from an AMPL .nl file through the ASL methods
bool CouenneAmplInterface::readnl() {

  std::string nlfile;
  options->GetStringValue("nlfile", nlfile, "");
  problem -> setProblemName (nlfile);

  // number of defined variables (aka common expressions)
  problem -> setNDefVars(como + comc + comb + como1 + comc1);

  // see "hooking your solver to AMPL", by David M. Gay, tables 3, 4, and 5

  // nonlinear in both objectives and constraints
  if (nlvb >= 0) {
    for (int i = 0; i < nlvb - nlvbi; i++) problem -> addVariable (false, problem -> domain ());
    for (int i = 0; i < nlvbi;        i++) problem -> addVariable (true,  problem -> domain ());
  }

  // nonlinear in either objectives or constraints
  if (nlvo > nlvc) {
    for (int i = 0; i < nlvc - (nlvb + nlvci); i++) problem -> addVariable (false, problem -> domain ());
    for (int i = 0; i < nlvci;                 i++) problem -> addVariable (true,  problem -> domain ());
    for (int i = 0; i < nlvo - (nlvc + nlvoi); i++) problem -> addVariable (false, problem -> domain ());
    for (int i = 0; i < nlvoi;                 i++) problem -> addVariable (true,  problem -> domain ());
  } else {
    for (int i = 0; i < nlvo - (nlvb + nlvoi); i++) problem -> addVariable (false, problem -> domain ());
    for (int i = 0; i < nlvoi;                 i++) problem -> addVariable (true,  problem -> domain ());
    for (int i = 0; i < nlvc - (nlvo + nlvci); i++) problem -> addVariable (false, problem -> domain ());
    for (int i = 0; i < nlvci;                 i++) problem -> addVariable (true,  problem -> domain ());
  }

  for (int i = 0; i < nwv; i++)                                  problem -> addVariable(false, problem -> domain ());//arc
  for (int i = n_var - (CoinMax (nlvc,nlvo) +niv+nbv+nwv); i--;) problem -> addVariable(false, problem -> domain ());//other
  for (int i = 0; i < nbv; i++)                                  problem -> addVariable(true,  problem -> domain ());//binary
  for (int i = 0; i < niv; i++)                                  problem -> addVariable(true,  problem -> domain ());//int.

  // add space for common expressions
  for (int i = problem->nDefVars(); i--;)  problem -> addVariable(false, problem -> domain ());

  // common expressions (or defined variables) ///////////////////////////////////////

#ifdef DEBUG
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

  try {
  for (int i = 0; i < como + comc + comb; i++) {

    struct cexp *common = ((const ASL_fg *) asl) -> I.cexps_ + i;
    expression *nle = nl2e (common -> e);

#ifdef DEBUG
    printf ("cexp  %d [%d]: ", i, problem -> nVars()); nle -> print ();  printf (" ||| ");
#endif

    int nlin = common -> nlin;  // Number of linear terms

    if (nlin > 0) {

      int       *indexL = new int       [nlin+1];
      CouNumber *coeff  = new CouNumber [nlin];

      linpart *L = common -> L;

      for (int j = 0; j < nlin; j++) {
	//vp = (expr_v *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
	//Printf( " %-g x[%-d]", L->fac, (int)(vp - VAR_E) );	
	coeff [j] = L [j]. fac;
	indexL [j] = ((uintptr_t) (L [j].v.rp) - (uintptr_t) VAR_E) / sizeof (expr_v);
#ifdef DEBUG
	Printf( " %+g x_%-3d", L [j]. fac, 
		(expr_v *) (L [j].v.rp) - VAR_E //((const ASL_fg *) asl) -> I.cexps_
		//L [j]. v.i
		);
#endif
      }

      indexL [nlin] = -1;

      expression **al = new expression * [1];
      *al = nle;

      std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      problem -> indcoe2vector (indexL, coeff, lcoeff);

      expression *eg = exprGroup::genExprGroup (0, lcoeff, al, 1);
      problem -> commonExprs (). push_back (eg);
    } 
    else problem -> commonExprs () . push_back (nle);
#ifdef DEBUG
    printf ("\n");
#endif
  }

  for (int i = 0; i < como1 + comc1; i++) {

    struct cexp1 *common = ((const ASL_fg *) asl) -> I.cexps1_ + i;
    expression *nle = nl2e (common -> e);

#ifdef DEBUG
    printf ("cexp1 %d [%d]: ", i, variables_ . size ()); nle -> print ();  printf (" ||| ");
#endif

    int nlin = common -> nlin;  // Number of linear terms

    if (nlin > 0) {

      int       *indexL = new int       [nlin+1];
      CouNumber *coeff  = new CouNumber [nlin];

      linpart *L = common -> L;

      for (int j = 0; j < nlin; j++) {
	//vp = (expr_v *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
	coeff  [j] = L [j]. fac;
	indexL [j] = ((uintptr_t) (L [j].v.rp) - (uintptr_t) VAR_E) / sizeof (expr_v);
#ifdef DEBUG
	Printf( " %+g x_%-3d", L [j]. fac, 
		(expr_v *) (L [j].v.rp) - VAR_E //((const ASL_fg *) asl) -> I.cexps_
		//L [j]. v.i
		);
#endif
      }

      indexL [nlin] = -1;

      expression **al = new expression * [1];
      *al = nle;

      std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      problem -> indcoe2vector (indexL, coeff, lcoeff);

      expression *eg = exprGroup::genExprGroup (0, lcoeff, al, 1);
      problem -> commonExprs () . push_back (eg);
    } 
    else problem -> commonExprs () . push_back (nle);
#ifdef DEBUG
    printf ("\n");
#endif
    //    addAuxiliary (nl2e (((const ASL_fg *) asl) -> I.cexps1_ [i] . e));
  }

  // objective functions /////////////////////////////////////////////////////////////

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
      *nl = nl2e (OBJ_DE [i] . e);

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
      problem -> indcoe2vector (indexL, coeff, lcoeff);

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

    expression *subst = Simplified (body);//  -> simplify ();

    // if (subst) {
    //   delete body; // VALGRIND
    //   body = subst;
    // }

    // ThirdParty/ASL/solvers/asl.h, line 336: 0 is minimization, 1 is maximization
    problem -> addObjective (body, (OBJ_sense [i] == 0) ? "min" : "max");
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
    if (lb > negInfinity)
      if (ub < Infinity) sign = COUENNE_RNG;
      else               sign = COUENNE_GE;
    else                 sign = COUENNE_LE;

    // this is an equality constraint  
    if (fabs (lb - ub) < COUENNE_EPS)
      sign = COUENNE_EQ;

    expression *body;

    expression **nll = new expression * [1];
    *nll = nl2e (CON_DE [i] . e);

    if (indexL [i] && (*(indexL [i]) >= 0)) {

      int code = (*nll) -> code ();

      std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
      problem -> indcoe2vector (indexL [i], coeff [i], lcoeff);

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

    expression *subst = Simplified (body);

    // -> simplify ();
    // if (subst) {
    //   delete body; // VALGRIND
    //   body = subst;
    // }

    // add them (and set lower-upper bound)
    switch (sign) {

    case COUENNE_EQ:  problem -> addEQConstraint  (body, new exprConst (ub)); break;
    case COUENNE_LE:  problem -> addLEConstraint  (body, new exprConst (ub)); break;
    case COUENNE_GE:  problem -> addGEConstraint  (body, new exprConst (lb)); break;
    case COUENNE_RNG: problem -> addRNGConstraint (body, new exprConst (lb), 
					           new exprConst (ub)); break;
    default: jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_INITIALIZATION, "Error: could not recognize constraint\n"); return false;
    }

    delete [] indexL [i];
    delete [] coeff  [i];
  }

  delete [] indexL;
  delete [] coeff;
  delete [] nterms;
  
  } catch (...) {
  	return false;
  }

  // create room for problem's variables and bounds
  CouNumber 
    *x  = (CouNumber *) malloc ((n_var + problem -> nDefVars() ) * sizeof (CouNumber)),
    *lb = (CouNumber *) malloc ((n_var + problem -> nDefVars() ) * sizeof (CouNumber)),
    *ub = (CouNumber *) malloc ((n_var + problem -> nDefVars() ) * sizeof (CouNumber));

  for (int i = n_var + problem -> nDefVars(); i--;) {
    x  [i] =  0.;
    lb [i] = -COUENNE_INFINITY;
    ub [i] =  COUENNE_INFINITY;
  }

  problem -> domain () -> push (n_var + problem -> nDefVars(), x, lb, ub);
  free (x); free (lb); free (ub);

  // suggested:
  // problem -> domain () -> push (n_var + problem -> nDefVars(), x, lb, ub, false);
  // //free (x); free (lb); free (ub);
  // saves three allocations (default last parameter is true, which copies x,l,b)

  // lower and upper bounds ///////////////////////////////////////////////////////////////

  if (LUv) {

    real *Uvx_copy = Uvx;

    if (!Uvx_copy)
      for (int i=0; i<n_var; i++) {

	int j = 2*i;

        problem -> Lb (i) = LUv[j]   <= -COUENNE_INFINITY ? -COUENNE_INFINITY : LUv[j]  ;
        problem -> Ub (i) = LUv[j+1] >=  COUENNE_INFINITY ?  COUENNE_INFINITY : LUv[j+1];
      }
    else
      for (int i=n_var; i--;) {
	problem -> Lb (i) = LUv [i]      <= -COUENNE_INFINITY ? -COUENNE_INFINITY : LUv[i];
	problem -> Ub (i) = Uvx_copy [i] >=  COUENNE_INFINITY ?  COUENNE_INFINITY : Uvx_copy[i];
      }

  } else
    for (int i=n_var; i--;) {
    	problem -> Lb (i) = - COUENNE_INFINITY;
    	problem -> Ub (i) =   COUENNE_INFINITY;
    }

  // initial x ////////////////////////////////////////////////////////////////////

  for (int i=n_var; i--;) 

    if (X0 && havex0 [i]) problem -> X (i) = X0 [i]; 

    else {

      CouNumber x, l = problem -> Lb (i), u = problem -> Ub (i);

      if      (l < - COUENNE_INFINITY)
	if    (u >   COUENNE_INFINITY)  x = 0.;
	else                            x = u;
      else if (u >   COUENNE_INFINITY)  x = l;
      else                              x = 0.5 * (l+u);

      problem -> X (i) = x;
    }

  for (int i=n_var; i < problem -> nDefVars() ; i++) {  //FIXME: shouldn't this loop go until n_var + problem -> nDefVars() ?

  	problem -> X  (i) =  0.;
  	problem -> Lb (i) = -COUENNE_INFINITY;
  	problem -> Ub (i) =  COUENNE_INFINITY;
  }

  return true;
}


// warning for non-implemented functions -- return 0 constant expression
//expression *notimpl (const std::string &fname) {
//static void notimpl (const std::string &fname) {
//  std::cerr << "*** Error: " << fname << " not implemented" << std::endl;
//  exit (-1);
//}

// converts an AMPL expression (sub)tree into an expression* (sub)tree
expression *CouenneAmplInterface::nl2e(expr *e) {

  switch (getOperator (e -> op)) {

  case OPPLUS:  return new exprSum (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMINUS: return new exprSub (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMULT:  return new exprMul (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPDIV:   return new exprDiv (nl2e (e -> L.e), nl2e (e -> R.e));
    //case OPREM:   notimpl ("remainder");
  case OPPOW:   return new exprPow (nl2e (e -> L.e), nl2e (e -> R.e));
    //case OPLESS:  notimpl ("less");
    //case MINLIST: notimpl ("min");
    //case MAXLIST: notimpl ("max");
    //case FLOOR:   notimpl ("floor");
    //case CEIL:    notimpl ("ceil");
  case ABS:     return new exprAbs (nl2e (e -> L.e));
  case OPUMINUS:return new exprOpp (nl2e (e -> L.e));
    //          return new exprOpp (nl2e (e -> L.e -> L.e));
    //case OPIFnl:  { notimpl ("ifnl");

    // see ASL/solvers/rops.c, IfNL
    //}

  case OP_tanh: return new exprDiv 
      (new exprSub (new exprExp (nl2e (e -> L.e)),
		    new exprExp (new exprOpp (nl2e (e->L.e)))),
       new exprSum (new exprExp (nl2e (e -> L.e)),
		    new exprExp (new exprOpp (nl2e (e->L.e)))));

  case OP_tan: {
    expression *arg;
    arg = nl2e (e -> L.e);
    return new exprDiv (new exprSin (arg), new exprCos (new exprClone (arg)));
  }
  case OP_sqrt:    return new exprPow (nl2e (e -> L.e), new exprConst (0.5));
  case OP_sinh:    return new exprMul (new exprConst (0.5),
				       new exprSub (new exprExp (nl2e (e -> L.e)),
						    new exprExp (new exprOpp (nl2e (e->L.e)))));
  case OP_sin:     return new exprSin (nl2e (e -> L.e));
  case OP_log10:   return new exprMul (new exprConst (1.0 / log (10.0)), 
				       new exprLog (nl2e (e -> L.e)));
  case OP_log:     return new exprLog (nl2e (e -> L.e));
  case OP_exp:     return new exprExp (nl2e (e -> L.e));
  case OP_cosh:    return new exprMul (new exprConst (0.5),
				       new exprSum (new exprExp (nl2e (e -> L.e)),
						    new exprExp (new exprOpp (nl2e (e->L.e)))));

  case OP_cos:   return new exprCos (nl2e (e -> L.e));
    //case OP_atanh: notimpl ("atanh");
    //case OP_atan2: notimpl ("atan2");
    //case OP_atan:  notimpl ("atan");
    //case OP_asinh: notimpl ("asinh");
    //case OP_asin:  notimpl ("asin");
    //case OP_acosh: notimpl ("acosh");
    //case OP_acos:  notimpl ("acos");

  case OPSUMLIST: {
    int i=0;
    expression **al = new expression * [(e->R.ep - e->L.ep)];
    for (expr **ep = e->L.ep; ep < e->R.ep; ep++)
      al [i++] = nl2e (*ep);
    return new exprSum (al, i);
  }
    //case OPintDIV: notimpl ("intdiv");
    //case OPprecision: notimpl ("precision");
    //case OPround:  notimpl ("round");
    //case OPtrunc:  notimpl ("trunc");

  case OP1POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (((expr_n *)e->R.e)->v));
  case OP2POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (2.));
  case OPCPOW: return new exprPow (new exprConst (((expr_n *)e->L.e)->v),
				   nl2e (e -> R.e));
    //case OPFUNCALL: notimpl ("function call");
  case OPNUM:     return new exprConst (((expr_n *)e)->v);
    //case OPPLTERM:  notimpl ("plterm");
    //case OPIFSYM:   notimpl ("ifsym");
    //case OPHOL:     notimpl ("hol");
  case OPVARVAL:  {

    int j = ((expr_v *) e) -> a;

    if (j >= problem -> nOrigVars()) // common expression
      // use base pointer otherwise the .a field returns an awkward, out-of-bound index
      j = ((expr_v *) e) - ((const ASL_fg *) asl) -> I.var_e_; 

    if (j >= problem -> nOrigVars() + problem -> nDefVars()) {
      jnlst -> Printf (Ipopt::J_ERROR, Ipopt::J_INITIALIZATION, "Error: unknown variable x_%d\n", j);
      throw -1;
    }

    return new exprClone (problem -> Variables() [j]);
  }

  default:
    jnlst -> Printf (Ipopt::J_ERROR, Ipopt::J_INITIALIZATION, "ERROR: unknown operator (address %p), aborting.\n", Intcast (e -> op));
    throw -2;
  }

  return new exprConst (0.);
}
