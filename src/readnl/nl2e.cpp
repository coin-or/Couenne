/*
 *
 * Name:    nl2e.cpp
 * Author:  Pietro Belotti
 * Purpose: converts a nl expression into a Couenne expression
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneTypes.hpp"

#include "CouenneExprVar.hpp"
#include "CouenneExprAbs.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprDiv.hpp"
#include "CouenneExprInv.hpp"
#include "CouenneExprSin.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprClone.hpp"
#include "CouenneExprLog.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprCos.hpp"
#include "CouenneExprExp.hpp"

#include "asl.h"
#include "nlp.h"
#include "opcode.hd"

using namespace Couenne;

// get ASL op. code relative to function pointer passed as parameter
size_t getOperator (efunc *);


// warning for non-implemented functions -- return 0 constant expression
//expression *notimpl (const std::string &fname) {
void notimpl (const std::string &fname) {
  std::cerr << "*** Error: " << fname << " not implemented" << std::endl;
  exit (-1);
}


// converts an AMPL expression (sub)tree into an expression* (sub)tree
expression *CouenneProblem::nl2e (expr *e, const ASL *asl) {

  switch (getOperator (e -> op)) {

  case OPPLUS:  return new exprSum (nl2e (e -> L.e, asl), nl2e (e -> R.e, asl));
  case OPMINUS: return new exprSub (nl2e (e -> L.e, asl), nl2e (e -> R.e, asl));
  case OPMULT:  return new exprMul (nl2e (e -> L.e, asl), nl2e (e -> R.e, asl));
  case OPDIV:   return new exprDiv (nl2e (e -> L.e, asl), nl2e (e -> R.e, asl));
  case OPREM:   notimpl ("remainder");
  case OPPOW:   return new exprPow (nl2e (e -> L.e, asl), nl2e (e -> R.e, asl));
  case OPLESS:  notimpl ("less");
  case MINLIST: notimpl ("min");
  case MAXLIST: notimpl ("max");
  case FLOOR:   notimpl ("floor");
  case CEIL:    notimpl ("ceil");
  case ABS:     return new exprAbs (nl2e (e -> L.e, asl));
  case OPUMINUS:return new exprOpp (nl2e (e -> L.e, asl));
    //          return new exprOpp (nl2e (e -> L.e -> L.e, asl));
  case OPIFnl:  { notimpl ("ifnl");

    // see ASL/solvers/rops.c, IfNL
  }

  case OP_tanh: return new exprDiv
      (new exprSub (new exprExp (nl2e (e -> L.e, asl)),
		    new exprExp (new exprOpp (nl2e (e->L.e, asl)))),
       new exprSum (new exprExp (nl2e (e -> L.e, asl)),
		    new exprExp (new exprOpp (nl2e (e->L.e, asl)))));

  case OP_tan:
    return new exprDiv (new exprSin (nl2e (e -> L.e, asl)), new exprCos (new exprClone (nl2e (e -> L.e, asl))));
  case OP_sqrt:    return new exprPow (nl2e (e -> L.e, asl), new exprConst (0.5));
  case OP_sinh:    return new exprMul (new exprConst (0.5),
				       new exprSub (new exprExp (nl2e (e -> L.e, asl)),
						    new exprExp (new exprOpp (nl2e (e->L.e, asl)))));
  case OP_sin:     return new exprSin (nl2e (e -> L.e, asl));
  case OP_log10:   return new exprMul (new exprConst (1.0 / log (10.0)),
				       new exprLog (nl2e (e -> L.e, asl)));
  case OP_log:     return new exprLog (nl2e (e -> L.e, asl));
  case OP_exp:     return new exprExp (nl2e (e -> L.e, asl));
  case OP_cosh:    return new exprMul (new exprConst (0.5),
				       new exprSum (new exprExp (nl2e (e -> L.e, asl)),
						    new exprExp (new exprOpp (nl2e (e->L.e, asl)))));

  case OP_cos:   return new exprCos (nl2e (e -> L.e, asl));
  case OP_atanh: notimpl ("atanh");
  case OP_atan2: notimpl ("atan2");
  case OP_atan:  notimpl ("atan");
  case OP_asinh: notimpl ("asinh");
  case OP_asin:  notimpl ("asin");
  case OP_acosh: notimpl ("acosh");
  case OP_acos:  notimpl ("acos");

  case OPSUMLIST: {
    int i=0;
    expression **al = new expression * [(e->R.ep - e->L.ep)];
    for (expr **ep = e->L.ep; ep < e->R.ep; ep++)
      al [i++] = nl2e (*ep, asl);
    return new exprSum (al, i);
  }
  case OPintDIV: notimpl ("intdiv");
  case OPprecision: notimpl ("precision");
  case OPround:  notimpl ("round");
  case OPtrunc:  notimpl ("trunc");

  case OP1POW: return new exprPow (nl2e (e -> L.e, asl), 		   new exprConst (((expr_n *)e->R.e)->v));
  case OP2POW: return new exprPow (nl2e (e -> L.e, asl), 		   new exprConst (2.));
  case OPCPOW: return new exprPow (new exprConst (((expr_n *)e->L.e)->v),  nl2e (e -> R.e, asl));
  case OPFUNCALL: notimpl ("function call");
  case OPNUM:     return new exprConst (((expr_n *)e)->v);
  case OPPLTERM:  notimpl ("plterm");
  case OPIFSYM:   notimpl ("ifsym");
  case OPHOL:     notimpl ("hol");
  case OPVARVAL:  {

    int j = ((expr_v *) e) -> a;

    if (j >= nOrigVars_) // common expression
      // use base pointer otherwise the .a field returns an awkward, out-of-bound index
      // TODO: fix! In itointqor.nl should return v51=y44 but returns v52=y44
      //                                          v??=y39 but returns v79=y39
      j = ((expr_v *) e) - ((const ASL_fg *) asl) -> I.var_e_;

    if (j >= nOrigVars_ + ndefined_) {
      printf ("error: unknown variable x_%d\n", j);
      //return new exprClone (variables_ [0]);
      exit (-1);
    }

    return new exprClone (variables_ [j]);
  }

  default:
    printf ("Couenne error: unknown operator (address %p), aborting.\n", Intcast (e -> op));
    exit (-1);
    //return new exprConst (0);
  }

  return new exprConst (0.);
}
