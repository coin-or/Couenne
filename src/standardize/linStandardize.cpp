/*
 *
 * Name:    linStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize sum expressions (expr{Sum,Sub,Quad,Group,Opp})
 *
 * (C) Carnegie-Mellon University, 2007-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>

#include "CouenneExprQuad.hpp"

#include "CouenneProblem.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprSum.hpp"
#include "CouenneExprSub.hpp"
#include "CouenneExprOpp.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprPow.hpp"
#include "CouenneExprGroup.hpp"
#include "CouenneLQelems.hpp"

using namespace Couenne;

/// standardization of linear exprOp's
exprAux *CouenneProblem::linStandardize (bool addAux,
					 CouNumber c0,
					 LinMap  &lmap,
					 QuadMap &qmap) {

  // check if quadratic forms are dense enough ///////////////////////////////////////////

  analyzeSparsity (c0, lmap, qmap);

  ////////////////////////////////////////////////////////////////////////////////////////

  int  nq = qmap.Map().size (),     /// data for exprQuad
      *qi = new int [nq+1],
      *qj = new int [nq+1];

  CouNumber *qc = new CouNumber [nq];

  int
     nl = lmap.Map().size(),      /// data for exprGroup
    *li = new int [nl+1];

  CouNumber *lc = new CouNumber [nl];

  // terminate arrays with negative index
  qi [nq] = li [nl] = -1;

  std::map <int, CouNumber>::iterator lit = lmap.Map().begin ();

  // fill in arrays for linear part
  for (int i=0; i<nl; i++, lit++) {

    li [i] = lit -> first;
    lc [i] = lit -> second;
  }

  std::map <std::pair <int, int>, CouNumber>::iterator qit = qmap.Map().begin ();

  // fill in arrays for quadratic part
  for (int i=0; i < nq; i++, qit++) {
    qi [i] = qit -> first. first;
    qj [i] = qit -> first. second;
    qc [i] = qit -> second;
  }

  nl = lmap.Map().size ();
  nq = qmap.Map().size ();

  // particular cases ///////////////////////////////////////////////////////////

  expression *ret;

  // a constant
  if ((nq==0) && (nl==0))

    ret = new exprConst (c0); // a constant auxiliary?

  else if ((nq==0) && (fabs (c0) < COUENNE_EPS) && (nl==1)) { // a linear monomial, cx

    if      (fabs (*lc - 1.) < COUENNE_EPS) ret = new exprClone (Var (*li));
    else if (fabs (*lc + 1.) < COUENNE_EPS) ret = new exprOpp (new exprClone (Var (*li)));
    else                                    ret = new exprMul (new exprConst (*lc),
				 			       new exprClone (Var (*li)));

  } else if ((nl==0) && (fabs (c0) < COUENNE_EPS) && (nq==1)) {

    // a bilinear/quadratic monomial, cx^2 or cxy

    expression *quad;

    if (*qi == *qj) quad = new exprPow (new exprClone (Var (*qi)), new exprConst (2.));
    else            quad = new exprMul (new exprClone (Var (*qi)),
					new exprClone (Var (*qj)));

    if (fabs (*qc - 1) < COUENNE_EPS)
      ret    = quad;
    else {
      quad = addAuxiliary (quad);
      ret  = new exprMul (new exprConst (*qc), new exprClone (quad));
    }

  } else {

    // general case ///////////////////////////////////////////////////////////////

    std::vector <std::pair <exprVar *, CouNumber> > lcoeff;
    indcoe2vector (li, lc, lcoeff);

    std::vector <quadElem> qcoeff;
    indcoe2vector (qi, qj, qc, qcoeff);

    if (!nq) ret = new exprGroup (c0, lcoeff); // create exprGroup
    else     ret = new exprQuad  (c0, lcoeff, qcoeff); // create exprQuad
  }

  delete [] li;
  delete [] lc;
  delete [] qi;
  delete [] qj;
  delete [] qc;

  if (jnlst_ -> ProduceOutput (Ipopt::J_ALL, J_REFORMULATE)) {
    printf ("\nlinstand (addaux = %d) ==> ", addAux);
    ret -> print (); printf ("\n");
    //printf (": "); fflush (stdout);
    //ret -> Image () -> print ();
  }

  //if (ret -> Type () == AUX)
  //return ret;

  return (addAux ? addAuxiliary (ret) : new exprAux (ret, &domain_));
}
