/* $Id$
 *
 * Name:    exprVar.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class for defining variables
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneExprAux.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprBound.hpp"
#include "CouenneDepGraph.hpp"

using namespace Couenne;

// Get lower and upper bound of a variable expression (if any)
void exprVar::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLowerBound (varIndex_, domain_); 
  ub = new exprUpperBound (varIndex_, domain_);
}


// Get value of lower and upper bound of a variable expression (if any)
void exprVar::getBounds (CouNumber &lb, CouNumber &ub) {

  lb = domain_ -> lb (varIndex_);
  ub = domain_ -> ub (varIndex_);
}


// generate convexification cut for constraint w = this
void exprVar::generateCuts (expression *w, //const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) {
  if (cg -> isFirst ())
    cg -> createCut (cs, 0., cg -> Problem () -> Var (w -> Index ()) -> sign (), w -> Index (), 1., varIndex_, -1);
}


/// implied bound processing. Expression w = x, upon change in lower
/// or upper bound of w, whose index is wind
bool exprVar::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  bool res = false;

  if (updateBound (-1, l + varIndex_, sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind])) {res = true; chg [varIndex_].setLower(t_chg_bounds::CHANGED);}
  if (updateBound (+1, u + varIndex_, sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind])) {res = true; chg [varIndex_].setUpper(t_chg_bounds::CHANGED);}

  return res;
}


/// update dependence set with index of this variable
void exprVar::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) 
{dep -> insert (g -> lookup (varIndex_));}


expression *exprVar::Lb () {return new exprLowerBound (varIndex_, domain_);}///< lower bound
expression *exprVar::Ub () {return new exprUpperBound (varIndex_, domain_);}///< upper bound


// dummy function -- the real one is in exprAux and returns either a
// CouenneObject or a CouenneComplObject
CouenneObject *exprVar::properObject (CouenneCutGenerator *c,
				      CouenneProblem *p, 
				      Bonmin::BabSetupBase *base, 
				      JnlstPtr jnlst_) {
  return new CouenneObject ();
}
