/* $Id: CouenneSOS.cpp 154 2009-06-16 18:52:53Z pbelotti $ */
/*
 * Name:    CouenneSOS.cpp
 * Author:  Pietro Belotti
 * Purpose: find SOS in problem 
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CouenneProblem.hpp"
#include "exprGroup.hpp"

#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareActual.hpp"
#include "CouenneSOSObject.hpp"

/// find SOS objects
int CouenneProblem::findSOS (OsiSolverInterface *solver,
			     OsiObject **objects) {

  // check auxiliaries defined as 
  // x_i binary. Disable it and add relative SOS to array "objects"

  int nSOS = 0;

  for (std::vector <exprVar *>::const_iterator v = variables_.begin ();
       v != variables_.end (); ++v) 

    if (((*v) -> Multiplicity () >    0) &&
	((*v) -> Type         () == AUX) &&
	((*v) -> Image () -> code () == COU_EXPRGROUP)) {

      expression *img = (*v) -> Image ();

      exprGroup *group = dynamic_cast <exprGroup *> (img -> isaCopy () ? 
						     img -> Copy () :
						     img);
      if (!group)
	continue;

      int wind = (*v) -> Index ();
      CouNumber cterm = group -> getc0 ();
      bool defVar = true;

      // now check if this is 
      //
      // 1) an auxiliary fixed to one  ==> it's a SOS if its image is x1+x2+...+xk
      // 2) an auxiliary fixed to zero ==> it's a SOS if its image is -x1-x2-...-xk+1

      if      (fabs (cterm - 1.) < COUENNE_EPS) defVar = false;
      else if (fabs (cterm)      > COUENNE_EPS) continue; // and defVar is true

      if (defVar &&
	  ((fabs (Lb (wind) - 1.) > COUENNE_EPS) ||
	   (fabs (Ub (wind) - 1.) > COUENNE_EPS)))
	continue;

      if (!defVar &&
	  ((fabs (Lb (wind)) > COUENNE_EPS)))
	continue;

      int lsz = group -> lcoeff (). size ();

      if (((lsz <= 2) &&  defVar) ||
	  ((lsz <= 1) && !defVar))
	continue;

      // there are two possibilities:
      //
      // 1) w is defined as w = 1 - \sum_{i=0}^n x_i  -- defvar = false
      //
      // 2) w is defined as \sum_{i=0}^n x_i and w \in [1,1] -- defvar = true

      bool
	intSOS = (*v) -> isInteger (),
	isSOS  = true,
	onlyOrigVars = true; // if SOS constraint only contains
			     // original variables, it has been
			     // spotted by Cbc already

      exprGroup::lincoeff &lcoe = group -> lcoeff ();

      for (exprGroup::lincoeff::iterator l = lcoe. begin (); 
	   l != lcoe. end (); ++l) {

	if ((fabs (l -> second - (defVar ? 1. : -1.)) > COUENNE_EPS) || // wrong coefficient?
	    (fabs (Lb (l -> first -> Index ()))       > COUENNE_EPS)) { // positive lower bound?

	  isSOS = false;
	  break;
	} else 
	  if (!(l -> first -> isInteger ()))
	    intSOS = false;

	if (l -> first -> Index () > nOrigVars_) // 
	  onlyOrigVars = false;
      }

      if (!isSOS || !intSOS)// || onlyOrigVars) 
	continue;

      printf ("----- found SOS: ");
      (*v) -> print (); printf (" := ");
      (*v) -> Image () -> print (); printf ("\n");

      // it is a SOS -- if intSOS==true, it's also integer

      int
	indStart = defVar ? 0 : 1,
	nelem    = indStart + lcoe. size (), 
	*indices = new int [nelem];

      if (!defVar)
	indices [0] = (*v) -> Index ();

      for (int i=indStart, j=0; i<nelem; i++)
	indices [i] = lcoe [j++]. first -> Index ();

      // TODO: if use Cbc, add CbcSOSBranchingObject

      /*CouenneSOSObject *newsos = new CouenneSOSObject (solver, nelem, indices, NULL, 1,
	jnlst_, true, true);*/

      OsiSOS *newsos = new OsiSOS (solver, nelem, indices, NULL, 1);

      objects [nSOS] = newsos;
      // as in BonBabSetupBase.cpp:675
      newsos -> setPriority (10);
      newsos -> setIntegerValued (intSOS);

      nSOS++;
    }

  return nSOS;
}
