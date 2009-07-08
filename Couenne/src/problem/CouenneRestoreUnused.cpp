/* $Id$
 *
 * Name:    CouenneRestoreUnused.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Restore a consistent value of duplicate variables 
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"


/// Some originals may be unused due to their zero multiplicity
/// (that happens when they are duplicates). This procedure creates
/// a structure for quickly checking and restoring their value after
/// solving.
void CouenneProblem::createUnusedOriginals () {

  if (nUnusedOriginals_ < 0) { // no structure yet, identify and store

    nUnusedOriginals_ = 0;

    int 
      nOrig = nOrigVars (),
      nvars = nVars     ();

    unusedOriginalsIndices_ = (int *) malloc (nvars * sizeof (int)); // will trim it later

    for (int i=0; i<nvars; i++) {

      int indVar = numbering_ [i];

      if ((indVar < nOrig) && 
	  (variables_ [indVar] -> Multiplicity () <= 0)) // found neglected variable!
	unusedOriginalsIndices_ [nUnusedOriginals_++] = indVar;
    }

    if (nUnusedOriginals_)
      unusedOriginalsIndices_ = (int *) realloc (unusedOriginalsIndices_, 
						 nUnusedOriginals_ * sizeof (int));
    else {
      free (unusedOriginalsIndices_);
      unusedOriginalsIndices_ = NULL;
    }
  }
}


/// Some originals may be unused due to their zero multiplicity (that
/// happens when they are duplicates). This procedure restores their
/// value after solving
void CouenneProblem::restoreUnusedOriginals (CouNumber *x) const {

  if (nUnusedOriginals_ <= 0) return;

  if (x)
    domain_.push (nVars(), x, NULL, NULL, false); // no need for another copy

  for (int i=0; i<nUnusedOriginals_; i++) {

    int indVar = unusedOriginalsIndices_ [i];
    expression *img = variables_ [indVar] -> Image ();

    if (img) {

      CouNumber value = (*img) ();

      X  (indVar) = value;

      if (x)
	x  [indVar] = value;
    }
  }

  if (x)
    domain_. pop ();
}
