/* $Id$
 *
 * Name:    TwoImpliedGenCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate cuts using two inequalities from the LP relaxation
 * 
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Common Public License (CPL)
 */

#include <stdlib.h>

#include "CglCutGenerator.hpp"
#include "CouenneTwoImplied.hpp"
#include "CoinPackedMatrix.hpp"
#include "CouennePrecisions.hpp"

using namespace Couenne;

enum signum {NEG, POS, DN, UP};


// type for alphas ////////////////////////////////////////////////

typedef struct {

  double alpha;
  int index;
  enum signum sign; // defines behavior of convex combination within
                    // [0,1]: always negative, always positive, + ->
                    // -, and - -> +
} threshold;


// compare two threshold values ///////////////////////////////////

int compthres (const void *t1, 
	       const void *t2) {

  double
    a1 = (*(threshold **) t1) -> alpha,
    a2 = (*(threshold **) t2) -> alpha;

  return (a1 < a2 ? -1 : 
	  a1 > a2 ?  1 : 0);
}


//
// Combines two inequalities whose indices, bounds, and coefficients
// are defined by the arguments ind, l and u, and a, respectively.
//

int combine (OsiCuts &cs, 
	     int n1, int n2, 
	     const int *ind1, // indices
	     const int *ind2, 
	     double *sa1, // coeff (sparse array) 
	     double *sa2, // DO NOT INVERT (already done at caller)
	     const double *a1,  // coeff
	     const double *a2, 
	     double *clb, // variable bounds
	     double *cub,
	     double l1, // constraint bounds
	     double l2,  
	     double u1, 
	     double u2, 
	     int sign) { // invert second constraint? -1: yes, +1: no

  printf ("here are the two ineqs:\n");
  printf ("%g <=", l1);
  for (int i=0; i<n1; i++) printf (" %+g x%d", a1 [i], ind1 [i]);
  printf ("<= %g\n%g <=", u1, l2);
  for (int i=0; i<n2; i++) printf (" %+g x%d", a2 [i], ind2 [i]);
  printf ("<= %g\n", u2);

  double signA;

  if (sign < 0) {

    double tmp = u2;
    u2 = -l2;
    l2 = -tmp;

    signA = -1.;

  } else signA = 1.;

  threshold
     *alphas   = new threshold  [n1 + n2], // contains all alphas (there might be at most n1+n2)
    **inalphas = new threshold* [n1 + n2], // points to those in [0,1[ (will be sorted)
     *curalpha = alphas;

  int 
    i1    = 0, 
    i2    = 0, 
    cnt   = 0,
    incnt = 0;

  while (i1 < n1 ||
	 i2 < n2) {

    if        (i1 < n1 && (i2 == n2 || ind1 [i1] < ind2 [i2])) {

      curalpha -> alpha = 1.;
      curalpha -> sign  = a1 [i1] < 0. ? NEG : POS;
      curalpha -> index = ind1 [i1];

      cnt++;
      i1++;

    } else if (i2 < n2 && (i1 == n1 || ind2 [i2] < ind1 [i1])) {

      curalpha -> alpha = 0.;
      curalpha -> sign  = signA * a2 [i2] < 0. ? NEG : POS;
      curalpha -> index = ind2 [i2];

      cnt++;
      i2++;

    } else {

      // the two indices are equal, < n1, n2 respectively, we may have
      // an alpha in ]0,1[

      double a2i2 = signA * a2 [i2];

      curalpha -> alpha = 
	(a1 [i1] == a2i2) ? 
	-1. : 
	- a2i2 / (a1 [i1] - a2i2);

      if (curalpha -> alpha <= 0. || 
	  curalpha -> alpha >= 1.)

	curalpha -> sign = a1 [i1] > 0. ? POS : NEG;

      else 

	curalpha -> sign = 
	  a1 [i1] == a2i2 ? 
	  (a1 [i1] < 0. ? NEG : POS) :
	  (a1 [i1] < 0. ? UP  : DN);

      curalpha -> index = ind1 [i1];

      if (curalpha -> alpha > 0. &&
	  curalpha -> alpha < 1.)
	inalphas [incnt++] = curalpha;

      cnt++;
      i1++;
      i2++;
    }

    printf ("(%d,%g,%s) ", 
    	    curalpha -> index, 
    	    curalpha -> alpha, 
    	    curalpha -> sign == NEG ? "NEG" :
    	    curalpha -> sign == POS ? "POS" :
    	    curalpha -> sign == DN  ? "DN"  : "UP");

    curalpha++;
  }

  printf ("\n");

  int ntightened = 0;

  if (incnt) {

    if (incnt > 1)
      qsort (inalphas, incnt, sizeof (threshold *), compthres);

    printf ("sorted: ");
    for (int i = 0; i < incnt; i++)
      printf ("(%d,%g,%s) ", 
    	      inalphas [i] -> index, 
    	      inalphas [i] -> alpha, 
    	      inalphas [i] -> sign == NEG ? "NEG" :
    	      inalphas [i] -> sign == POS ? "POS" :
    	      inalphas [i] -> sign == DN  ? "DN"  : "UP");    
    printf ("\n");

    // if none of them has an alpha in ]0,1[, nothing needs to be
    // done.

    // Otherwise, this is the actual procedure. All we need now is
    // alphas, inalphas, rlb/rub, clb/cub and a lot of debugging

    // Compute 
    //
    // L+ = sum {i in N+} c_i xL_i + sum {i in N-} c_i xU_i 
    // L- = sum {i in N+} c_i xU_i + sum {i in N-} c_i xL_i 
    //
    // with all c_i computed at alpha=0, therefore we consider all
    // coefficients of the second constraint

    double 
      minSum1 = 0., minSum2 = 0.,
      maxSum1 = 0., maxSum2 = 0.;

    int 
      mInfs = 0, // number of -inf summed to minSum 1 and 2 (both)
      pInfs = 0; //           +inf           maxSum

    // compute for a_2 ////////////////////////////////////////////////////////

    for (int i=0; i<n2; i++) {

      int index2 = ind2 [i];

      double 
	a2i  = a2  [i] * signA,
	clb2 = clb [index2],
	cub2 = cub [index2];

      if (a2i < 0.) {

	if (cub2 >   COUENNE_INFINITY) mInfs ++; else minSum2 += a2i * cub2;
	if (clb2 < - COUENNE_INFINITY) pInfs ++; else maxSum2 += a2i * clb2;

      } else {

	if (clb2 < - COUENNE_INFINITY) mInfs ++; else minSum2 += a2i * clb2;
	if (cub2 >   COUENNE_INFINITY) pInfs ++; else maxSum2 += a2i * cub2;
      }
    }

    // copy into 1's data. We are at alpha = 0, where everything looks like 2

    maxSum1 = maxSum2;
    minSum1 = minSum2;

    //
    // scan all alphas in ]0,1[ ///////////////////////////////////////////////
    //

    // for all variables x[i], look for tighter lower/upper bounds on x[i]

    // denominator of all variables --- SEPARATED PER CONSTRAINT, and sparse
    double *den = new double [n1 + n2];
    CoinFillN (den, n1, 0.);
    for (int i=0; i<n2; i++) 
      den [n1 + i] = signA * a2 [i];

    std::vector <std::pair <int, double> > 
      newLB, // pairs (index, value) of new bounds
      newUB;

    for (int i = 0; i < incnt;) {

      printf ("  looking at %d: (%d,%g,%s); ", i,
    	      inalphas [i] -> index, 
    	      inalphas [i] -> alpha, 
    	      inalphas [i] -> sign == NEG ? "NEG" :
    	      inalphas [i] -> sign == POS ? "POS" :
    	      inalphas [i] -> sign == DN  ? "DN"  : "UP");

      for (int j=0; j < n1; j++) printf ("x%d [%g,%g] ", ind1 [j], clb [ind1 [j]], cub [ind1 [j]]);
      for (int j=0; j < n2; j++) printf ("x%d [%g,%g] ", ind2 [j], clb [ind2 [j]], cub [ind2 [j]]);

      printf ("\n");

      // look at each inalphas, those are the only breakpoints where
      // bounds can improve. For alpha in {0,1} there is already a
      // procedure (single constraint implied bound).

      // check for improved bounds at all variables

      double alpha = inalphas [i] -> alpha;

      for (int j=0; j<n1+n2; j++) {

	if (j >= n1 && sa1 [ind2 [j - n1]] != 0.) 
	  continue; // scan all nonzeros, those on one side and those
		    // on both sides. This has already been scanned.

	int index = ((j < n1) ? ind1 [j] : ind2 [j - n1]);

	double
	  den_j = den [j] + alpha * (sa1 [index] - sa2 [index]),

	  newL = - COUENNE_INFINITY,
	  newU =   COUENNE_INFINITY,

	  clbi = clb [index],
	  cubi = cub [index],

	  sa1i = sa1 [index],
	  sa2i = sa2 [index],

	  ci = alpha * sa1i + (1 - alpha) * sa2i, // combination of coefficients

	  subMin1 = 0., subMin2 = 0.,
	  subMax1 = 0., subMax2 = 0.;

	int 
	  tickMin = 0, // one if subtracting this term removes one infinity from minSum
 	  tickMax = 0; //                                                        max   

	// m[in,ax]Sum[1,2]A must be updated by subtracting the term
	// that is now at the denominator (no need to multiply by
	// signA, sa2 is already adjusted)

	if (clbi < - COUENNE_INFINITY) {

	  // we are deleting an infinite bound from the numerator, and
	  // whether it goes towards max- or min- Sum? depends on each
	  // coefficient

	  if (ci > 0.) tickMin++; 
	  else         tickMax++;

	} else {

	  if (ci > 0.) {

	    subMin1 = sa1i * clbi; 
	    subMin2 = sa2i * clbi; 

	  } else {

	    subMax1 = sa1i * clbi;
	    subMax2 = sa2i * clbi;
	  }
	}

	if (cubi >   COUENNE_INFINITY) {

	  // we are deleting an infinite bound from the numerator, and
	  // whether it goes towards max- or min- Sum?A it depends on each coefficient

	  if (ci > 0.) tickMax++; 
	  else         tickMin++;

	} else {

	  if (ci > 0.) {

	    subMax1 = sa1i * cubi; 
	    subMax2 = sa2i * cubi; 

	  } else {

	    subMin1 = sa1i * cubi;
	    subMin2 = sa2i * cubi;
	  }
	}

	if (fabs (den_j) < 1.e-50) { // must check that the numerator behaves

	} else if (den_j > 0.) {

	  if ((l1 > -COUENNE_INFINITY) && 
	      (l2 > -COUENNE_INFINITY) &&
	      (pInfs == tickMax))
	    newL = ((l1 - l2 - (maxSum1 - maxSum2) + (subMax1 - subMax2)) * alpha + l2 - maxSum2 + subMax2) / den_j;

	  if ((u1 <  COUENNE_INFINITY) && 
	      (u2 <  COUENNE_INFINITY) &&
	      (mInfs == tickMin))
	    newU = ((u1 - u2 - (minSum1 - minSum2) + (subMin1 - subMin2)) * alpha + u2 - minSum2 + subMin2) / den_j;

	} else {

	  if ((u1 <  COUENNE_INFINITY) && 
	      (u2 <  COUENNE_INFINITY) &&
	      (mInfs == tickMin))
	    newL = ((u1 - u2 - (minSum1 - minSum2) + (subMin1 - subMin2)) * alpha + u2 - minSum2 + subMin2) / den_j;

	  if ((l1 > -COUENNE_INFINITY) && 
	      (l2 > -COUENNE_INFINITY) &&
	      (pInfs == tickMax))
	    newU = ((l1 - l2 - (maxSum1 - maxSum2) + (subMax1 - subMax2)) * alpha + l2 - maxSum2 + subMax2) / den_j;
	}

	if (newL > clbi) {

	  ntightened++;
	  newLB. push_back (std::pair <int, double> (index, newL));
	  printf ("  new bound: x%d >= %g [%g]\n", index, newL, clb [index]);
	}

	if (newU < cubi) {

	  ntightened++;
	  newUB. push_back (std::pair <int, double> (index, newU));
	  printf ("  new bound: x%d <= %g [%g]\n", index, newU, cub [index]);
	}
      }

      // enter next interval, update dynamic max/min sums

      do {

	int 
	  indalpha  = inalphas [i] -> index,
	  signalpha = inalphas [i] -> sign;

	if (clb [indalpha] < - COUENNE_INFINITY) {

	  if (signalpha == DN) {pInfs++; mInfs--;} 
	  else                 {pInfs--; mInfs++;} // sign is UP

	} else {

	  if (signalpha == DN) { // means sa1 [indalpha] > 0 while sa2 [indalpha] < 0

	    minSum1 -= sa1 [indalpha] * clb [indalpha];
	    minSum2 -= sa2 [indalpha] * clb [indalpha];
	    maxSum1 += sa1 [indalpha] * clb [indalpha];
	    maxSum2 += sa2 [indalpha] * clb [indalpha];

	  } else { // sign must be DN

	    maxSum1 -= sa1 [indalpha] * clb [indalpha];
	    maxSum2 -= sa2 [indalpha] * clb [indalpha];
	    minSum1 += sa1 [indalpha] * clb [indalpha];
	    minSum2 += sa2 [indalpha] * clb [indalpha];
	  }
	}

	// same for upper bound

	if (cub [indalpha] >   COUENNE_INFINITY) {

	  if (signalpha == DN) {pInfs--; mInfs++;} 
	  else                 {pInfs++; mInfs--;} // sign is UP

	} else {

	  if (signalpha == DN) { // means sa1 [indalpha] > 0 while sa2 [indalpha] < 0

	    maxSum1 -= sa1 [indalpha] * cub [indalpha];
	    maxSum2 -= sa2 [indalpha] * cub [indalpha];
	    minSum1 += sa1 [indalpha] * cub [indalpha];
	    minSum2 += sa2 [indalpha] * cub [indalpha];

	  } else { // sign must be DN

	    minSum1 -= sa1 [indalpha] * cub [indalpha];
	    minSum2 -= sa2 [indalpha] * cub [indalpha];
	    maxSum1 += sa1 [indalpha] * cub [indalpha];
	    maxSum2 += sa2 [indalpha] * cub [indalpha];
	  }
	}

      } while (++i < incnt && alpha == inalphas [i] -> alpha);
    }

    for (std::vector <std::pair <int, double> >::iterator i = newLB. begin (); i != newLB. end (); ++i) if (i -> second > clb [i -> first]) clb [i -> first] = i -> second;
    for (std::vector <std::pair <int, double> >::iterator i = newUB. begin (); i != newUB. end (); ++i) if (i -> second < cub [i -> first]) cub [i -> first] = i -> second;
  }

  delete [] inalphas;
  delete [] alphas;

  return ntightened;
}
