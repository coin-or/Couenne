/* $Id$
 *
 * Name:    TwoImpliedCombine.cpp
 * Author:  Pietro Belotti
 * Purpose: Bound reduction using two inequalities from the LP relaxation
 * 
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdlib.h>

#include "CglCutGenerator.hpp"
#include "CouenneTwoImplied.hpp"
#include "CoinPackedMatrix.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprVar.hpp"

#define MIN_DENOM 1.e-10

namespace Couenne {

enum signum {NEG, POS, DN, UP};

  //#define DEBUG

// type for alphas ////////////////////////////////////////////////

typedef struct {

  double alpha;
  int indVar;
  enum signum sign; // defines behavior of convex combination within
                    // [0,1]: always negative, always positive, + ->
                    // -, and - -> +
} threshold;


// compare two threshold values ///////////////////////////////////

int compthres (register const void *t1, 
	       register const void *t2) {

  register double
    a1 = (*(threshold **) t1) -> alpha,
    a2 = (*(threshold **) t2) -> alpha;

  return ((a1 < a2) ? -1 : 
	  (a1 > a2) ?  1 : 0);
}

// structure to sort indices with respective values

struct indPosPair {

  int index;
  int position;
};

// compare two pairs (index,value) ///////////////////////////////////

int compPair (register const void *p1, 
	      register const void *p2) {

  register int
    i1 = ((struct indPosPair *) p1) -> index,
    i2 = ((struct indPosPair *) p2) -> index;

  return ((i1 < i2) ? -1 : 
	  (i1 > i2) ?  1 : 0);
}


//
// Combines two inequalities whose indices, bounds, and coefficients
// are defined by the arguments ind, l and u, and a, respectively.
//

int combine (CouenneProblem *p,
	     int n1, 
	     int n2, 
	     const int *ind1c, // indices
	     const int *ind2c, 
	     double *sa1, // coeff (sparse array) 
	     double *sa2, // DO NOT INVERT (already done at caller)
	     const double *a1c,  // coeff
	     const double *a2c, 
	     double *clb, // variable bounds
	     double *cub,
	     double l1, // constraint bounds
	     double l2,  
	     double u1, 
	     double u2,
	     bool *isInteger,
	     int sign) { // invert second constraint? -1: yes, +1: no

  // first, sort ind1/a1 and ind2/a2 w.r.t. indices. They may not be
  // sorted and this messes up the while loop below.

  int 
    *ind1 = new int [n1],
    *ind2 = new int [n2];

  double 
    *a1 = new double [n1],
    *a2 = new double [n2];

  CouNumber 
    *Lb = p -> Lb (),
    *Ub = p -> Ub ();

  struct indPosPair *pairs = new struct indPosPair [CoinMax (n1,n2)];

  // re-order ind1 and a1, ind2 and a2 ///////////////////////////////////////////

  for (int i=0; i<n1; i++) {
    pairs [i]. index    = ind1c [i];
    pairs [i]. position = i;
  }

  // This is not very wise: qsort behaves horribly on already sorted
  // arrays (if the pivot element is the last one).

  qsort (pairs, n1, sizeof (indPosPair), compPair);

  for (int i=0; i<n1; i++) {

    int rightpos = pairs [i]. position;

    ind1 [i] = ind1c [rightpos];
    a1   [i] = a1c   [rightpos];
  }

  ///////////////////////////

  for (int i=0; i<n2; i++) {
    pairs [i]. index    = ind2c [i];
    pairs [i]. position = i;
  }

  qsort (pairs, n2, sizeof (indPosPair), compPair);

  for (int i=0; i<n2; i++) {

    int rightpos = pairs [i]. position;

    ind2 [i] = ind2c [rightpos];
    a2   [i] = a2c   [rightpos];
  }

  delete [] pairs;

  // Set multiplier of second constraint, to be used with a2 but not
  // with sa2.

  double signA;

  if (sign < 0) {

    double tmp = u2;
    u2 = -l2;
    l2 = -tmp;

    signA = -1.;

  } else signA = 1.;

#ifdef DEBUG
  printf ("here are the two ineqs (sign=%d):\n", sign);
  printf (" 1: %g <=", l1);             for (int i=0; i<n1; i++) printf (" %+g x%d",         a1 [i], ind1 [i]);
  printf ("<= %g\n 2: %g <=", u1, l2);  for (int i=0; i<n2; i++) printf (" %+g x%d", signA * a2 [i], ind2 [i]); printf ("<= %g\n", u2);
#endif

  threshold
     *alphas   = new threshold  [n1 + n2], // contains all alphas (there might be at most n1+n2-1, but let's be flexible)
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
      curalpha -> sign  = (a1 [i1] < 0.) ? NEG : POS;
      curalpha -> indVar = ind1 [i1];

      cnt++;
      i1++;

    } else if (i2 < n2 && (i1 == n1 || ind2 [i2] < ind1 [i1])) {

      curalpha -> alpha = 0.;
      curalpha -> sign  = (signA * a2 [i2] < 0.) ? NEG : POS;
      curalpha -> indVar = ind2 [i2];

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

	curalpha -> sign = (a1 [i1] > 0.) ? POS : NEG;

      else 

	curalpha -> sign = 
	  (a1 [i1] == a2i2) ? 
	  ((a1 [i1] < 0.) ? NEG : POS) :
	  ((a1 [i1] < 0.) ? DN  : UP);

      curalpha -> indVar = ind1 [i1];

      if (curalpha -> alpha > 0. &&
	  curalpha -> alpha < 1.)
	inalphas [incnt++] = curalpha;

      cnt++;
      i1++;
      i2++;
    }

#ifdef DEBUG
    printf ("(%d,%g,%s) ", 
    	    curalpha -> indVar, 
    	    curalpha -> alpha, 
    	    (curalpha -> sign == NEG) ? "NEG" :
	    (curalpha -> sign == POS) ? "POS" :
	    (curalpha -> sign == DN)  ? "DN"  : "UP");
#endif

    curalpha++;
  }

#ifdef DEBUG
  printf ("\n");
#endif

  // If none of them has an alpha in ]0,1[, nothing needs to be done.

  if (!incnt) {

    delete [] inalphas;
    delete [] alphas;
    delete [] ind1;
    delete [] ind2;
    delete [] a1;
    delete [] a2;

    return 0;
  }

  //----------------------------------------------------------------------
  // Done setting up zeros of (alpha a'_i + (1-alpha) a''_i) x_i
  //----------------------------------------------------------------------

  // Now all thresholds are defined. For all threshold in ]0,1[,
  // determine if they define a convex combination of the two
  // inequalities that allows to tighten the bounds of all variables
  // with at least one nonzero coefficient in the two ineqs.

  int ntightened = 0;

  if (incnt > 1)
    qsort (inalphas, incnt, sizeof (threshold *), compthres);

#ifdef DEBUG
  printf ("--sorted: ");
  for (int i = 0; i < incnt; i++)
    printf ("(%d,%g,%s) ", 
	    inalphas [i] -> indVar, 
	    inalphas [i] -> alpha, 
	    (inalphas [i] -> sign == NEG) ? "NEG" :
	    (inalphas [i] -> sign == POS) ? "POS" :
	    (inalphas [i] -> sign == DN)  ? "DN"  : "UP");    
  printf ("\n");
#endif

  // Compute the max/min of the two constraint bodies. Consider bounds
  // according to the sign of the CONVEX COMBINATION OF THE
  // COEFFICIENT:
  //
  // a_i = alpha a1_i + (1-alpha) a2_i
  //
  // L1_min = sum {i in N} a1_i f-(i)
  // L1_max = sum {i in N} a1_i f+(i)
  //
  // L2_min = sum {i in N} a2_i f-(i)
  // L2_max = sum {i in N} a2_i f+(i)
  //
  // with f-(i) = l_i if a_i > 0 and a*_i > 0
  //                  or a_i < 0 and a*_i < 0
  //              u_i otherwise;
  //
  //      f+(i) = u_i if a_i > 0 and a*_i > 0
  //                  or a_i < 0 and a*_i < 0
  //              l_i otherwise
  //
  // and a* is a1 or a2 depending on the equation.

  double 
    minSum1 = 0., minSum2 = 0.,
    maxSum1 = 0., maxSum2 = 0.;

  int 
    mInfs = 0, // number of -inf summed to minSum 1 and 2 (both)
    pInfs = 0; //           +inf           maxSum

  // The maximum sum in the constraint's body depends on the
  // coefficients arising from the convex combination, but alpha=0
  // at initialization.

  // compute for a_1 ////////////////////////////////////////////////////////

  for (int i=0; i<n1; i++) {

    int indVar1 = ind1 [i];

    double 
      a1i  = a1  [i],
      a2i  = sa2 [indVar1],
      clb1 = clb [indVar1],
      cub1 = cub [indVar1];

    // if no corresponding term on other side, use coefficient of this
    // inequality

    if (a2i == 0.)
      a2i = a1i;

#ifdef DEBUG
    else if (fabs (a2i) < 1e-10)
      printf ("\nNumerics in a2i\n\n");
#endif

    if        (a2i < 0.) {

      bool zeroA = (fabs (sa2 [indVar1]) == 0.);

      if (cub1 >   COUENNE_INFINITY/10) {if (zeroA) mInfs ++;} else minSum1 += a1i * cub1;
      if (clb1 < - COUENNE_INFINITY/10) {if (zeroA) pInfs ++;} else maxSum1 += a1i * clb1;
				                                                         
    } else if (a2i > 0.) {	                                                         
				                                                         
      bool zeroA = (fabs (sa2 [indVar1]) == 0.);

      if (clb1 < - COUENNE_INFINITY/10) {if (zeroA) mInfs ++;} else minSum1 += a1i * clb1;
      if (cub1 >   COUENNE_INFINITY/10) {if (zeroA) pInfs ++;} else maxSum1 += a1i * cub1;
    }
  }

  // compute for a_2 ////////////////////////////////////////////////////////

  for (int i=0; i<n2; i++) {

    int indVar2 = ind2 [i];

    double 
      a2i  = a2  [i] * signA,
      clb2 = clb [indVar2],
      cub2 = cub [indVar2];

    if (a2i < 0.) {

      if (cub2 >   COUENNE_INFINITY/10) mInfs ++; else minSum2 += a2i * cub2;
      if (clb2 < - COUENNE_INFINITY/10) pInfs ++; else maxSum2 += a2i * clb2;

    } else {

      if (clb2 < - COUENNE_INFINITY/10) mInfs ++; else minSum2 += a2i * clb2;
      if (cub2 >   COUENNE_INFINITY/10) pInfs ++; else maxSum2 += a2i * cub2;
    }
  }

  // At this point, m[in,ax]Sum2 contain the finite max/min sums of
  // the second term ax, and pInfs and mInfs the number of bounds
  // that make the min/max sum infinite. We are at alpha=0,
  // therefore this is the max/min sum of the whole convex
  // combination of the two constraint bodies.

#ifdef DEBUG
  printf ("  at alpha=zero, m[in,ax]sum[12] = (1:(%g,%g), 2:(%g,%g)), %d mInf, %d pInf\n",
	  minSum1, maxSum1, 
	  minSum2, maxSum2,
	  mInfs, pInfs);
#endif

  // scan all alphas in ]0,1[ ///////////////////////////////////////////////

  // for all variables x[i], look for tighter lower/upper bounds on x[i]

  std::vector <std::pair <int, double> > 
    newLB, // pairs (indVar, value) of new bounds
    newUB;

  delete [] a1;
  delete [] a2;

#ifdef DEBUG
  printf ("  ");
  for (int j=0; j < n1; j++) printf ("x%d [%g,%g] ", ind1 [j], clb [ind1 [j]], cub [ind1 [j]]); printf ("--- ");
  for (int j=0; j < n2; j++) printf ("x%d [%g,%g] ", ind2 [j], clb [ind2 [j]], cub [ind2 [j]]); printf ("\n");
#endif

  for (int i = 0; i < incnt; i++) {

    threshold *inalpha = inalphas [i];

#ifdef DEBUG
    printf ("  looking at %d: (%d,%g,%s)\n", i,
	    inalphas [i] -> indVar, 
	    inalphas [i] -> alpha, 
	    (inalphas [i] -> sign == NEG) ? "NEG" :
	    (inalphas [i] -> sign == POS) ? "POS" :
	    (inalphas [i] -> sign == DN)  ? "DN"  : "UP"); 

    fflush (stdout);
#endif

    // look at each inalphas, those are the only breakpoints where
    // bounds can improve. For alpha in {0,1} there is already a
    // procedure (single constraint implied bound).

    // check for improved bounds at all variables

    double alpha = inalpha -> alpha;

    // enter next interval, update dynamic max/min sums for this
    // value of alpha

    while (true) {

      inalpha = inalphas [i];

      int
	signalpha = inalpha -> sign,   // either UP or DN (not POS or NEG)
	indVar    = inalpha -> indVar; // the index of the variable whose a_i(alpha) changes sign

      double
	clbi = clb [indVar],
	cubi = cub [indVar];

      // as a variable is changing sign of (alpha a'_i + (1-alpha)
      // a''_i), so must m[in|ax]Sum[12]

      if (fabs (clbi) > 0.) {

	if (clbi < - COUENNE_INFINITY/10) {

	  if (signalpha == DN) {pInfs++; mInfs--;} 
	  else                 {pInfs--; mInfs++;} // sign is UP

	} else {

	  double 
	    sa1ic = sa1 [indVar] * clbi,
	    sa2ic = sa2 [indVar] * clbi;

	  if (signalpha == DN) { // means sa1 [indVar] < 0 while sa2 [indVar] > 0

	    minSum1 -= sa1ic;    minSum2 -= sa2ic;
	    maxSum1 += sa1ic;    maxSum2 += sa2ic;

	  } else {               // sign must be UP

	    minSum1 += sa1ic;    minSum2 += sa2ic;
	    maxSum1 -= sa1ic;    maxSum2 -= sa2ic;
	  }
	}
      }

      // same for upper bound

      if (fabs (cubi) > 0.) {

	if (cubi >   COUENNE_INFINITY/10) {

	  if (signalpha == DN) {pInfs--; mInfs++;} 
	  else                 {pInfs++; mInfs--;} // sign is UP

	} else {

	  double 
	    sa1ic = sa1 [indVar] * cubi,
	    sa2ic = sa2 [indVar] * cubi;

	  if (signalpha == DN) { // means sa1 [indVar] < 0 while sa2 [indVar] > 0

	    minSum1 += sa1ic;    minSum2 += sa2ic;
	    maxSum1 -= sa1ic;    maxSum2 -= sa2ic;

	  } else { // sign must be UP

	    minSum1 -= sa1ic;    minSum2 -= sa2ic;
	    maxSum1 += sa1ic;    maxSum2 += sa2ic;
	  }
	}
      }

      if ((i < incnt-1) && 
	  (inalphas [i+1] -> alpha == alpha)) ++i;
      else break;
    }

#ifdef DEBUG
    printf ("  at alpha=%g, m[in,ax]sum[12] = (1:(%g,%g), 2:(%g,%g)), %d mInf, %d pInf\n",
	    inalphas [i] -> alpha,
	    minSum1, maxSum1, 
	    minSum2, maxSum2,
	    mInfs, pInfs);
#endif

    for (int j=0; j<n1+n2; j++) {

      // Loop on all variables that might improve a bound as a
      // result of the convex combination with the current value of
      // alpha.

      // scan all nonzeros, those on one side and those on both
      // sides. This has already been scanned.
      if (j >= n1 && sa1 [ind2 [j - n1]] != 0.) 
	continue; 

      int
	indVar = ((j < n1) ? ind1 [j] : ind2 [j - n1]),
	tickMin = 0, // one if subtracting this term removes one infinity from minSum
	tickMax = 0; //                                                        max   

      double
	sa1i = sa1 [indVar],
	sa2i = sa2 [indVar],
	ci = alpha * sa1i + (1. - alpha) * sa2i; // combination of coefficients

      // ignore variables whose coefficient is cancelled by this
      // value of alpha

      if (fabs (ci) < MIN_DENOM)
	continue;

      double
	clbi = clb [indVar],
	cubi = cub [indVar],

	newL = - COUENNE_INFINITY,
	newU =   COUENNE_INFINITY,

	subMin1 = 0., subMin2 = 0.,
	subMax1 = 0., subMax2 = 0.;

#ifdef DEBUG
      printf ("  now tightening x%d\n    denominator: %g * %g + %g * %g = %g\n", 
	      (j < n1) ? ind1 [j] : ind2 [j-n1],
	      alpha, sa1i, 1-alpha, sa2i, ci);
#endif

      // m[in,ax]Sum[1,2]A must be updated by subtracting the term
      // that is now at the denominator (no need to multiply by
      // signA, sa2 is already adjusted)

      if (clbi <= - COUENNE_INFINITY/10) {

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

      if (cubi >=   COUENNE_INFINITY/10) {

	// we are deleting an infinite bound from the numerator, and
	// whether it goes towards max- or min- Sum?A it depends on
	// each coefficient

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

      if (fabs (ci) < MIN_DENOM) { // must check that the numerator behaves

      } else {

	if ((l1 > -COUENNE_INFINITY/10) && 
	    (l2 > -COUENNE_INFINITY/10) &&
	    (pInfs == tickMax)) {

	  newL = ((l1 - l2 - (maxSum1 - maxSum2) + (subMax1 - subMax2)) * alpha + l2 - maxSum2 + subMax2) / ci;

#ifdef DEBUG
	  printf ("    \
attempting newL = ((l1 - l2 - (maxSum1 - maxSum2) + (subMax1 - subMax2)) * alpha + l2 - maxSum2 + subMax2) / ci\n    \
((%g - %g - (%g - %g) + (%g - %g)) * %g + %g - %g + %g) / %g = %g\n",
		  l1, l2, maxSum1, maxSum2, subMax1, subMax2, alpha, l2, maxSum2, subMax2, ci, newL);
#endif
	}

	if ((u1 <  COUENNE_INFINITY/10) && 
	    (u2 <  COUENNE_INFINITY/10) &&
	    (mInfs == tickMin)) {

	  newU = ((u1 - u2 - (minSum1 - minSum2) + (subMin1 - subMin2)) * alpha + u2 - minSum2 + subMin2) / ci;

#ifdef DEBUG
	  printf ("    \
attempting newU = ((u1 - u2 - (minSum1 - minSum2) + (subMin1 - subMin2)) * alpha + u2 - minSum2 + subMin2) / ci\n    \
((%g - %g - (%g - %g) + (%g - %g)) * %g + %g - %g + %g) / %g = %g\n",
		  u1, u2, minSum1, minSum2, subMin1, subMin2, alpha, u2, minSum2, subMin2, ci, newU);
#endif
	}

	if (ci < 0.) { // should have done the opposite assignment -- just swap them

#ifdef DEBUG
	  printf ("    swap'em: %g, %g\n", newL, newU);
#endif

	  register double tmp = newL <= - COUENNE_INFINITY / 10 ?   COUENNE_INFINITY : newL;
	  newL                = newU >=   COUENNE_INFINITY / 10 ? - COUENNE_INFINITY : newU;
	  newU                = tmp;
	}

#ifdef DEBUG
	printf ("    final: %g, %g\n", newL, newU);
#endif
      }

#ifdef DEBUG
      printf ("    bounds for x_%d: [%g,%g] vs [%g,%g]\n", indVar, newL, newU, clb [indVar], cub [indVar]);
#endif

      if ((newL > cubi + COUENNE_EPS) || 
	  (newU < clbi - COUENNE_EPS)) {

	delete [] inalphas;
	delete [] alphas;
	delete [] ind1;
	delete [] ind2;

#ifdef DEBUG
	printf ("infeasible\n");
#endif
	return -1;
      }

#ifdef DEBUG
      if (p -> bestSol () &&
	  ((p -> bestSol () [indVar] > newU + COUENNE_EPS) ||
	   (p -> bestSol () [indVar] < newL - COUENNE_EPS)) &&
	  (p -> bestSol () [indVar] >= clbi) &&
	  (p -> bestSol () [indVar] <= cubi))

	printf ("optimum violated: %g not in [%g,%g]\n", p -> bestSol () [indVar], newL, newU);
#endif

      double 
	&lbi = Lb [indVar],
	&ubi = Ub [indVar];

      if ((newL > lbi + COUENNE_EPS) && (newL > -COUENNE_INFINITY / 10)) {
		   
	ntightened++;
	lbi = newL;
	newLB. push_back (std::pair <int, double> (indVar, newL));

#ifdef DEBUG
	printf ("    new lower bound: x%d >= %g [%g]\n", indVar, newL, clb [indVar]);
#endif
      }

      if ((newU < ubi - COUENNE_EPS) && (newU < COUENNE_INFINITY / 10)) {

	ntightened++;
	ubi = newU;
	newUB. push_back (std::pair <int, double> (indVar, newU));

#ifdef DEBUG
	printf ("    new upper bound: x%d <= %g [%g]\n", indVar, newU, cub [indVar]);
#endif
      }
    }
  }

#ifdef DEBUG
  printf ("===================================\n");
#endif

  for (std::vector <std::pair <int, double> >::iterator i = newLB. begin (); i != newLB. end (); ++i) {

    if (isInteger [i -> first]) i -> second = ceil (i -> second - COUENNE_EPS);

    if (i -> second > clb [i -> first]) 
      Lb [i -> first] = clb [i -> first] = i -> second;
  }

  for (std::vector <std::pair <int, double> >::iterator i = newUB. begin (); i != newUB. end (); ++i) {

    if (isInteger [i -> first]) i -> second = floor (i -> second + COUENNE_EPS);

    if (i -> second < cub [i -> first]) 
      Ub [i -> first] = cub [i -> first] = i -> second;
  }

  delete [] inalphas;
  delete [] alphas;

  delete [] ind1;
  delete [] ind2;

  return ntightened;
}

}
