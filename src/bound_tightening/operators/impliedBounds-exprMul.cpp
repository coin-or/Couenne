/* $Id$
 *
 * Name:    impliedBounds-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bounds for multiplications
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <stdio.h>
#include <stdlib.h>

#include "CouenneExprMul.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneConfig.h"
#include "CoinFinite.hpp"

//#define FM_CHECK
//#define FM_MOD

using namespace Couenne;


/// implied bound processing for expression w = x*y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprMul::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg, enum auxSign sign) {

  bool resL, resU = resL = false;
  int ind;

  CouNumber 
    wl = sign == expression::AUX_GEQ ? -COIN_DBL_MAX : l [wind],
    wu = sign == expression::AUX_LEQ ?  COIN_DBL_MAX : u [wind];

  if ((arglist_ [ind=0] -> Type () <= CONST) || 
      (arglist_ [ind=1] -> Type () <= CONST)) {

    // at least one constant in product w=cx:
    //
    // wl/c <= x <= wu/c, if c is positive
    // wu/c <= x <= wl/c, if c is negative

    CouNumber c = arglist_ [ind] -> Value ();

    bool argInt = arglist_ [1-ind] -> isInteger ();

    // get the index of the nonconstant part
    ind = arglist_ [1-ind] -> Index ();

    if (ind==-1) // should not happen, it is a product of constants
      return false;

    if        (c >   COUENNE_EPS) {

      resL = (wl > - COUENNE_INFINITY) && updateBound (-1, l + ind, argInt ? ceil  (wl / c - COUENNE_EPS) : (wl / c));
      resU = (wu <   COUENNE_INFINITY) && updateBound ( 1, u + ind, argInt ? floor (wu / c + COUENNE_EPS) : (wu / c));

    } else if (c < - COUENNE_EPS) {

      resL = (wu <   COUENNE_INFINITY) && updateBound (-1, l + ind, argInt ? ceil  (wu / c - COUENNE_EPS) : (wu / c));
      resU = (wl > - COUENNE_INFINITY) && updateBound ( 1, u + ind, argInt ? floor (wl / c + COUENNE_EPS) : (wl / c));
    } 

    if (resL) chg [ind].setLower(t_chg_bounds::CHANGED);
    if (resU) chg [ind].setUpper(t_chg_bounds::CHANGED);

      /*printf ("w_%d [%g,%g] -------> x_%d in [%g,%g] ", 
	      wind, l [wind], u [wind], 
	      ind,  l [ind],  u [ind]);*/
  } else {

    // these bounds would be implied by McCormick's convexification,
    // however we write them explicitly for internal use within bound
    // tightening, as otherwise they would be known to Clp only.

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber
      *xl = l + xi, *yl = l + yi,
      *xu = u + xi, *yu = u + yi;

    // w's lower bound ///////////////////////////////////////////

    bool resxL,  resxU,  resyL, resyU = 
         resxL = resxU = resyL = false;

    bool xlIsZero = (fabs(*xl) < COUENNE_EPS);
    bool xuIsZero = (fabs(*xu) < COUENNE_EPS);
    bool ylIsZero = (fabs(*yl) < COUENNE_EPS);
    bool yuIsZero = (fabs(*yu) < COUENNE_EPS);
    bool wlIsZero = (fabs(wl)  < COUENNE_EPS);
    bool wuIsZero = (fabs(wu)  < COUENNE_EPS);

#ifdef FM_MOD

    if(wlIsZero) {
      if ((!xuIsZero) && (!yuIsZero) && (*xu * *yu < wl)) {
	if (!ylIsZero) {
	  resxU = (*xu * *yl < wl) && updateBound (+1, xu, wl / *yl);
	  xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false);                
	} 
	if (!xlIsZero) {
	  resyU = (*xl * *yu < wl) && updateBound (+1, yu, wl / *xl);
	  yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false);
	}
      }

      // point C in central infeasible area

      if ((!xlIsZero) && (!ylIsZero) && (*xl * *yl < wl)) {
	if (!yuIsZero) {
	  resxL = (*xl * *yu < wl) && updateBound (-1, xl, wl / *yu);
	  xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false); 
	}
	if (!xuIsZero) {
	  resyL = (*xu * *yl < wl) && updateBound (-1, yl, wl / *xu);
	  ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	}
      }
    }
    else // wl is not zero

#endif
    {
      if (wl >= 0.) {
	// point B in central infeasible area
	if (*xu * *yu < wl) {

#ifdef FM_CHECK  
          if((*xu * *yl < wl) && (ylIsZero)) {
	    printf("ylIsZero (A): %g\n", *yl);
	    exit(1);
	  }
          if((*xl * *yu < wl) && (xlIsZero)) {
	    printf("xlIsZero (B): %g\n", *xl);
	    exit(1);
          }
#endif

	  resxU = (*xu * *yl < wl) && updateBound (+1, xu, wl / *yl);
          xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false); 
	  resyU = (*xl * *yu < wl) && updateBound (+1, yu, wl / *xl);
	  yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false);
	}
	
	// point C in central infeasible area
	
	if (*xl * *yl < wl) {

#ifdef FM_CHECK
          if((*xl * *yu < wl) && (yuIsZero)) {
            printf("yuIsZero (C): %g\n", *yu);
            exit(1);
          }
          if((*xu * *yl < wl) && (xuIsZero)) {
            printf("xuIsZero (D): %g\n", *xu);
            exit(1);
          }
#endif
	  resxL = (*xl * *yu < wl) && (!yuIsZero) && updateBound (-1, xl, wl / *yu);
          xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false);
	  resyL = (*xu * *yl < wl) && (!xuIsZero) && updateBound (-1, yl, wl / *xu);
          ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	}
	
      } else if (wl > -COUENNE_INFINITY) {
	
	// the infeasible set is a hyperbola with two branches
	
	// upper left

#ifdef FM_CHECK
        if(((*xl * *yl < wl) && (*yl > 0)) && (ylIsZero)) {
	  printf("ylIsZero (E): %g\n", *yl);
	}
        if(((*xu * *yu < wl) && (*yu > 0)) && (xuIsZero)) {
	  printf("xuIsZero (F): %g\n", *xu);                                  
	  exit(1);
	}
        if(((*xl * *yl < wl) && (*yl < 0)) && (xlIsZero)) {
	  printf("xlIsZero (G): %g\n", *xl);                                  
	  exit(1);
	}
        if(( (*xu * *yu < wl) && (*yu < 0)) && (yuIsZero)) {
	  printf("yuIsZero (H): %g\n", *yu);                                  
	}
#endif

	resxL = (*xl * *yl < wl) && (*yl > 0.) && (!ylIsZero) && updateBound (-1, xl, wl / *yl); // point C
        xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false);
	resyU = (*xu * *yu < wl) && (*yu > 0.) && (!xuIsZero) && updateBound (+1, yu, wl / *xu); // point B
        yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false); 
	
	// lower right
	resyL = (*xl * *yl < wl) && (*yl < 0.) && (!xlIsZero) && updateBound (-1, yl, wl / *xl); // point C
	ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	resxU = (*xu * *yu < wl) && (*yu < 0.) && (!yuIsZero) && updateBound (+1, xu, wl / *yu); // point B
        xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false);
      }
    }

    bool
      xInt = arglist_ [0] -> isInteger (),
      yInt = arglist_ [1] -> isInteger ();
    
    if (resxL && xInt) *xl = ceil  (*xl - COUENNE_EPS);
    if (resxU && xInt) *xu = floor (*xu + COUENNE_EPS);
    if (resyL && yInt) *yl = ceil  (*yl - COUENNE_EPS);
    if (resyU && yInt) *yu = floor (*yu + COUENNE_EPS);
    
    // w's upper bound ///////////////////////////////////////////
#ifdef FM_MOD
    
    if(wuIsZero) {

      if((!xuIsZero) && (!ylIsZero) && (*xu * *yl > wu)) {
	if(!yuIsZero) {
	  resxU = (*xu * *yu > wu) && updateBound (+1, xu, wu / *yu) || resxU;
	  xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false);
	}
	if(!xlIsZero) {
	  resyL = (*xl * *yl > wu) && updateBound (-1, yl, wu / *xl) || resyL;
	  ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	}
      }

      if((!xlIsZero) && (!yuIsZero) && (*xl * *yu > wu)) {
	if(!ylIsZero) {
	  resxL = (*xl * *yl > wu) && updateBound (-1, xl, wu / *yl) || resxL;
	  xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false);
	}
	if(!xuIsZero) {
	  resyU = (*xu * *yu > wu) && updateBound (+1, yu, wu / *xu) || resyU;
	  yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false);
	}
      }
    }

    else // wu is not zero 

#endif /* FM_MOD */

    {
      if (wu >= 0.) {
	if (wu < COUENNE_INFINITY) {
	  // the infeasible set is a hyperbola with two branches

#ifdef FM_CHECK
          if(((*xu * *yl > wu) && (*yl > 0)) && (ylIsZero)) {
            printf("ylIsZero (A2): yl: %g xu: %g wu: %g\n", *yl, *xu, wu);
            exit(1);
          }
          if(((*xl * *yu > wu) && (*yu > 0)) && (xlIsZero)) {
	    printf("xlIsZero (B2): %g\n", *xl);
	    exit(1);
          }
#endif
	  
	  // upper right
	  resxU = ((*xu * *yl > wu) && (*yl > 0.) && (!ylIsZero) && updateBound (+1, xu, wu / *yl)) || resxU; // point D
          xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false);
	  resyU = ((*xl * *yu > wu) && (*yu > 0.) && (!xlIsZero) && updateBound (+1, yu, wu / *xl)) || resyU; // point A
          yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false);
	  
#ifdef FM_CHECK
          if(((*xl * *yu > wu) && (*yu < 0)) && (yuIsZero)) {
            printf("yuIsZero (C2): %g\n", *yu);
           }
          if(((*xu * *yl > wu) && (*yl < 0)) && (xuIsZero)) {
            printf("xuIsZero (D2): %g\n", *xu);
            exit(1);
          }
#endif

	  // lower left
	  resxL = ((*xl * *yu > wu) && (*yu < 0.) && (!yuIsZero) && updateBound (-1, xl, wu / *yu)) || resxL; // point A
          xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false); 
	  resyL = ((*xu * *yl > wu) && (*yl < 0.) && (!xuIsZero) && updateBound (-1, yl, wu / *xu)) || resyL; // point D
          ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	}
	
      } else {
	
	// point D in central infeasible area
	
	if (*xu * *yl > wu) {

#ifdef FM_CHECK
          if((*xu * *yu > wu) && (yuIsZero)) {
            printf("yuIsZero (E2): %g\n", *yu);
            exit(1);
          }
          if((*xl * *yl > wu) && (xlIsZero)) {
            printf("xlIsZero (F2): %g\n", *xl);
            exit(1);
          }
#endif

	  resxU = ((*xu * *yu > wu) && (!yuIsZero) && updateBound (+1, xu, wu / *yu)) || resxU;
          xuIsZero = (fabs(*xu) < COUENNE_EPS ? true : false);
	  resyL = ((*xl * *yl > wu) && (!xlIsZero) && updateBound (-1, yl, wu / *xl)) || resyL;
          ylIsZero = (fabs(*yl) < COUENNE_EPS ? true : false);
	}
	
	// point A in central infeasible area
	
	if (*xl * *yu > wu) {

#ifdef FM_CHECK
          if((*xl * *yl > wu) && (ylIsZero)) {
            printf("ylIsZero (G2): %g\n", *yl);
            exit(1);
          }
          if((*xu * *yu > wu) && (xuIsZero)) {
            printf("xuIsZero (H2): %g\n", *xu);
            exit(1);
          }
#endif

	  resxL = ((*xl * *yl > wu) && (!ylIsZero) && updateBound (-1, xl, wu / *yl)) || resxL;
          xlIsZero = (fabs(*xl) < COUENNE_EPS ? true : false);
	  resyU = ((*xu * *yu > wu) && (!xuIsZero) && updateBound (+1, yu, wu / *xu)) || resyU;
          yuIsZero = (fabs(*yu) < COUENNE_EPS ? true : false);
	}
      }
    }

    // extra integrality check

    if (resxL) {chg [xi].setLower(t_chg_bounds::CHANGED); if (xInt) *xl = ceil  (*xl - COUENNE_EPS);}
    if (resxU) {chg [xi].setUpper(t_chg_bounds::CHANGED); if (xInt) *xu = floor (*xu + COUENNE_EPS);}
    if (resyL) {chg [yi].setLower(t_chg_bounds::CHANGED); if (yInt) *yl = ceil  (*yl - COUENNE_EPS);}
    if (resyU) {chg [yi].setUpper(t_chg_bounds::CHANGED); if (yInt) *yu = floor (*yu + COUENNE_EPS);}

    resL = resxL || resyL;
    resU = resxU || resyU;
  }

  return (resL || resU);
}
