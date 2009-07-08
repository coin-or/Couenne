/* $Id$
 *
 * Name:    CouenneLPtightenBoundsCLP.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: adaptation of OsiClpSolverInterface::tightenBounds()
 *
 * (C) Carnegie-Mellon University, 2009.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"
#include "exprVar.hpp"

//#define COIN_DEVELOP 4

// Tighten bounds. Returns -1 if infeasible, otherwise number of
// variables tightened.
template <class T> 
int CouenneSolverInterface<T>::tightenBoundsCLP (int lightweight) {

  // Copied from OsiClpSolverInterface::tightenBounds

  int
    numberRows    = T::getNumRows(),
    numberColumns = T::getNumCols(),
    iRow, iColumn;

  const double * columnUpper = T::getColUpper();
  const double * columnLower = T::getColLower();
  const double * rowUpper = T::getRowUpper();
  const double * rowLower = T::getRowLower();

  // Column copy of matrix
  const double * element = T::getMatrixByCol()->getElements();
  const int * row = T::getMatrixByCol()->getIndices();
  const CoinBigIndex * columnStart = T::getMatrixByCol()->getVectorStarts();
  const int * columnLength = T::getMatrixByCol()->getVectorLengths();
  const double *objective = T::getObjCoefficients() ;

  double direction = T::getObjSense();
  double * down = new double [numberRows];

  if (lightweight)
    return tightenBoundsCLP_Light (lightweight);

  // NOT LIGHTWEIGHT /////////////////////////////////////////////////////

  double * up = new double [numberRows];
  double * sum = new double [numberRows];
  int * type = new int [numberRows];
  CoinZeroN(down,numberRows);
  CoinZeroN(up,numberRows);
  CoinZeroN(sum,numberRows);
  CoinZeroN(type,numberRows);
  double infinity = T::getInfinity();

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    double lower = columnLower[iColumn];
    double upper = columnUpper[iColumn];
    if (lower==upper) {
      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	double value = element[j];
	sum[iRow]+=2.0*fabs(value*lower);
	if ((type[iRow]&1)==0)
	  down[iRow] += value*lower;
	if ((type[iRow]&2)==0)
	  up[iRow] += value*lower;
      }
    } else {
      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	double value = element[j];
	if (value>0.0) {
	  if ((type[iRow]&1)==0) {
	    if (lower!=-infinity) {
	      down[iRow] += value*lower;
	      sum[iRow]+=fabs(value*lower);
	    } else {
	      type[iRow] |= 1;
	    }
	  }
	  if ((type[iRow]&2)==0) {
	    if (upper!=infinity) {
	      up[iRow] += value*upper;
	      sum[iRow]+=fabs(value*upper);
	    } else {
	      type[iRow] |= 2;
	    }
	  }
	} else {
	  if ((type[iRow]&1)==0) {
	    if (upper!=infinity) {
	      down[iRow] += value*upper;
	      sum[iRow]+=fabs(value*upper);
	    } else {
	      type[iRow] |= 1;
	    }
	  }
	  if ((type[iRow]&2)==0) {
	    if (lower!=-infinity) {
	      up[iRow] += value*lower;
	      sum[iRow]+=fabs(value*lower);
	    } else {
	      type[iRow] |= 2;
	    }
	  }
	}
      }
    }
  }

  int nTightened = 0;
  double tolerance = 1.0e-6;

  for (iRow=0;iRow<numberRows;iRow++) {
    if ((type[iRow]&1)!=0)
      down[iRow]=-infinity;
    if (down[iRow]>rowUpper[iRow]) {
      if (down[iRow]>rowUpper[iRow]+tolerance+1.0e-8*sum[iRow]) {
	// infeasible
#ifdef COIN_DEVELOP
	printf("infeasible on row %d\n",iRow);
#endif
	nTightened=-1;
	break;
      } else {
	down[iRow]=rowUpper[iRow];
      }
    }
    if ((type[iRow]&2)!=0)
      up[iRow]=infinity;
    if (up[iRow]<rowLower[iRow]) {
      if (up[iRow]<rowLower[iRow]-tolerance-1.0e-8*sum[iRow]) {
	// infeasible
#ifdef COIN_DEVELOP
	printf("infeasible on row %d\n",iRow);
#endif
	nTightened=-1;
	break;
      } else {
	up[iRow]=rowLower[iRow];
      }
    }
  }

  if (nTightened)
    numberColumns = 0; // so will skip

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double lower = columnLower[iColumn];
    double upper = columnUpper[iColumn];
    double gap = upper-lower;

    if (!gap)
      continue;

    int canGo=0;

    CoinBigIndex
      start =         columnStart  [iColumn],
      end   = start + columnLength [iColumn];

    if (lower < -1.0e8 && upper > 1.0e8)
      continue; // Could do severe damage to accuracy


    // there was an ifInteger condition here. We do like tightened
    // bounds for continuous variables too, so we don't test for
    // integrality.

    std::vector <exprVar *> &vars = cutgen_ -> Problem () -> Variables ();

    {
      if (vars [iColumn] -> isInteger ()) {

	if (lower < ceil (lower - COUENNE_EPS) - COUENNE_EPS) {
#ifdef COIN_DEVELOP
	  printf("increasing lower bound on %d from %e to %e\n",iColumn,
		 lower,ceil(lower - COUENNE_EPS));
#endif
	  lower=ceil(lower - COUENNE_EPS);
	  gap=upper-lower;
	  T::setColLower(iColumn,lower);
	}

	if (upper > floor(upper + COUENNE_EPS) + COUENNE_EPS) {
#ifdef COIN_DEVELOP
	  printf("decreasing upper bound on %d from %e to %e\n",iColumn,
		 upper,floor(upper + COUENNE_EPS));
#endif
	  upper=floor(upper + COUENNE_EPS);
	  gap=upper-lower;
	  T::setColUpper(iColumn,upper);
	}
      }

      double newLower=lower;
      double newUpper=upper;

      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	double value = element[j];
	if (value>0.0) {
	  if ((type[iRow]&1)==0) {
	    // has to be at most something
	    if (down[iRow] + value*gap > rowUpper[iRow]+tolerance) {
	      double newGap = (rowUpper[iRow]-down[iRow])/value;
	      // adjust
	      newGap += 1.0e-10*sum[iRow];
	      if (vars [iColumn] -> isInteger ())
		newGap = floor(newGap);
	      if (lower+newGap<newUpper)
		newUpper=lower+newGap;
	    }
	  }
	  if (down[iRow]<rowLower[iRow])
	    canGo |=1; // can't go down without affecting result
	  if ((type[iRow]&2)==0) {
	    // has to be at least something
	    if (up[iRow] - value*gap < rowLower[iRow]-tolerance) {
	      double newGap = (up[iRow]-rowLower[iRow])/value;
	      // adjust
	      newGap += 1.0e-10*sum[iRow];
	      if (vars [iColumn] -> isInteger ())
		newGap = floor(newGap);
	      if (upper-newGap>newLower)
		newLower=upper-newGap;
	    }
	  }
	  if (up[iRow]>rowUpper[iRow])
	    canGo |=2; // can't go up without affecting result
	} else {
	  if ((type[iRow]&1)==0) {
	    // has to be at least something
	    if (down[iRow] - value*gap > rowUpper[iRow]+tolerance) {
	      double newGap = -(rowUpper[iRow]-down[iRow])/value;
	      // adjust
	      newGap += 1.0e-10*sum[iRow];
	      if (vars [iColumn] -> isInteger ())
		newGap = floor(newGap);
	      if (upper-newGap>newLower)
		newLower=upper-newGap;
	    }
	  }
	  if (up[iRow]>rowUpper[iRow])
	    canGo |=1; // can't go down without affecting result
	  if ((type[iRow]&2)==0) {
	    // has to be at most something
	    if (up[iRow] + value*gap < rowLower[iRow]-tolerance) {
	      double newGap = -(up[iRow]-rowLower[iRow])/value;
	      // adjust
	      newGap += 1.0e-10*sum[iRow];
	      if (vars [iColumn] -> isInteger ())
		newGap = floor(newGap);
	      if (lower+newGap<newUpper)
		newUpper=lower+newGap;
	    }
	  }
	  if (down[iRow]<rowLower[iRow])
	    canGo |=2; // can't go up without affecting result
	}
      }

      if (newUpper<upper || newLower>lower) {
	nTightened++;
	if (newLower>newUpper) {
	  // infeasible
#if COIN_DEVELOP>1
	  printf("infeasible on column %d\n",iColumn);
#endif
	  nTightened=-1;
	  break;
	} else {
	  T::setColLower(iColumn,newLower);
	  T::setColUpper(iColumn,newUpper);
	}
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  double value = element[j];
	  if (value>0.0) {
	    if ((type[iRow]&1)==0) down [iRow] += value*(newLower-lower);
	    if ((type[iRow]&2)==0) up   [iRow] += value*(newUpper-upper);
	  } else {
	    if ((type[iRow]&1)==0) down [iRow] += value*(newUpper-upper);
	    if ((type[iRow]&2)==0) up   [iRow] += value*(newLower-lower);
	  }
	}
      } else {

	if (canGo!=3) {

	  double objValue = direction*objective[iColumn];

	  if (objValue>=0.0&&(canGo&1)==0) {
#if COIN_DEVELOP>2
	    printf("dual fix down on column %d\n",iColumn);
#endif
	    nTightened++;
	    T::setColUpper(iColumn,lower);
	  } else if (objValue<=0.0 && (canGo&2)==0) {
#if COIN_DEVELOP>2
	    printf("dual fix up on column %d\n",iColumn);
#endif
	    nTightened++;
	    T::setColLower(iColumn,upper);
	  }
	}	    
      }
    }

//     else {

//       // CONTINUOUS //////////////////////////////////////////

//       // just do dual tests
//       for (CoinBigIndex j=start;j<end;j++) {
// 	int iRow = row[j];
// 	double value = element[j];
// 	if (value>0.0) {
// 	  if (down [iRow] < rowLower [iRow]) canGo |=1; // can't go down without affecting result
// 	  if (up   [iRow] > rowUpper [iRow]) canGo |=2; // can't go up   without affecting result
// 	} else {
// 	  if (up   [iRow] > rowUpper [iRow]) canGo |=1; // can't go down without affecting result
// 	  if (down [iRow] < rowLower [iRow]) canGo |=2; // can't go up   without affecting result
// 	}
//       }

//       if (canGo!=3) {
// 	double objValue = direction*objective[iColumn];
// 	if (objValue>=0.0&&(canGo&1)==0) {
// #if COIN_DEVELOP>2
// 	  printf("dual fix down on continuous column %d lower %g\n",
// 		 iColumn,lower);
// #endif
// 	  // Only if won't cause numerical problems
// 	  if (lower>-1.0e10) {
// 	    nTightened++;;
// 	    setColUpper(iColumn,lower);
// 	  }
// 	} else if (objValue<=0.0&&(canGo&2)==0) {
// #if COIN_DEVELOP>2
// 	  printf("dual fix up on continuous column %d upper %g\n",
// 		 iColumn,upper);
// #endif
// 	  // Only if won't cause numerical problems
// 	  if (upper<1.0e10) {
// 	    nTightened++;;
// 	    setColLower(iColumn,upper);
// 	  }
// 	}
//       }
//     }

  }

  delete [] type;
  delete [] down;
  delete [] up;
  delete [] sum;

  return nTightened;
}
  
