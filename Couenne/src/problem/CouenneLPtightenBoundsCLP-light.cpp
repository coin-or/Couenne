/* $Id$
 *
 * Name:    CouenneLPtightenBoundsCLP-light.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: adaptation of OsiClpSolverInterface::tightenBounds() -- light version
 *
 * (C) Carnegie-Mellon University, 2009.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneSolverInterface.hpp"

#define COIN_DEVELOP 4

// Tighten bounds - lightweight. Returns -1 if infeasible, otherwise
// number of variables tightened.
int CouenneSolverInterface::tightenBoundsCLP_Light (int lightweight) {

  // Copied from OsiClpSolverInterface::tightenBounds

  int
    numberRows    = getNumRows(),
    numberColumns = getNumCols();

  const double * columnUpper = getColUpper();
  const double * columnLower = getColLower();
  const double * rowUpper = getRowUpper();
  const double * rowLower = getRowLower();

  // Column copy of matrix
  const double * element = getMatrixByCol()->getElements();
  const int * row = getMatrixByCol()->getIndices();
  const CoinBigIndex * columnStart = getMatrixByCol()->getVectorStarts();
  const int * columnLength = getMatrixByCol()->getVectorLengths();

  double * down = new double [numberRows];

  int * first = new int[numberRows];
  CoinZeroN(first,numberRows);
  CoinZeroN(down,numberRows);
  double * sum = new double [numberRows];
  CoinZeroN(sum,numberRows);
  int numberTightened=0;

  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    double lower = columnLower[iColumn];
    double upper = columnUpper[iColumn];
    if (lower==upper) {
      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	double value = element[j];
	down[iRow] += value*lower;
	sum[iRow] += fabs(value*lower);
      }
    } else {
      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	int n=first[iRow];
	if (n==0&&element[j])
	  first[iRow]=-iColumn-1;
	else if (n<0) 
	  first[iRow]=2;
      }
    }
  }

  double tolerance = 1.0e-6;

  const char * integerInformation = modelPtr_ -> integerInformation ();

  for (int iRow=0;iRow<numberRows;iRow++) {

    int iColumn = first[iRow];

    if (iColumn<0) {

      iColumn = -iColumn-1;

      double lowerRow = rowLower[iRow];
      if (lowerRow>-1.0e20)
	lowerRow -= down[iRow];
      double upperRow = rowUpper[iRow];
      if (upperRow<1.0e20)
	upperRow -= down[iRow];
      double lower = columnLower[iColumn];
      double upper = columnUpper[iColumn];
      double value=0.0;
      for (CoinBigIndex j = columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	if (iRow==row[j]) {
	  value=element[j];
	  break;
	}
      }

      assert (value);

      // convert rowLower and Upper to implied bounds on column

      double
	newLower = -COIN_DBL_MAX,
	newUpper =  COIN_DBL_MAX;

      if (value > 0.0) {
	if (lowerRow > -1.0e20) newLower = lowerRow / value;
	if (upperRow <  1.0e20) newUpper = upperRow / value;
      } else {
	if (upperRow <  1.0e20) newLower = upperRow / value;
	if (lowerRow > -1.0e20) newUpper = lowerRow / value;
      }

      double tolerance2 = 1.0e-6 + 1.0e-8 * sum [iRow];

      if (integerInformation && integerInformation[iColumn]) {

	newLower = (newLower-floor(newLower)<tolerance2) ?
	  floor (newLower) :
	  ceil  (newLower);

	newUpper = (ceil(newUpper)-newUpper<tolerance2) ?
	  ceil  (newUpper) : 
	  floor (newUpper);
      }

      if (newLower>lower+10.0*tolerance2||
	  newUpper<upper-10.0*tolerance2) {
	numberTightened++;
	newLower = CoinMax(lower,newLower);
	newUpper = CoinMin(upper,newUpper);
	if (newLower>newUpper+tolerance) {
	  //printf("XXYY inf on bound\n");
	  numberTightened=-1;
	  break;
	}
	setColLower(iColumn,newLower);
	setColUpper(iColumn,CoinMax(newLower,newUpper));
      }
    }
  }

  delete [] first;
  delete [] down;
  delete [] sum;
  return numberTightened;
}
