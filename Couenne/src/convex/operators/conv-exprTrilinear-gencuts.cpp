/* $Id$
 *
 * Name:       conv-exprTrilinear-gencuts.cpp.cpp
 * Source:     GNU C++
 * Author:     Sonia Cafieri
 * Purpose:    generate inequalities defining the convex envelope of a 
 *             trilinear monomial
 * History:    Nov 2010 work started
 */

#include "CouenneCutGenerator.hpp"

#include "CouenneTypes.hpp"
#include "CouenneExprMul.hpp"
#include "CouenneExprTrilinear.hpp"
#include "CouenneProblem.hpp"
#include "CouenneExprAux.hpp"

#include <vector>

//#define DEBUG
using namespace Couenne;

#define EPSILONT 1.e-6

//typedef CouNumber double;

//// permutations of 3 elements
void permutation3(int **ind,int *ibnd)
{
  ind[0][0] = ibnd[0]; ind[0][1] = ibnd[1]; ind[0][2] = ibnd[2];
  ind[1][0] = ibnd[0]; ind[1][1] = ibnd[2]; ind[1][2] = ibnd[1];
  ind[2][0] = ibnd[1]; ind[2][1] = ibnd[0]; ind[2][2] = ibnd[2];
  ind[3][0] = ibnd[1]; ind[3][1] = ibnd[2]; ind[3][2] = ibnd[0];
  ind[4][0] = ibnd[2]; ind[4][1] = ibnd[0]; ind[4][2] = ibnd[1];
  ind[5][0] = ibnd[2]; ind[5][1] = ibnd[1]; ind[5][2] = ibnd[0];
}


/// generate convexification cuts for constraint w = x*y*z
void TriLinCuts (double *vlb, double *vub, int *varIndices,
		 std::vector <std::vector <int> >    &cutIndices,
		 std::vector <std::vector <double> > &cutCoeff,
		 std::vector <double>                &cutLb,
		 std::vector <double>                &cutUb) {

  // var indices
  int v1, v2, v3, v4;
  // number of cuts
  int defcons_size = 20;
  // bounds on cuts
  double *bnd = new double [12];    

  CouNumber cf = 1.;

  int **ind;
  ind = new int*[6];
  for(int i=0; i<6; i++) {
    ind[i] = new int[6];
  }
     
  int *ibnd; 
  ibnd = new int[3];
  ibnd[0] = varIndices[0]; ibnd[1] = varIndices[1]; ibnd[2] = varIndices[2];
#ifdef DEBUG
std::cout << "ibnd[0] =" << ibnd[0] << "  ibnd[1] =" << ibnd[1] << "  ibnd[2] =" << ibnd[2] << std::endl;
std::cout << "vlb[ibnd[0]] =" << vlb[ibnd[0]] << "  vub[ibnd[0]] =" << vub[ibnd[0]] << std::endl;
std::cout << "vlb[ibnd[1]] =" << vlb[ibnd[1]] << "  vub[ibnd[1]] =" << vub[ibnd[1]] << std::endl;
std::cout << "vlb[ibnd[2]] =" << vlb[ibnd[2]] << "  vub[ibnd[2]] =" << vub[ibnd[2]] << std::endl;
#endif
 
  // compute the 6 permutations of the 3 variables 
  permutation3(ind,ibnd);

  int i, flag=0, idx=0;
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] <=0 && vub[ind[i][2]] >=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 1;  // this is case 1
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
       && vub[ind[i][1]] >=0 && vub[ind[i][2]] >=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 2;  // this is case 2
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
       && vub[ind[i][0]] >=0 && vub[ind[i][1]] >=0 && vub[ind[i][2]] >=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 3;  // this is case 3
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
       && vub[ind[i][1]] >=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 4;  // this is case 4
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
       && vub[ind[i][0]] >=0 && vub[ind[i][1]] >=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 5;  // this is case 5
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] <=0 && vub[ind[i][0]] >=0 && vub[ind[i][1]] <=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 6;  // this is case 6
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] >=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 7;  // this is case 7
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 8;  // this is case 8
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vub[ind[i][1]] <=0 
       && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 9;  // this is case 9
    }
    i++; 
  }
  i = 0;
  while(i < 6 && flag == 0) {
    if(vlb[ind[i][0]] <=0 && vub[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vub[ind[i][1]] <=0 
       && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
      idx = i;   // store the index of the permutation satisfying the condition
      flag = 10;  // this is case 10
    }
    i++; 
  }

  if (flag==0) {
    std::cout << "ERROR: case not implemented" << std::endl;
    exit(0);
  }

  // var indices
  v1 = ind[idx][0]; 
  v2 = ind[idx][1]; 
  v3 = ind[idx][2];
  v4 = varIndices [3];

  // lower and upper bound on variables
  double xL1 = cf*vlb[v1];  double  xU1 = cf*vub[v1];
  double xL2 = vlb[v2];  double xU2 = vub[v2];
  double xL3 = vlb[v3];  double xU3 = vub[v3];

#define prepareVectors(a) {			        \
							\
    int size = (int) (cutIndices.size ());		\
      							\
    for (int i=0; i<a; i++) {				\
							\
      cutIndices. push_back (std::vector <int>    ());	\
      cutCoeff.   push_back (std::vector <double> ());	\
      							\
      cutLb.      push_back (-COUENNE_INFINITY);	\
      cutUb.      push_back ( COUENNE_INFINITY);	\
							\
      for (int j=0; j<4; j++) {				\
							\
	cutIndices [size+i].push_back (-1);	       	\
	cutCoeff   [size+i].push_back (0.);		\
      }							\
    }							\
}

  /*----------------------------------------------------------------------------------------*/

  // case 1
  if(flag == 1) {
#ifdef DEBUG
    std::cout << " -- case 1 --" << std::endl;
#endif

    double theta  = xL1*xU2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3; 
    double theta1 = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3; 

    defcons_size = 12;

    prepareVectors(defcons_size);

    for(int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     

      cutCoeff [ii][3] = 1.;
    }
 
    cutCoeff [0][0] = -xU2*xU3; cutCoeff [0][1] = -xU1*xU3; cutCoeff [0][2] = -xU1*xU2;  bnd[0] = - 2.*xU1*xU2*xU3;
    cutCoeff [1][0] = -xU2*xL3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xU2;  bnd[1] = - xL1*xU2*xL3 - xL1*xU2*xU3;
    cutCoeff [2][0] = -xU2*xL3; cutCoeff [2][1] = -xL1*xL3; cutCoeff [2][2] = -xL1*xL2;  bnd[2] = - xL1*xU2*xL3 - xL1*xL2*xL3;
    cutCoeff [3][0] = -xL2*xU3; cutCoeff [3][1] = -xU1*xL3; cutCoeff [3][2] = -xU1*xL2;  bnd[3] = - xU1*xL2*xU3 - xU1*xL2*xL3; 
    cutCoeff [4][0] = -xL2*xL3; cutCoeff [4][1] = -xU1*xL3; cutCoeff [4][2] = -xL1*xL2;  bnd[4] = - xU1*xL2*xL3 - xL1*xL2*xL3;
    cutCoeff [5][0] = -xL2*xU3; cutCoeff [5][1] = -xL1*xU3; cutCoeff [5][2] = -(theta/(xU3-xL3));
    bnd[5] = (-(theta*xL3)/(xU3-xL3)) - xL1*xU2*xU3 - xU1*xL2*xU3 + xU1*xU2*xL3; 

    cutCoeff [6][0] = -xU2*xL3; cutCoeff [6][1] = -xU1*xL3; cutCoeff [6][2] = -xU1*xU2;  bnd[6] = - 2.*xU1*xU2*xL3;
    cutCoeff [7][0] = -xL2*xL3; cutCoeff [7][1] = -xU1*xU3; cutCoeff [7][2] = -xU1*xL2;  bnd[7] = - xU1*xL2*xU3 - xU1*xL2*xL3;
    cutCoeff [8][0] = -xU2*xU3; cutCoeff [8][1] = -xL1*xU3; cutCoeff [8][2] = -xL1*xL2;  bnd[8] = - xL1*xU2*xU3 - xL1*xL2*xU3;
    cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xL1*xL3; cutCoeff [9][2] = -xL1*xU2;  bnd[9] = - xL1*xU2*xU3 - xL1*xU2*xL3;
    cutCoeff [10][0] = -xL2*xU3; cutCoeff [10][1] = -xU1*xU3; cutCoeff [10][2] = -xL1*xL2;  bnd[10] = - xU1*xL2*xU3 - xL1*xL2*xU3;
    cutCoeff [11][0] = -xL2*xL3; cutCoeff [11][1] = -xL1*xL3; cutCoeff [11][2] = -(theta1/(xL3-xU3)); cutCoeff [11][3] = 1.;
    bnd[11] = (-(theta1*xU3)/(xL3-xU3)) - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;
  } // end if case 1

  /*----------------------------------------------------------------------------------------*/

  // case 2
  if(flag == 2) {
#ifdef DEBUG
    std::cout << " -- case 2 --" << std::endl;
#endif

    defcons_size = 12;

    prepareVectors (defcons_size);
 
    // compute the 6 permutations of the 3 variables 
    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
    permutation3(ind,ibnd);
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 6 && flagg == 0) {
      if(vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
	 <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]])
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
      std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

    double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
    double xL2(vlb[v2]); double xU2(vub[v2]);
    double xL3(vlb[v3]); double xU3(vub[v3]);

    for(int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     
      cutCoeff [ii][3] = 1.;
    }

    double theta1 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;
    double theta2 = xU1*xL2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xU2*xU3;

    cutCoeff [0][0] = -xU2*xU3; cutCoeff [0][1] = -xU1*xU3; cutCoeff [0][2] = -xU1*xU2;  bnd[0] = - 2.*xU1*xU2*xU3; 
    cutCoeff [1][0] = -xL2*xL3; cutCoeff [1][1] = -xU1*xL3; cutCoeff [1][2] = -xU1*xL2;  bnd[1] = - 2.*xU1*xL2*xL3;
    cutCoeff [2][0] = -xU2*xL3; cutCoeff [2][1] = -xL1*xU3; cutCoeff [2][2] = -xL1*xU2;  bnd[2] = - xL1*xU2*xL3 - xL1*xU2*xU3;
    cutCoeff [3][0] = -xU2*xL3; cutCoeff [3][1] = -xL1*xL3; cutCoeff [3][2] = -xL1*xL2;  bnd[3] = - xL1*xU2*xL3 - xL1*xL2*xL3;
    cutCoeff [4][0] = -xL2*xU3; cutCoeff [4][1] = -(theta1/(xL2-xU2)); cutCoeff [4][2] = -xL1*xL2;  
    bnd[4] = (-(theta1*xU2)/(xL2-xU2)) - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;
    cutCoeff [5][0] = -xL2*xU3; cutCoeff [5][1] = -xL1*xU3; cutCoeff [5][2] = -(theta2/(xU3-xL3));
    bnd[5] = (-(theta2*xL3)/(xU3-xL3)) - xU1*xL2*xU3 - xL1*xU2*xU3 + xU1*xU2*xL3;
      
    if (  vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
	  >= vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]) {

      double theta1c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
      double theta2c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xL2*xU3;

      cutCoeff [6][0] = -xL2*xU3; cutCoeff [6][1] = -xU1*xU3; cutCoeff [6][2] = -xU1*xL2;  bnd[6] = - 2.*xU1*xL2*xU3;
      cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xU1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xL3;
      cutCoeff [8][0] = -xU2*xU3; cutCoeff [8][1] = -xL1*xU3; cutCoeff [8][2] = -xL1*xL2;  bnd[8] = - xL1*xU2*xU3 - xL1*xL2*xU3;
      cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xL1*xL3; cutCoeff [9][2] = -xL1*xU2;  bnd[9] = - xL1*xU2*xU3 - xL1*xU2*xL3;
      cutCoeff [10][0] = -xL2*xL3; cutCoeff [10][1] = -xL1*xL3; cutCoeff [10][2] = -(theta1c/(xL3-xU3));  
      bnd[10] = (-(theta1c*xU3)/(xL3-xU3)) - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;
      cutCoeff [11][0] = -xL2*xL3; cutCoeff [11][1] = -(theta2c/(xL2-xU2)); cutCoeff [11][2] = -xL1*xL2;
      bnd[11] = (-(theta2c*xU2)/(xL2-xU2)) - xU1*xL2*xL3 - xL1*xL2*xU3 + xU1*xU2*xU3;

    } else {
      double theta1c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
      double theta2c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

      cutCoeff [6][0] = -xL2*xU3; cutCoeff [6][1] = -xU1*xU3; cutCoeff [6][2] = -xU1*xL2;  bnd[6] = - 2.*xU1*xL2*xU3;
      cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xU1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xL3;
      cutCoeff [8][0] = -xL2*xL3; cutCoeff [8][1] = -xL1*xL3; cutCoeff [8][2] = -xL1*xU2;  bnd[8] = - xL1*xL2*xL3 - xL1*xU2*xL3;
      cutCoeff [9][0] = -xL2*xL3; cutCoeff [9][1] = -xL1*xU3; cutCoeff [9][2] = -xL1*xL2;  bnd[9] = - xL1*xL2*xL3 - xL1*xL2*xU3; 
      cutCoeff [10][0] = -xU2*xU3; cutCoeff [10][1] = -xL1*xU3; cutCoeff [10][2] = -(theta1c/(xU3-xL3));  
      bnd[10] = (-(theta1c*xL3)/(xU3-xL3)) - xU1*xU2*xU3 - xL1*xL2*xU3 + xU1*xL2*xL3; 
      cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -(theta2c/(xU2-xL2)); cutCoeff [11][2] = -xL1*xU2;
      bnd[11] = (-(theta2c*xL2)/(xU2-xL2)) - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3; 
    }
      
  } // end if case 2

  /*----------------------------------------------------------------------------------------*/

  // case 3
  if(flag == 3) {
#ifdef DEBUG
    std::cout << " -- case 3 --" << std::endl;
#endif

    int last;

    if(vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3]
       <= vlb[v1]*vlb[v2]*vlb[v3] + 2.*vub[v1]*vub[v2]*vub[v3] &&
       vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]
       <= vlb[v1]*vub[v2]*vub[v3] + 2.*vub[v1]*vlb[v2]*vlb[v3] &&
       vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
       <= vub[v1]*vlb[v2]*vub[v3] + 2.*vlb[v1]*vub[v2]*vlb[v3] &&
       vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3]
       <= vub[v1]*vub[v2]*vlb[v3] + 2.*vlb[v1]*vlb[v2]*vub[v3] ) {
    
      double theta3x = 0.5*(xL1*xU2*xU3 + xL1*xL2*xL3 - xU1*xU2*xL3 - xU1*xL2*xU3)/(xL1-xU1);
      double theta3y = 0.5*(xU1*xL2*xU3 + xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xU2*xU3)/(xL2-xU2);
      double theta3z = 0.5*(xU1*xU2*xL3 + xL1*xL2*xL3 - xU1*xL2*xU3 - xL1*xU2*xU3)/(xL3-xU3);
      double theta3c = xL1*xL2*xL3 - theta3x*xL1 - theta3y*xL2 - theta3z*xL3;

      defcons_size = 5;
      prepareVectors (defcons_size);

      for(int ii = 0; ii < 5; ii++) {

	cutIndices [ii][0] = v1; 
	cutIndices [ii][1] = v2; 
	cutIndices [ii][2] = v3; 
	cutIndices [ii][3] = v4;     

	cutCoeff [ii][3] = 1.;
      }

      cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xL3;
      cutCoeff [1][0] = -xL2*xU3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xL2;  bnd[1] = - 2.*xL1*xL2*xU3;
      cutCoeff [2][0] = -xU2*xU3; cutCoeff [2][1] = -xU1*xU3; cutCoeff [2][2] = -xU1*xU2;  bnd[2] = - 2.*xU1*xU2*xU3;
      cutCoeff [3][0] = -xL2*xL3; cutCoeff [3][1] = -xU1*xL3; cutCoeff [3][2] = -xU1*xL2;  bnd[3] = - 2.*xU1*xL2*xL3;
      cutCoeff [4][0] = -theta3x; cutCoeff [4][1] = -theta3y; cutCoeff [4][2] = -theta3z;  bnd[4] = theta3c;
 
      last=4;

    } else if (vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3] 
	       >= vlb[v1]*vlb[v2]*vlb[v3] + 2.*vub[v1]*vub[v2]*vub[v3]) {
#ifdef DEBUG
std::cout << "else if " << std::endl;
#endif
         
      double theta1 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
      double theta2 = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;
      double theta3 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

      defcons_size = 6;
      prepareVectors (defcons_size);

      for(int ii = 0; ii < 6; ii++) {

	cutIndices [ii][0] = v1; 
	cutIndices [ii][1] = v2; 
	cutIndices [ii][2] = v3; 
	cutIndices [ii][3] = v4;     

	cutCoeff [ii][3] = 1.;
      }

      cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xL3;
      cutCoeff [1][0] = -xL2*xU3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xL2;  bnd[1] = - 2.*xL1*xL2*xU3;
      cutCoeff [2][0] = -xL2*xL3; cutCoeff [2][1] = -xU1*xL3; cutCoeff [2][2] = -xU1*xL2;  bnd[2] = - 2.*xU1*xL2*xL3;
      cutCoeff [3][0] = -(theta1/(xU1-xL1)); cutCoeff [3][1] = -xU1*xU3; cutCoeff [3][2] = -xU1*xU2;  
      bnd[3] = (-(theta1*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;
      cutCoeff [4][0] = -xU2*xU3; cutCoeff [4][1] = -xU1*xU3; cutCoeff [4][2] = -(theta2/(xU3-xL3));  
      bnd[4] = (-(theta2*xL3)/(xU3-xL3)) - xU1*xL2*xU3 - xL1*xU2*xU3 + xL1*xL2*xL3;
      cutCoeff [5][0] = -xU2*xU3; cutCoeff [5][1] = -(theta3/(xU2-xL2)); cutCoeff [5][2] = -xU1*xU2;
      bnd[5] = (-(theta3*xL2)/(xU2-xL2)) - xU1*xU2*xL3 - xL1*xU2*xU3 + xL1*xL2*xL3;

      last=5;

    } else {
#ifdef DEBUG
std::cout << "else " << std::endl;
#endif
      // compute the 6 permutations of the 3 variables 
      ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
      permutation3(ind,ibnd);
      int i, flagg=0, idx=0;
      i = 0;
      while(i < 6 && flagg == 0) {
	if (((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] 
	     >= vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] + 2.*vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]]) ))
	  {
	    idx = i;   // store the index of the permutation satisfying the condition
	    flagg = 1;  // condition is satisfied
	  }
	i++; 
      }
      if (flagg==0) {
	std::cout << "ERROR!!!" << std::endl; exit(0);
      }
      v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

      double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
      double xL2(vlb[v2]); double xU2(vub[v2]);
      double xL3(vlb[v3]); double xU3(vub[v3]);

      //} else if (vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] 
      //         >= vlb[v1]*vub[v2]*vub[v3] + 2.*vub[v1]*vlb[v2]*vlb[v3]) {

      defcons_size = 6;
      prepareVectors (defcons_size);

      for(int ii = 0; ii < 6; ii++) {

	cutIndices [ii][0] = v1; 
	cutIndices [ii][1] = v2; 
	cutIndices [ii][2] = v3; 
	cutIndices [ii][3] = v4;     

	cutCoeff [ii][3] = 1.;
      }

      double theta1 = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
      double theta2 = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
      double theta3 = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xU2*xL3;

      cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xL3;
      cutCoeff [1][0] = -xL2*xU3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xL2;  bnd[1] = - 2.*xL1*xL2*xU3;
      cutCoeff [2][0] = -xU2*xU3; cutCoeff [2][1] = -xU1*xU3; cutCoeff [2][2] = -xU1*xU2;  bnd[2] = - 2.*xU1*xU2*xU3;
      cutCoeff [3][0] = -xL2*xL3; cutCoeff [3][1] = -(theta1/(xL2-xU2)); cutCoeff [3][2] = -xU1*xL2;  
      bnd[3] = (-(theta1*xU2)/(xL2-xU2)) - xL1*xL2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;
      cutCoeff [4][0] = -(theta2/(xU1-xL1)); cutCoeff [4][1] = -xU1*xL3; cutCoeff [4][2] = -xU1*xL2;  
      bnd[4] = (-(theta2*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;
      cutCoeff [5][0] = -xL2*xL3; cutCoeff [5][1] = -xU1*xL3; cutCoeff [5][2] = -(theta3/(xL3-xU3));
      bnd[5] = (-(theta3*xU3)/(xL3-xU3)) - xL1*xL2*xL3 - xU1*xU2*xL3 + xL1*xU2*xU3;

      last=5;
    }

    if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
       >= vub[v1]*vub[v2]*vub[v3] + 2.*vlb[v1]*vlb[v2]*vlb[v3] &&
       vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
       >= vub[v1]*vlb[v2]*vlb[v3] + 2.*vlb[v1]*vub[v2]*vub[v3] &&
       vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
       >= vlb[v1]*vub[v2]*vlb[v3] + 2.*vub[v1]*vlb[v2]*vub[v3] &&
       vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       >= vlb[v1]*vlb[v2]*vub[v3] + 2.*vub[v1]*vub[v2]*vlb[v3] ) {
   
#ifdef DEBUG
std::cout << "2 - if " << std::endl;
#endif
      double theta3x = 0.5*(xU1*xL2*xL3 + xU1*xU2*xU3 - xL1*xL2*xU3 - xL1*xU2*xL3)/(xU1-xL1);
      double theta3y = 0.5*(xL1*xU2*xL3 + xU1*xU2*xU3 - xL1*xL2*xU3 - xU1*xL2*xL3)/(xU2-xL2);
      double theta3z = 0.5*(xL1*xL2*xU3 + xU1*xU2*xU3 - xL1*xU2*xL3 - xU1*xL2*xL3)/(xU3-xL3);
      double theta3c = xU1*xU2*xU3 - theta3x*xU1 - theta3y*xU2 - theta3z*xU3;

      prepareVectors (5);

      for(int ii = last+1; ii <= last+5; ii++) {

	cutIndices [ii][0] = v1; 
        cutIndices [ii][1] = v2; 
        cutIndices [ii][2] = v3; 
        cutIndices [ii][3] = v4;     
	cutCoeff [ii][3] = 1.;
      }

      cutCoeff [last+1][0] = -xL2*xL3; cutCoeff [last+1][1] = -xL1*xL3; cutCoeff [last+1][2] = -xL1*xL2;  bnd[last+1] = - 2.*xL1*xL2*xL3;
      cutCoeff [last+2][0] = -xU2*xL3; cutCoeff [last+2][1] = -xU1*xL3; cutCoeff [last+2][2] = -xU1*xU2;  bnd[last+2] = - 2.*xU1*xU2*xL3;
      cutCoeff [last+3][0] = -xL2*xU3; cutCoeff [last+3][1] = -xU1*xU3; cutCoeff [last+3][2] = -xU1*xL2;  bnd[last+3] = - 2.*xU1*xL2*xU3;
      cutCoeff [last+4][0] = -xU2*xU3; cutCoeff [last+4][1] = -xL1*xU3; cutCoeff [last+4][2] = -xL1*xU2;  bnd[last+4] = - 2.*xL1*xU2*xU3;
      cutCoeff [last+5][0] = -theta3x; cutCoeff [last+5][1] = -theta3y; cutCoeff [last+5][2] = -theta3z;  bnd[last+5] = theta3c;

      defcons_size = last+6;

    } else if (vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] 
	       <= vub[v1]*vub[v2]*vub[v3] + 2.*vlb[v1]*vlb[v2]*vlb[v3]) {
#ifdef DEBUG
std::cout << "2 - else if" << std::endl;
#endif

      double theta1 = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
      double theta2 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
      double theta3 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xU1*xL2*xL3;

      prepareVectors (6);

      for(int ii = last+1; ii <= last+6; ii++) {

	cutIndices [ii][0] = v1; 
	cutIndices [ii][1] = v2; 
	cutIndices [ii][2] = v3; 
	cutIndices [ii][3] = v4;     

	cutCoeff [ii][3] = 1.;
      }
        
      cutCoeff [last+1][0] = -xU2*xL3; cutCoeff [last+1][1] = -xU1*xL3; cutCoeff [last+1][2] = -xU1*xU2;  bnd[last+1] = - 2.*xU1*xU2*xL3;
      cutCoeff [last+2][0] = -xL2*xU3; cutCoeff [last+2][1] = -xU1*xU3; cutCoeff [last+2][2] = -xU1*xL2;  bnd[last+2] = - 2.*xU1*xL2*xU3;       
      cutCoeff [last+3][0] = -xU2*xU3; cutCoeff [last+3][1] = -xL1*xU3; cutCoeff [last+3][2] = -xL1*xU2;  bnd[last+3] = - 2.*xL1*xU2*xU3;
      cutCoeff [last+4][0] = -xL2*xL3; cutCoeff [last+4][1] = -xL1*xL3; cutCoeff [last+4][2] = -(theta1/(xL3-xU3));  
      bnd[last+4] = (-(theta1*xU3)/(xL3-xU3)) - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;
      cutCoeff [last+5][0] = -(theta2/(xL1-xU1)); cutCoeff [last+5][1] = -xL1*xL3; cutCoeff [last+5][2] = -xL1*xL2;
      bnd[last+5] = (-(theta2*xU1)/(xL1-xU1)) - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xU2*xU3;
      cutCoeff [last+6][0] = - xL2*xL3; cutCoeff [last+6][1] = -(theta3/(xL2-xU2)); cutCoeff [last+6][2] = -xL1*xL2;
      bnd[last+6] = (-(theta3*xU2)/(xL2-xU2)) - xL1*xL2*xU3 - xU1*xL2*xL3 + xU1*xU2*xU3;

      defcons_size = last+7;

    } else //if (vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
	   //    <= vub[v1]*vlb[v2]*vlb[v3] + 2.*vlb[v1]*vub[v2]*vub[v3]) 
      {
#ifdef DEBUG
std::cout << "2 - another else if" << std::endl;
std::cout << "v1 = " << v1 << " v2 =" << v2 << "  v3 =" << v3 << std::endl;
#endif
      // compute the 6 permutations of the 3 variables 
      //ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
      permutation3(ind,ibnd);
      int i, flagg=0, idx=0;
      i = 0;
      while(i < 6 && flagg == 0) {
        if (vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
               <= vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + 2.*vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]])
	  {
	    idx = i;   // store the index of the permutation satisfying the condition
	    flagg = 1;  // condition is satisfied
	  }
	i++; 
      }
      if (flagg==0) {
	std::cout << "ERROR!!!" << std::endl; exit(0);
      }
      v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

      double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
      double xL2(vlb[v2]); double xU2(vub[v2]);
      double xL3(vlb[v3]); double xU3(vub[v3]);

      double theta1 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;
      double theta2 = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;
      double theta3 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

      prepareVectors (6);

      for(int ii = last+1; ii <= last+6; ii++) {

	cutIndices [ii][0] = v1; 
	cutIndices [ii][1] = v2; 
	cutIndices [ii][2] = v3; 
	cutIndices [ii][3] = v4;     

	cutCoeff [ii][3] = 1.;
      }
       
      cutCoeff [last+1][0] = -xL2*xL3; cutCoeff [last+1][1] = -xL1*xL3; cutCoeff [last+1][2] = -xL1*xL2;  bnd[last+1] = - 2.*xL1*xL2*xL3;
      cutCoeff [last+2][0] = -xU2*xL3; cutCoeff [last+2][1] = -xU1*xL3; cutCoeff [last+2][2] = -xU1*xU2;  bnd[last+2] = - 2.*xU1*xU2*xL3;
      cutCoeff [last+3][0] = -xL2*xU3; cutCoeff [last+3][1] = -xU1*xU3; cutCoeff [last+3][2] = -xU1*xL2;  bnd[last+3] = - 2.*xU1*xL2*xU3;
      cutCoeff [last+4][0] = -(theta1/(xL1-xU1)); cutCoeff [last+4][1] = -xL1*xU3; cutCoeff [last+4][2] = -xL1*xU2;  
      bnd[last+4] = (-(theta1*xU1)/(xL1-xU1)) - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3; 
      cutCoeff [last+5][0] = -xU2*xU3; cutCoeff [last+5][1] = -(theta2/(xU2-xL2)); cutCoeff [last+5][2] = -xL1*xU2;  
      bnd[last+5] = (-(theta2*xL2)/(xU2-xL2)) - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3; 
      cutCoeff [last+6][0] = -xU2*xU3; cutCoeff [last+6][1] = -xL1*xU3; cutCoeff [last+6][2] = -(theta3/(xU3-xL3));  
      bnd[last+6] = (-(theta3*xL3)/(xU3-xL3)) - xL1*xL2*xU3 - xU1*xU2*xU3 + xU1*xL2*xL3;

      defcons_size = last+7;
    } 
    
  } // end if case 3

  /*----------------------------------------------------------------------------------------*/

  // case 4
  if(flag == 4) {
#ifdef DEBUG
    std::cout << " -- case 4 --" << std::endl;
#endif

    double theta  = xU1*xL2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xL2*xL3;
    double theta1 = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

    defcons_size = 12;

    prepareVectors (defcons_size);

    for(int ii = 0; ii < defcons_size; ii++) {
      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     
      cutCoeff [ii][3] = 1.;
    }

    cutCoeff [0][0] = -xL2*xL3; cutCoeff [0][1] = -xU1*xL3; cutCoeff [0][2] = -xU1*xL2;  bnd[0] = - 2.*xU1*xL2*xL3;
    cutCoeff [1][0] = -xU2*xL3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xU2;  bnd[1] = - xL1*xU2*xL3 - xL1*xU2*xU3;
    cutCoeff [2][0] = -xU2*xL3; cutCoeff [2][1] = -xL1*xL3; cutCoeff [2][2] = -xL1*xL2;  bnd[2] = - xL1*xU2*xL3 - xL1*xL2*xL3;
    cutCoeff [3][0] = -xL2*xU3; cutCoeff [3][1] = -xU1*xU3; cutCoeff [3][2] = -xU1*xU2;  bnd[3] = - xU1*xL2*xU3 - xU1*xU2*xU3;
    cutCoeff [4][0] = -xU2*xU3; cutCoeff [4][1] = -xL1*xU3; cutCoeff [4][2] = -xU1*xU2;  bnd[4] = - xL1*xU2*xU3 - xU1*xU2*xU3; 
    cutCoeff [5][0] = -xL2*xU3; cutCoeff [5][1] = -(theta/(xL2-xU2)); cutCoeff [5][2] = -xL1*xL2;
    bnd[5] = (-(theta*xU2)/(xL2-xU2)) - xU1*xL2*xU3 - xL1*xL2*xL3 + xU1*xU2*xL3;
 
    cutCoeff [6][0] = -xU2*xL3; cutCoeff [6][1] = -xU1*xL3; cutCoeff [6][2] = -xU1*xU2;  bnd[6] = - 2.*xU1*xU2*xL3;
    cutCoeff [7][0] = -xU2*xU3; cutCoeff [7][1] = -xU1*xU3; cutCoeff [7][2] = -xU1*xL2;  bnd[7] = - xU1*xU2*xU3 - xU1*xL2*xU3;
    cutCoeff [8][0] = -xL2*xL3; cutCoeff [8][1] = -xL1*xL3; cutCoeff [8][2] = -xL1*xU2;  bnd[8] = - xL1*xU2*xL3 - xL1*xL2*xL3;
    cutCoeff [9][0] = -xL2*xL3; cutCoeff [9][1] = -xL1*xU3; cutCoeff [9][2] = -xL1*xL2;  bnd[9] = - xL1*xL2*xU3 - xL1*xL2*xL3;
    cutCoeff [10][0] = -xL2*xU3; cutCoeff [10][1] = -xL1*xU3; cutCoeff [10][2] = -xU1*xL2;  bnd[10] = - xL1*xL2*xU3 - xU1*xL2*xU3; 
    cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -(theta1/(xU2-xL2)); cutCoeff [11][2] = -xL1*xU2;
    bnd[11] = (-(theta1*xL2)/(xU2-xL2)) - xL1*xU2*xL3 - xU1*xU2*xU3 + xU1*xL2*xL3;
      
  } // end if case 4

  /*----------------------------------------------------------------------------------------*/

  // case 5
  if(flag == 5) {
#ifdef DEBUG
    std::cout << " -- case 5 --" << std::endl;
 std::cout << "v1 = " << v1 << " v2 =" << v2 << "  v3 =" << v3 << std::endl;
#endif

    defcons_size = 12;
    prepareVectors (defcons_size);

    // compute the permutations of the 3 variables 
    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
   ind[0][0] = ibnd[0]; ind[0][1] = ibnd[1]; ind[0][2] = ibnd[2];
   ind[1][0] = ibnd[1]; ind[1][1] = ibnd[0]; ind[1][2] = ibnd[2];
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 2 && flagg == 0) {
      if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	 >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) 
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
      std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
#ifdef DEBUG
 std::cout << "v1 = " << v1 << " v2 =" << v2 << "  v3 =" << v3 << std::endl;
#endif

    double xL1 = cf*vlb[v1]; double xU1 = cf*vub[v1];
    double xL2 = vlb[v2]; double xU2 = vub[v2];
    double xL3 = vlb[v3]; double xU3 = vub[v3];

    for(int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     
      cutCoeff [ii][3] = 1.;
    }

    if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {
#ifdef DEBUG
    std::cout << " -- 5 if --" << std::endl;
#endif

      double theta1 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
      double theta2 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

      cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xL3; 
      cutCoeff [1][0] = -xL2*xL3; cutCoeff [1][1] = -xU1*xL3; cutCoeff [1][2] = -xU1*xL2;  bnd[1] = - 2.*xU1*xL2*xL3;
      cutCoeff [2][0] = -xU2*xU3; cutCoeff [2][1] = -xL1*xU3; cutCoeff [2][2] = -xL1*xL2;  bnd[2] = - xL1*xU2*xU3 - xL1*xL2*xU3;
      cutCoeff [3][0] = -xL2*xU3; cutCoeff [3][1] = -xU1*xU3; cutCoeff [3][2] = -xL1*xL2;  bnd[3] = - xU1*xL2*xU3 - xL1*xL2*xU3;
      cutCoeff [4][0] = -(theta1/(xU1-xL1)); cutCoeff [4][1] = -xU1*xU3; cutCoeff [4][2] = -xU1*xU2;  
      bnd[4] = (-(theta1*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;
      cutCoeff [5][0] = -xU2*xU3; cutCoeff [5][1] = -(theta2/(xU2-xL2)); cutCoeff [5][2] = -xU1*xU2;
      bnd[5] = (-(theta2*xL2)/(xU2-xL2)) - xU1*xU2*xL3 - xL1*xU2*xU3 + xL1*xL2*xL3;

    } else {
#ifdef DEBUG
    std::cout << " -- 5 else --" << std::endl;
#endif
      double theta1 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xU2*xU3;
      double theta2 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;

      cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xL3;
      cutCoeff [1][0] = -xL2*xL3; cutCoeff [1][1] = -xU1*xL3; cutCoeff [1][2] = -xU1*xL2;  bnd[1] = - 2.*xU1*xL2*xL3;
      cutCoeff [2][0] = -xL2*xU3; cutCoeff [2][1] = -xU1*xU3; cutCoeff [2][2] = -xU1*xU2;  bnd[2] = - xU1*xL2*xU3 - xU1*xU2*xU3; 
      cutCoeff [3][0] = -xU2*xU3; cutCoeff [3][1] = -xL1*xU3; cutCoeff [3][2] = -xU1*xU2;  bnd[3] = - xL1*xU2*xU3 - xU1*xU2*xU3; 
      cutCoeff [4][0] = -(theta1/(xL1-xU1)); cutCoeff [4][1] = -xL1*xU3; cutCoeff [4][2] = -xL1*xL2;  
      bnd[4] = (-(theta1*xU1)/(xL1-xU1)) - xL1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xL3;
      cutCoeff [5][0] = -xL2*xU3; cutCoeff [5][1] = -(theta2/(xL2-xU2)); cutCoeff [5][2] = -xL1*xL2;
      bnd[5] = (-(theta2*xU2)/(xL2-xU2)) - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;
    }

    double theta1c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
    double theta2c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

    cutCoeff [6][0] = -xL2*xL3; cutCoeff [6][1] = -xL1*xL3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - 2.*xL1*xL2*xL3;
    cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xU1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xL3;
    cutCoeff [8][0] = -xL2*xU3; cutCoeff [8][1] = -xL1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xL1*xL2*xU3 - xU1*xL2*xU3;
    cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xU1*xU3; cutCoeff [9][2] = -xU1*xL2;  bnd[9] = - xU1*xU2*xU3 - xU1*xL2*xU3;
    cutCoeff [10][0] = -(theta1c/(xL1-xU1)); cutCoeff [10][1] = -xL1*xU3; cutCoeff [10][2] = -xL1*xU2;  
    bnd[10] = (-(theta1c*xU1)/(xL1-xU1)) - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;
    cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -(theta2c/(xU2-xL2)); cutCoeff [11][2] = -xL1*xU2;
    bnd[11] = (-(theta2c*xL2)/(xU2-xL2)) - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;

  } // end if case 5

  /*----------------------------------------------------------------------------------------*/

  // case 6
  if(flag == 6) {
#ifdef DEBUG
    std::cout << " -- case 6 --" << std::endl;
#endif

    double theta = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
    double theta1 = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;

    defcons_size = 12;
    prepareVectors (defcons_size);

    for (int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     

      cutCoeff [ii][3] = 1.;
    }
    
    cutCoeff [0][0] = -xL2*xL3; cutCoeff [0][1] = -xU1*xL3; cutCoeff [0][2] = -xU1*xL2;  bnd[0] = - 2.*xU1*xL2*xL3;
    cutCoeff [1][0] = -xU2*xU3; cutCoeff [1][1] = -xL1*xL3; cutCoeff [1][2] = -xL1*xU2;  bnd[1] = - xL1*xU2*xU3 - xL1*xU2*xL3;
    cutCoeff [2][0] = -xU2*xU3; cutCoeff [2][1] = -xL1*xU3; cutCoeff [2][2] = -xL1*xL2;  bnd[2] = - xL1*xU2*xU3 - xL1*xL2*xU3;
    cutCoeff [3][0] = -xL2*xU3; cutCoeff [3][1] = -xU1*xU3; cutCoeff [3][2] = -xL1*xL2;  bnd[3] = - xU1*xL2*xU3 - xL1*xL2*xU3;
    cutCoeff [4][0] = -xU2*xL3; cutCoeff [4][1] = -xL1*xL3; cutCoeff [4][2] = -xU1*xU2;  bnd[4] = - xL1*xU2*xL3 - xU1*xU2*xL3;
    cutCoeff [5][0] = -(theta/(xU1-xL1)); cutCoeff [5][1] = -xU1*xU3; cutCoeff [5][2] = -xU1*xU2;
    bnd[5] = (-(theta*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;

    cutCoeff [6][0] = -xL2*xL3; cutCoeff [6][1] = -xL1*xL3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - 2.*xL1*xL2*xL3;
    cutCoeff [7][0] = -xL2*xU3; cutCoeff [7][1] = -xL1*xU3; cutCoeff [7][2] = -xU1*xL2;  bnd[7] = - xL1*xL2*xU3 - xU1*xL2*xU3; 
    cutCoeff [8][0] = -xU2*xU3; cutCoeff [8][1] = -xU1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xU1*xU2*xU3 - xU1*xL2*xU3; 
    cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xU1*xL3; cutCoeff [9][2] = -xU1*xU2;  bnd[9] = - xU1*xU2*xU3 - xU1*xU2*xL3;
    cutCoeff [10][0] = -xU2*xL3; cutCoeff [10][1] = -xU1*xL3; cutCoeff [10][2] = -xL1*xU2;  bnd[10] = - xL1*xU2*xL3 - xU1*xU2*xL3; 
    cutCoeff [11][0] = -(theta1/(xL1-xU1)); cutCoeff [11][1] = -xL1*xU3; cutCoeff [11][2] = -xL1*xU2;
    bnd[11] = (-(theta1*xU1)/(xL1-xU1)) - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;
  } // end if case 6
    
  /*----------------------------------------------------------------------------------------*/

  // case 7
  if(flag == 7) {
#ifdef DEBUG
    std::cout << " -- case 7 --" << std::endl;
#endif

    defcons_size = 12;
    prepareVectors (defcons_size);

    if ((vlb[v1]<=EPSILONT && vlb[v2]<=EPSILONT && vlb[v3]<=EPSILONT) ||
        (vlb[v1]==vlb[v2] && vlb[v1]==vlb[v3] && vub[v1]==vub[v2] && vub[v1]==vub[v3])) {
#ifdef DEBUG
      std::cout << " -- epsilonT --" << std::endl; 
#endif
    } else {
    // compute the 6 permutations of the 3 variables 
    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
    permutation3(ind,ibnd);
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 6 && flagg == 0) {
      if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	 <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
	 vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	 <= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) 
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
      std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

    double xL1 = cf*vlb[v1]; double xU1 = cf*vub[v1];
    double xL2 = vlb[v2]; double xU2 = vub[v2];
    double xL3 = vlb[v3]; double xU3 = vub[v3];

    //if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
    // <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] &&
    // vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
    // <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {

    double theta1 = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
    double theta2 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

    for(int ii = 0; ii < defcons_size; ii++) {
      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     

      cutCoeff [ii][3] = 1.;
    }
    
    cutCoeff [0][0] = -xL2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xL2;  bnd[0] = - 2.*xL1*xL2*xL3; 
    cutCoeff [1][0] = -xU2*xU3; cutCoeff [1][1] = -xU1*xU3; cutCoeff [1][2] = -xU1*xU2;  bnd[1] = - 2.*xU1*xU2*xU3;
    cutCoeff [2][0] = -xL2*xU3; cutCoeff [2][1] = -xL1*xU3; cutCoeff [2][2] = -xU1*xL2;  bnd[2] = - xL1*xL2*xU3 - xU1*xL2*xU3;
    cutCoeff [3][0] = -xU2*xL3; cutCoeff [3][1] = -xU1*xL3; cutCoeff [3][2] = -xL1*xU2;  bnd[3] = - xU1*xU2*xL3 - xL1*xU2*xL3;
    cutCoeff [4][0] = -(theta1/(xU1-xL1)); cutCoeff [4][1] = -xU1*xL3; cutCoeff [4][2] = -xU1*xL2;  
    bnd[4] = (-(theta1*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;
    cutCoeff [5][0] = -(theta2/(xL1-xU1)); cutCoeff [5][1] = -xL1*xU3; cutCoeff [5][2] = -xL1*xU2;
    bnd[5] = (-(theta2*xU1)/(xL1-xU1)) - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3; 
    //}
    cutCoeff [6][0] = -xL2*xL3; cutCoeff [6][1] = -xU1*xL3; cutCoeff [6][2] = -xU1*xU2;  bnd[6] = - xU1*xU2*xL3 - xU1*xL2*xL3; 
    cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xL1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - xU1*xU2*xL3 - xL1*xU2*xL3;
    cutCoeff [8][0] = -xL2*xL3; cutCoeff [8][1] = -xU1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xU1*xL2*xU3 - xU1*xL2*xL3; 
    cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xL1*xL3; cutCoeff [9][2] = -xL1*xU2;  bnd[9] = - xL1*xU2*xU3 - xL1*xU2*xL3;
    cutCoeff [10][0] = -xL2*xU3; cutCoeff [10][1] = -xU1*xU3; cutCoeff [10][2] = -xL1*xL2;  bnd[10] = - xU1*xL2*xU3 - xL1*xL2*xU3;
    cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -xL1*xU3; cutCoeff [11][2] = -xL1*xL2;  bnd[11] = - xL1*xU2*xU3 - xL1*xL2*xU3; 
  } // end if case 7

  /*----------------------------------------------------------------------------------------*/
      
  // case 8
  if(flag == 8) {
#ifdef DEBUG
    std::cout << " -- case 8 --" << std::endl;
#endif

    defcons_size = 12;
    prepareVectors (defcons_size);

    // compute the 6 permutations of the 3 variables 
    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
    permutation3(ind,ibnd);
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 6 && flagg == 0) {
      if(vub[ind[i][2]] <=0  &&
	 ((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	   >= vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] &&
	   vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	   >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) ||
	  (vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	   >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]])) ) 
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
      std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

    double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
    double xL2(vlb[v2]); double xU2(vub[v2]);
    double xL3(vlb[v3]); double xU3(vub[v3]);

    for(int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     
      cutCoeff [ii][3] = 1.;
    }
    
    // compute the 6 permutations of the 3 variables 
    //if(vub[v3]<=0) {
    cutCoeff [0][0] = -xU2*xL3; cutCoeff [0][1] = -xL1*xL3; cutCoeff [0][2] = -xL1*xL2;  bnd[0] = - xL1*xU2*xL3 - xL1*xL2*xL3; 
    cutCoeff [1][0] = -xU2*xL3; cutCoeff [1][1] = -xL1*xU3; cutCoeff [1][2] = -xL1*xU2;  bnd[1] = - xL1*xU2*xL3 - xL1*xU2*xU3; 
    cutCoeff [2][0] = -xL2*xU3; cutCoeff [2][1] = -xU1*xL3; cutCoeff [2][2] = -xU1*xL2;  bnd[2] = - xU1*xL2*xU3 - xU1*xL2*xL3;
    cutCoeff [3][0] = -xL2*xU3; cutCoeff [3][1] = -xU1*xU3; cutCoeff [3][2] = -xU1*xU2;  bnd[3] = - xU1*xL2*xU3 - xU1*xU2*xU3;
    cutCoeff [4][0] = -xL2*xL3; cutCoeff [4][1] = -xU1*xL3; cutCoeff [4][2] = -xL1*xL2;  bnd[4] = - xU1*xL2*xL3 - xL1*xL2*xL3;
    cutCoeff [5][0] = -xU2*xU3; cutCoeff [5][1] = -xL1*xU3; cutCoeff [5][2] = -xU1*xU2;  bnd[5] = - xU1*xU2*xU3 - xL1*xU2*xU3;

    if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       >= vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3] &&
       vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

      double theta1c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
      double theta2c = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

      cutCoeff [6][0] = -xL2*xU3; cutCoeff [6][1] = -xL1*xU3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - 2.*xL1*xL2*xU3;
      cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xU1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xL3;
      cutCoeff [8][0] = -xL2*xL3; cutCoeff [8][1] = -xU1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xU1*xL2*xU3 - xU1*xL2*xL3;
      cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xL1*xL3; cutCoeff [9][2] = -xL1*xU2;  bnd[9] = - xL1*xU2*xU3 - xL1*xU2*xL3;
      cutCoeff [10][0] = -xL2*xL3; cutCoeff [10][1] = -xL1*xL3; cutCoeff [10][2] = -(theta1c/(xL3-xU3));  
      bnd[10] = (-(theta1c*xU3)/(xL3-xU3)) - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;
      cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -xU1*xU3; cutCoeff [11][2] = -(theta2c/(xU3-xL3));
      bnd[11] = (-(theta2c*xL3)/(xU3-xL3)) - xU1*xL2*xU3 - xL1*xU2*xU3 + xL1*xL2*xL3;

    } else if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
	      >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

      double theta1c = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
      double theta2c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

      cutCoeff [6][0] = -xL2*xU3; cutCoeff [6][1] = -xL1*xU3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - 2.*xL1*xL2*xU3;
      cutCoeff [7][0] = -xU2*xL3; cutCoeff [7][1] = -xU1*xL3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xL3;
      cutCoeff [8][0] = -xL2*xL3; cutCoeff [8][1] = -xL1*xL3; cutCoeff [8][2] = -xL1*xU2;  bnd[8] = - xL1*xL2*xL3 - xL1*xU2*xL3;
      cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xU1*xU3; cutCoeff [9][2] = -xU1*xL2;  bnd[9] = - xU1*xL2*xU3 - xU1*xU2*xU3;
      cutCoeff [10][0] = -xL2*xL3; cutCoeff [10][1] = -(theta1c/(xL2-xU2)); cutCoeff [10][2] = -xU1*xL2;  
      bnd[10] = (-(theta1c*xU2)/(xL2-xU2)) - xL1*xL2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3; 
      cutCoeff [11][0] = -xU2*xU3; cutCoeff [11][1] = -(theta2c/(xU2-xL2)); cutCoeff [11][2] = -xL1*xU2;
      bnd[11] = (-(theta2c*xL2)/(xU2-xL2)) - xL1*xU2*xL3 - xU1*xU2*xU3 + xU1*xL2*xL3;
    }
    //} 
  } // end if case 8

  /*----------------------------------------------------------------------------------------*/

  // case 9
  if(flag == 9) {
#ifdef DEBUG
    std::cout << " -- case 9 --" << std::endl;
#endif

    defcons_size = 12;
    prepareVectors (defcons_size);

    // compute the 6 permutations of the 3 variables 
    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
    permutation3(ind,ibnd);
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 6 && flagg == 0) {
      if(vlb[ind[i][0]] >=0  &&
	 ((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	   <= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
	   vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	   <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) || 
	  (vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
	   <= vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] &&
	   vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
	   <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) ))
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
      std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

    double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
    double xL2(vlb[v2]); double xU2(vub[v2]);
    double xL3(vlb[v3]); double xU3(vub[v3]);

    for (int ii = 0; ii < defcons_size; ii++) {

      cutIndices [ii][0] = v1; 
      cutIndices [ii][1] = v2; 
      cutIndices [ii][2] = v3; 
      cutIndices [ii][3] = v4;     
      cutCoeff [ii][3] = 1.;
    }
    
    //if(vlb[v1]>=0) {
    if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] &&
       vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
       <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

      double theta1 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
      double theta2 = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xU2*xL3;

      cutCoeff [0][0] = -xU2*xU3; cutCoeff [0][1] = -xL1*xU3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xU3;
      cutCoeff [1][0] = -xL2*xL3; cutCoeff [1][1] = -xU1*xL3; cutCoeff [1][2] = -xU1*xL2;  bnd[1] = - 2.*xU1*xL2*xL3;
      cutCoeff [2][0] = -xL2*xU3; cutCoeff [2][1] = -xU1*xU3; cutCoeff [2][2] = -xL1*xL2;  bnd[2] = - xU1*xL2*xU3 - xL1*xL2*xU3;
      cutCoeff [3][0] = -xU2*xL3; cutCoeff [3][1] = -xL1*xL3; cutCoeff [3][2] = -xU1*xU2;  bnd[3] = - xU1*xU2*xL3 - xL1*xU2*xL3;
      cutCoeff [4][0] = -(theta1/(xL1-xU1)); cutCoeff [4][1] = -xL1*xL3; cutCoeff [4][2] = -xL1*xL2;  
      bnd[4] = (-(theta1*xU1)/(xL1-xU1)) - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xU2*xU3;
      cutCoeff [5][0] = -(theta2/(xU1-xL1)); cutCoeff [5][1] = -xU1*xU3; cutCoeff [5][2] = -xU1*xU2;
      bnd[5] = (-(theta2*xL1)/(xU1-xL1)) - xU1*xL2*xU3 - xU1*xU2*xL3 + xL1*xL2*xL3;

    } else if(vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
	      <= vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3] &&
	      vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
	      <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

      double theta1 = xU1*xU2*xU3 - xL1*xL2*xU3 - xU1*xU2*xL3 + xL1*xU2*xL3;
      double theta2 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;

      cutCoeff [0][0] = -xU2*xU3; cutCoeff [0][1] = -xL1*xU3; cutCoeff [0][2] = -xL1*xU2;  bnd[0] = - 2.*xL1*xU2*xU3;
      cutCoeff [1][0] = -xL2*xL3; cutCoeff [1][1] = -xU1*xL3; cutCoeff [1][2] = -xU1*xL2;  bnd[1] = - 2.*xU1*xL2*xL3;
      cutCoeff [2][0] = -xL2*xU3; cutCoeff [2][1] = -xU1*xU3; cutCoeff [2][2] = -xU1*xU2;  bnd[2] = - xU1*xL2*xU3 - xU1*xU2*xU3;
      cutCoeff [3][0] = -xU2*xL3; cutCoeff [3][1] = -xL1*xL3; cutCoeff [3][2] = -xL1*xL2;  bnd[3] = - xL1*xL2*xL3 - xL1*xU2*xL3;
      cutCoeff [4][0] = -xU2*xL3; cutCoeff [4][1] = -(theta1/(xU2-xL2)); cutCoeff [4][2] = -xU1*xU2;  
      bnd[4] = (-(theta1*xL2)/(xU2-xL2)) - xU1*xU2*xU3 - xL1*xU2*xL3 + xL1*xL2*xU3;
      cutCoeff [5][0] = -xL2*xU3; cutCoeff [5][1] = -(theta2/(xL2-xU2)); cutCoeff [5][2] = -xL1*xL2;
      bnd[5] = (-(theta2*xU2)/(xL2-xU2)) - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;
    }

    cutCoeff [6][0] = -xL2*xL3; cutCoeff [6][1] = -xL1*xU3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - xL1*xL2*xU3 - xL1*xL2*xL3;
    cutCoeff [7][0] = -xL2*xL3; cutCoeff [7][1] = -xL1*xL3; cutCoeff [7][2] = -xL1*xU2;  bnd[7] = - xL1*xU2*xL3 - xL1*xL2*xL3;
    cutCoeff [8][0] = -xL2*xU3; cutCoeff [8][1] = -xL1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xL1*xL2*xU3 - xU1*xL2*xU3;
    cutCoeff [9][0] = -xU2*xU3; cutCoeff [9][1] = -xU1*xU3; cutCoeff [9][2] = -xU1*xL2;  bnd[9] = - xU1*xU2*xU3 - xU1*xL2*xU3;
    cutCoeff [10][0] = -xU2*xU3; cutCoeff [10][1] = -xU1*xL3; cutCoeff [10][2] = -xU1*xU2;  bnd[10] = - xU1*xU2*xU3 - xU1*xU2*xL3;
    cutCoeff [11][0] = -xU2*xL3; cutCoeff [11][1] = -xU1*xL3; cutCoeff [11][2] = -xL1*xU2;  bnd[11] = - xL1*xU2*xL3 - xU1*xU2*xL3;
    //}
  } // end if case 9

  /*----------------------------------------------------------------------------------------*/

  // case 10
  if(flag == 10) {
#ifdef DEBUG
    std::cout << " -- case 10 --" << std::endl;
#endif

    defcons_size = 12;
    prepareVectors (defcons_size);

    ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
    permutation3(ind,ibnd);
    int i, flagg=0, idx=0;
    i = 0;
    while(i < 6 && flagg == 0) {
      if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	 >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
	 vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
	 >= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) 
	{
	  idx = i;   // store the index of the permutation satisfying the condition
	  flagg = 1;  // condition is satisfied
	}
      i++; 
    }
    if (flagg==0) {
       std::cout << "ERROR!!!" << std::endl; exit(0);
    }
    v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

      double xL1(cf*vlb[v1]); double xU1(cf*vub[v1]);
      double xL2(vlb[v2]); double xU2(vub[v2]);
      double xL3(vlb[v3]); double xU3(vub[v3]);

      for(int ii = 0; ii < defcons_size; ii++) {
 
        cutIndices [ii][0] = v1; 
        cutIndices [ii][1] = v2; 
        cutIndices [ii][2] = v3; 
        cutIndices [ii][3] = v4;     
        cutCoeff [ii][3] = 1.;
      }
    
    // compute the 6 permutations of the 3 variables 
      cutCoeff [0][0] = -xL2*xL3; cutCoeff [0][1] = -xU1*xL3; cutCoeff [0][2] = -xU1*xU2;  bnd[0] = - xU1*xU2*xL3 - xU1*xL2*xL3;
      cutCoeff [1][0] = -xU2*xU3; cutCoeff [1][1] = -xL1*xL3; cutCoeff [1][2] = -xL1*xU2;  bnd[1] = - xL1*xU2*xU3 - xL1*xU2*xL3;
      cutCoeff [2][0] = -xU2*xL3; cutCoeff [2][1] = -xL1*xL3; cutCoeff [2][2] = -xU1*xU2;  bnd[2] = - xU1*xU2*xL3 - xL1*xU2*xL3;
      cutCoeff [3][0] = -xU2*xU3; cutCoeff [3][1] = -xL1*xU3; cutCoeff [3][2] = -xL1*xL2;  bnd[3] = - xL1*xU2*xU3 - xL1*xL2*xU3;
      cutCoeff [4][0] = -xL2*xU3; cutCoeff [4][1] = -xU1*xU3; cutCoeff [4][2] = -xL1*xL2;  bnd[4] = - xU1*xL2*xU3 - xL1*xL2*xU3;
      cutCoeff [5][0] = -xL2*xL3; cutCoeff [5][1] = -xU1*xU3; cutCoeff [5][2] = -xU1*xL2;  bnd[5] = - xU1*xL2*xU3 - xU1*xL2*xL3;

      //if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
      // >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] &&
      // vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
      // >= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {

      double theta1c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
      double theta2c = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;

      cutCoeff [6][0] = -xL2*xL3; cutCoeff [6][1] = -xL1*xL3; cutCoeff [6][2] = -xL1*xL2;  bnd[6] = - 2.*xL1*xL2*xL3;
      cutCoeff [7][0] = -xU2*xU3; cutCoeff [7][1] = -xU1*xU3; cutCoeff [7][2] = -xU1*xU2;  bnd[7] = - 2.*xU1*xU2*xU3;
      cutCoeff [8][0] = -xL2*xU3; cutCoeff [8][1] = -xL1*xU3; cutCoeff [8][2] = -xU1*xL2;  bnd[8] = - xL1*xL2*xU3 - xU1*xL2*xU3;
      cutCoeff [9][0] = -xU2*xL3; cutCoeff [9][1] = -xU1*xL3; cutCoeff [9][2] = -xL1*xU2;  bnd[9] = - xL1*xU2*xL3 - xU1*xU2*xL3;
      cutCoeff [10][0] = -(theta1c/(xL1-xU1)); cutCoeff [10][1] = -xL1*xU3; cutCoeff [10][2] = -xL1*xU2;  
      bnd[10] = (-(theta1c*xU1)/(xL1-xU1)) - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;
      cutCoeff [11][0] = -(theta2c/(xU1-xL1)); cutCoeff [11][1] = -xU1*xL3; cutCoeff [11][2] = -xU1*xL2;
      bnd[11] = (-(theta2c*xL1)/(xU1-xL1)) - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;

  } // end if case 10

  /*----------------------------------------------------------------------------------------*/

  // get lower and upper bound on the constraints to be added to the problem
  for (int ii = 0; ii < defcons_size; ii++) {

    //printf ("flag = %d\n", flag);

    if (ii < defcons_size/2) {

      if (cf > 0) {cutLb[ii] =  bnd[ii];           cutUb[ii] = COUENNE_INFINITY;}
      if (cf < 0) {cutLb[ii] = -COUENNE_INFINITY;  cutUb[ii] = bnd[ii];}

    } else {

      if (cf > 0) {cutLb[ii] = -COUENNE_INFINITY;  cutUb[ii] = bnd[ii];}
      if (cf < 0) {cutLb[ii] =  bnd[ii];           cutUb[ii] = COUENNE_INFINITY;}
    }
#ifdef DEBUG
std::cout << ii << ") cutLb =" << cutLb[ii] << " " << "cutUb = " << cutUb[ii] << std::endl;
#endif
  }

  for (int i=0; i<6; i++)
    delete [] ind [i];

  delete [] ibnd;
  delete [] bnd;
  delete [] ind;
}


// generate cuts for trilinear expressions
void exprTrilinear::generateCuts (expression *w, 
				  OsiCuts &cs, const CouenneCutGenerator *cg,
				  t_chg_bounds *chg, int wind, 
				  CouNumber lbw, CouNumber ubw) {

  expression **args = w -> Image () -> ArgList ();

  int *varInd = new int [4];

  for (int i=0; i<3; i++)
    varInd [i] = args [i] -> Index (); 

  varInd [3] = w -> Index ();	

  std::vector <std::vector <int> >    cutIndices;
  std::vector <std::vector <double> > cutCoeff;
  std::vector <double>                cutLb, cutUb;

  TriLinCuts (cg -> Problem () -> Lb (),
	      cg -> Problem () -> Ub (),
	      varInd,
	      cutIndices, cutCoeff,
	      cutLb, cutUb);

  // sanity check on returned vectors
  assert (cutIndices.size () == cutCoeff.size () && 
	  cutIndices.size () == cutLb.size    () && 
	  cutIndices.size () == cutUb.size    ());

  //printf ("trilinear cuts:\n");

  for (int i = (int) cutIndices.size (); i--;) {

    int 
       size = (int) cutIndices [i].size (),
      *ind  = new int [size];

    double *coe = new double [size];

    std::copy (cutIndices [i].begin (), cutIndices [i].end (), ind);
    std::copy (cutCoeff   [i].begin (), cutCoeff   [i].end (), coe);

    OsiRowCut cut (cutLb [i], cutUb [i], 4, 4, ind, coe);
    //cut.print ();

    delete [] ind;
    delete [] coe;

    if (cg -> Problem () -> bestSol ()) {

      // check validity of cuts by verifying they don't cut the
      // optimum in the node

      double *sol = cg -> Problem () -> bestSol ();
      const double
	*lb = cg -> Problem () -> Lb (), 
	*ub = cg -> Problem () -> Ub ();

      int nVars = cg -> Problem () -> nVars ();

      bool optIn = true;

      for (int i=0; i< nVars; i++)
	if ((sol [i] < lb [i] - COUENNE_EPS) ||
	    (sol [i] > ub [i] + COUENNE_EPS)) {
	  optIn = false;
	  break;
	}

      if (optIn) {

	for (unsigned int i=0; i<cutIndices.size (); i++) {

	  double chs = 0.;

	  for (unsigned int j=0; j<cutIndices[i].size(); j++)
	    chs += cutCoeff [i] [j] * sol [cutIndices [i] [j]];

	  if ((chs < cutLb [i] - COUENNE_EPS) ||
	      (chs > cutUb [i] + COUENNE_EPS)) {

	    printf ("cut %d violates optimum:\n", i);

	    if (cutLb [i] > -COUENNE_INFINITY) printf ("%g <= ", cutLb [i]);
	    for (unsigned int j=0; j<cutIndices[i].size(); j++) printf ("%+g x%d ", cutCoeff [i] [j],       cutIndices [i] [j]);  printf ("\n = ");
	    for (unsigned int j=0; j<cutIndices[i].size(); j++) printf ("%+g *%g ", cutCoeff [i] [j],  sol [cutIndices [i] [j]]); printf ("\n = ");
	    for (unsigned int j=0; j<cutIndices[i].size(); j++) printf ("%+g ",     cutCoeff [i] [j] * sol [cutIndices [i] [j]]); printf (" = %g", chs);
	    if (cutUb [i] <  COUENNE_INFINITY) printf (" <= %g", cutUb [i]);
	    printf ("\n");

	  }
	}
      }
    }

    cs.insert (cut);
  }

  delete [] varInd;
}
