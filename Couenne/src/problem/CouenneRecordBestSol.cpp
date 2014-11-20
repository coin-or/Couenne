// (C) Copyright Francois Margot and Carnegie Mellon University 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Francois Margot, Tepper School of Business, Carnegie Mellon University,
//
// Date : 3/31/2011

#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneRecordBestSol.hpp"

using namespace Couenne;

//#define TRACE

/*************************************************************/
/** Default constructor. */
CouenneRecordBestSol::CouenneRecordBestSol() {

  cardInitDom = -1;
  initIsInt = NULL;
  initDomLb = NULL;
  initDomUb = NULL;

  hasSol = false;
  cardSol = -1;
  sol = NULL;
  val = -1;
  maxViol = -1;

  cardModSol = -1;
  modSol = NULL;
  modSolVal = -1;
  modSolMaxViol = -1;
}

/*************************************************************/
// copy constructor
CouenneRecordBestSol::CouenneRecordBestSol(const CouenneRecordBestSol &other) {

  cardInitDom = other.cardInitDom;
  if(cardInitDom > -1) {
    initIsInt = new bool[other.cardInitDom];
    initDomLb = new CouNumber[other.cardInitDom];
    initDomUb = new CouNumber[other.cardInitDom];

    CoinCopyN(other.initIsInt, cardInitDom, initIsInt);
    CoinCopyN(other.initDomLb, cardInitDom, initDomLb);
    CoinCopyN(other.initDomUb, cardInitDom, initDomUb);
  }
  else {
    initIsInt = NULL;
    initDomLb = NULL;
    initDomUb = NULL;
  }

  for(unsigned int i=0; i<other.listInt.size(); i++) {
    listInt.push_back(other.listInt[i]);
  }

  hasSol = other.hasSol;
  cardSol = other.cardSol;
  val = other.val;
  maxViol = other.maxViol;

  if(other.sol != NULL) {
    sol = new double[other.cardSol];
    CoinCopyN(other.sol, cardSol, sol);
  }
  else {
    sol = NULL;
  }

  if (other.modSol != NULL) {
    modSol = new double[other.cardSol];
    CoinCopyN(other.modSol, cardSol, modSol);
  }
  else {
    modSol = NULL;
  }
  cardModSol = other.cardModSol;
  modSolVal = other.modSolVal;
  modSolMaxViol = other.modSolMaxViol;
} 

/*************************************************************/
/** Destructor. */
CouenneRecordBestSol::~CouenneRecordBestSol(){

  if(cardInitDom > -1) {
    delete[] initIsInt;
    delete[] initDomLb;
    delete[] initDomUb;
  }

  if(sol != NULL) {
    delete[] sol;
  }

  if(modSol != NULL) {
    delete[] modSol;
  }
}

/*****************************************************************************/
void CouenneRecordBestSol::setInitIsInt(const bool *givenIsInt,
					const int givenCard) {

  if(initIsInt == NULL) {
    if(cardInitDom == -1) {
      cardInitDom = givenCard;
    }
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitIsInt(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
    initIsInt = new bool[givenCard];
  }
  else {
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitIsInt(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
  }
  CoinCopyN(givenIsInt, givenCard, initIsInt);

  listInt.empty();
  for(int i=0; i<givenCard; i++) {
    if(initIsInt[i]) {
      listInt.push_back(i);
    }
  }
} /* setInitIsInt */

/*****************************************************************************/
void CouenneRecordBestSol::setInitDomLb(const CouNumber *givenLb, 
					const int givenCard) {
  if(initDomLb == NULL) {
    if(cardInitDom == -1) {
      cardInitDom = givenCard;
    }
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitDomLb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
    initDomLb = new CouNumber[givenCard];
  }
  else {
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitDomLb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
  }
  CoinCopyN(givenLb, givenCard, initDomLb);
} /* setInitDomLb */

/*****************************************************************************/
void CouenneRecordBestSol::setInitDomUb(const CouNumber *givenUb, 
					const int givenCard) {
  if(initDomUb == NULL) {
    if(cardInitDom == -1) {
      cardInitDom = givenCard;
    }
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitDomUb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
    initDomUb = new CouNumber[givenCard];
  }
  else {
    if(givenCard != cardInitDom) {
      printf("### ERROR: CouenneRecordBestSol::setInitDomUb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
  }
  CoinCopyN(givenUb, givenCard, initDomUb);
} /* setInitDomUb */

/*****************************************************************************/
void CouenneRecordBestSol::setHasSol(const bool givenHasSol) {
  hasSol = givenHasSol;
}

/*****************************************************************************/
void CouenneRecordBestSol::setCardSol(const int givenCard) {
  cardSol = givenCard;
}

/*****************************************************************************/
void CouenneRecordBestSol::setSol(const double *givenSol, const int givenCard,
                                  const double givenMaxViol) {
  if(sol == NULL) {
    cardSol = givenCard;
    sol = new double[givenCard];
    if(modSol == NULL) {
      modSol = new double[givenCard];
    }
  }
  else {
    if(givenCard != cardSol) {
      //printf("CouenneRecordBestSol::setSol(): ### ERROR: givenCard: %d  cardSol: %d", givenCard, cardSol);
      //exit(1);

	double *newSol = new double [givenCard];
	CoinCopyN (givenSol, givenCard, newSol);
	delete [] modSol;
	modSol = newSol;
	cardSol = givenCard;
    }
  }
  CoinCopyN(givenSol, givenCard, sol);
  maxViol = givenMaxViol;

#ifdef TRACE
  printf("CouenneRecordBestSol::setSol(): New solution set\n");
#endif

} /* setSol */

/*****************************************************************************/
void CouenneRecordBestSol::setVal(const double givenVal) {

#ifdef TRACE
  printf("CouenneRecordBestSol::setVal(): set to %10.6f\n", givenVal);
#endif

  val = givenVal;
  hasSol = true;
}

/*****************************************************************************/
void CouenneRecordBestSol::update(const double *givenSol, const int givenCard, 
			   const double givenVal, const double givenMaxViol) {
  if (!hasSol || (givenVal < val)) {
    setSol(givenSol, givenCard, givenMaxViol);
    setVal(givenVal);
  }
} /* update */

/*****************************************************************************/
void CouenneRecordBestSol::update() {
  if(modSol == NULL) {
    printf(" CouenneRecordBestSol::update(): ### ERROR: modSol == NULL\n");
    exit(1);
  }

  update(modSol, cardModSol, modSolVal, modSolMaxViol);
} /* update */

/*****************************************************************************/
int CouenneRecordBestSol::compareAndSave(const double *solA, const double solAVal, const double solAMaxViol, const bool solAIsFeas,
					 const double *solB, const double solBVal, const double solBMaxViol, const bool solBIsFeas,
					 const int cardSol,
					 const double precision) {
  int retval = -2;
  if(solBIsFeas) {
    if(solAIsFeas) {
      if(solAVal < solBVal - precision) {
	retval = 0;
      }
      else {
	retval = 1;
      }
    }
    else {
      retval = 1;
    }
  }
  else {
    if(solAIsFeas) {
      retval = 0;
    }
    else { // both solutions are infeasible; select the one with min viol.
      if(solAVal < 1e49) {
	if(solBVal < 1e49) {
	  if(solAMaxViol < solBMaxViol) {
	    retval = 0;
	  }
	  else {
	    retval = 1;
	  }
	}
	else {
	  retval = 0;
	}
      }
      else {
	if(solBVal < 1e49) {
	  retval = 1;
	}
	else {
	  retval = -1;
	}
      }
    }
  }
  
  switch (retval) {
    case 0: update(solA, cardSol, solAVal, solAMaxViol); break;
    case 1: update(solB, cardSol, solBVal, solBMaxViol); break;
    case -1: break;
    default: printf("CouenneRecordBestSol::compareAndSave(): ### ERROR: retval: %d\n",
		    retval); break;
  }

  return(retval);
} /* compareAndSave */ 

/*****************************************************************************/
double * CouenneRecordBestSol::getModSol(const int expectedCard) { 
  if(modSol == NULL) {
    cardModSol = expectedCard;
    modSol = new double[expectedCard];
  }
  else {
    if(expectedCard != cardModSol) {
      printf("CouenneRecordBestSol::getModSol(): ### ERROR: expectedCard: %d  cardModSol: %d", expectedCard, cardModSol);
      exit(1);
    }
  }
  return modSol;
} /* getModSol */

/*****************************************************************************/
void CouenneRecordBestSol::setModSol(const double *givenModSol, 
			      const int givenModCard, 
			      const double givenModVal, 
			      const double givenModMaxViol) {
  
  if(givenModSol != NULL) {
    if(modSol == NULL) {
      cardModSol = givenModCard;
      modSol = new double[givenModCard];
    }
    else {
      if(givenModCard != cardModSol) {
	// printf("CouenneRecordBestSol::setModSol(): ### ERROR: givenModCard: %d  cardModSol: %d", givenModCard, cardModSol);
	// exit(1);

	double *newModSol = new double [givenModCard];
	CoinCopyN (givenModSol, givenModCard, newModSol);
	delete [] modSol;
	modSol = newModSol;
	cardModSol = givenModCard;
      }
    }
    CoinCopyN(givenModSol, givenModCard, modSol);
  }
  modSolVal = givenModVal;
  modSolMaxViol = givenModMaxViol;
} /* setModSol */

/*****************************************************************************/
void CouenneRecordBestSol::printSol(FILE *fsol) const {

  if(sol != NULL) {
    fprintf(fsol, "%d\n", cardSol);
    for(int i=0; i<cardSol; i++) {
      fprintf(fsol, " %12.8f", sol[i]);
      if(i % 10 == 9) {
	fprintf(fsol, "\n");
      }
    }
    if(cardSol % 10 != 0) {
      fprintf(fsol, "\n");	
    }
    fprintf(fsol, "Value: %16.14g\n", val);
    fprintf(fsol, "Tolerance: %16.14g\n", maxViol);
  }
} /* printSol */ 
