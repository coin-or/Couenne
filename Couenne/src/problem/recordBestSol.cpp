// $Id$
//
// (C) Copyright Francois Margot and Carnegie Mellon University 2011
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Francois Margot, Tepper School of Business, Carnegie Mellon University,
//
// Date : 3/31/2011

#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>

#include "CouenneProblem.hpp"
#include "recordBestSol.hpp"

/*************************************************************/
  /** Default constructor. */
recordBestSol::recordBestSol() {

  cardInitDom = -1;
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
recordBestSol::recordBestSol(const recordBestSol &other) {

  cardInitDom = other.cardInitDom;
  if(cardInitDom > -1) {
    initDomLb = new CouNumber[other.cardInitDom];
    initDomUb = new CouNumber[other.cardInitDom];

    CoinCopyN(other.initDomLb, cardInitDom, initDomLb);
    CoinCopyN(other.initDomUb, cardInitDom, initDomUb);
  }
  else {
    initDomLb = NULL;
    initDomUb = NULL;
  }

  hasSol = other.hasSol;
  cardSol = other.cardSol;
  val = other.val;
  maxViol = other.maxViol;

  if(hasSol) {
    sol = new double[other.cardSol];
    CoinCopyN(other.sol, cardSol, sol);
  }
  else {
    sol = NULL;
  }

  if(modSol != NULL) {
    modSol = new double[other.cardSol];
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
recordBestSol::~recordBestSol(){

  delete[] initDomLb;
  delete[] initDomUb;

  if(sol != NULL) {
    delete[] sol;
  }

  if(modSol != NULL) {
    delete[] modSol;
  }
}

/*****************************************************************************/
void recordBestSol::setInitDomLb(const CouNumber *givenLb, 
				 const int givenCard) {
  if(initDomLb == NULL) {
    cardInitDom = givenCard;
    initDomLb = new CouNumber[givenCard];
  }
  else {
    if(givenCard != cardInitDom) {
      printf("### ERROR: recordBestSol::setInitDomLb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
  }
  CoinCopyN(givenLb, givenCard, initDomLb);
} /* setInitDomLb */

/*****************************************************************************/
void recordBestSol::setInitDomUb(const CouNumber *givenUb, 
				 const int givenCard) {
  if(initDomUb == NULL) {
    cardInitDom = givenCard;
    initDomUb = new CouNumber[givenCard];
  }
  else {
    if(givenCard != cardInitDom) {
      printf("### ERROR: recordBestSol::setInitDomUb(): cardInitDom: %d  givenCard: %d\n", cardInitDom, givenCard);
      exit(1);
    }
  }
  CoinCopyN(givenUb, givenCard, initDomUb);
} /* setInitDomLb */

/*****************************************************************************/
void recordBestSol::setHasSol(const bool givenHasSol) {
  hasSol = givenHasSol;
}

/*****************************************************************************/
void recordBestSol::setCardSol(const int givenCard) {
  cardSol = givenCard;
}

/*****************************************************************************/
void recordBestSol::setSol(const double *givenSol, const int givenCard,
			   const double givenMaxViol) {
  if(!hasSol) {
    hasSol = true;
    cardSol = givenCard;
    sol = new double[givenCard];
    if(modSol == NULL) {
      modSol = new double[givenCard];
    }
  }
  else {
    if(givenCard != cardSol) {
      printf("recordBestSol::setSol(): ### ERROR: givenCard: %d  cardSol: %d", givenCard, cardSol);
      exit(1);
    }
  }
  CoinCopyN(givenSol, givenCard, sol);
  maxViol = givenMaxViol;
} /* setSol */

/*****************************************************************************/
void recordBestSol::setVal(const double givenVal) {

#ifdef TRACE
  printf("recordBestSol::setVal(): set to %g\n", givenVal);
#endif

  val = givenVal;
}

/*****************************************************************************/
void recordBestSol::update(const double *givenSol, const int givenCard, 
			   const double givenVal, const double givenMaxViol) {
  if((!hasSol) || ((hasSol) && (givenVal < val))) {
    setSol(givenSol, givenCard, givenMaxViol);
    setVal(givenVal);
  }
} /* update */

/*****************************************************************************/
void recordBestSol::update() {
  if(modSol == NULL) {
    printf(" recordBestSol::update(): ### ERROR: modSol == NULL\n");
    exit(1);
  }
  update(modSol, cardModSol, modSolVal, modSolMaxViol);
} /* update */

/*****************************************************************************/
double * recordBestSol::getModSol(const int expectedCard) { 
  if(modSol == NULL) {
    cardModSol = expectedCard;
    modSol = new double[expectedCard];
  }
  else {
    if(expectedCard != cardModSol) {
      printf("recordBestSol::getModSol(): ### ERROR: expectedCard: %d  cardModSol: %d", expectedCard, cardModSol);
      exit(1);
    }
  }
  return modSol;
} /* getModSol */

/*****************************************************************************/
void recordBestSol::setModSol(const double *givenModSol, 
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
	printf("recordBestSol::setModSol(): ### ERROR: givenModCard: %d  cardModSol: %d", givenModCard, cardModSol);
	exit(1);
      }
    }
    CoinCopyN(givenModSol, givenModCard, modSol);
  }
  modSolVal = givenModVal;
  modSolMaxViol = givenModMaxViol;
} /* setModSol */
