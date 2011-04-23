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

#ifndef RECBESTSOL_H
#define RECBESTSOL_H

class CouenneProblem;

// class to record best found feasible solution
class recordBestSol {

public:

  // size of initial domain
  int cardInitDom;
  // copy of initial domain lower bounds
  CouNumber *initDomLb;
  // copy of initial domain upper bounds
  CouNumber *initDomUb;

  // true if a solution has been recorded, false otherwise
  bool hasSol;
  // size of vector sol
  int cardSol;
  // recorded solution
  double *sol;
  // recorded value
  double val;
  // recorded maximum violation of a bound, integrality or constraint
  double maxViol;

  // used by checkNLP2 and for update
  int cardModSol;
  double *modSol; 
  double modSolVal; 
  double modSolMaxViol; 

public:
  /// Constructor
  recordBestSol();

  /// Copy constructor
  recordBestSol(const recordBestSol &other);

  /// Destructor
  ~recordBestSol();

  inline int getCardInitDom() const {return cardInitDom;};
  inline CouNumber *getInitDomLb() const {return initDomLb;};
  void setInitDomLb(const CouNumber *givenLb, const int givenCard);
  inline CouNumber *getInitDomUb() const {return initDomUb;};
  void setInitDomUb(const CouNumber *givenUb, const int givenCard);

  void setHasSol(const bool givenHasSol);
  inline bool getHasSol() const {return hasSol;};
  void setSol(const double *givenSol, const int givenCard, 
	      const double givenMaxViol);
  inline int getCardSol() const {return cardSol;};
  void setCardSol(const int givenCard);
  inline double *getSol() const {return sol;};
  inline double getMaxViol() const {return maxViol;};
  void setVal(const double givenVal);
  inline double getVal() {return val;};

  // record givenSol has best solution if givenVal is smaller
  // than val (or if no solution was recorded previously)
  void update(const double *givenSol, const int givenCard, 
	      const double givenVal, const double givenMaxViol);

  // use modSol, modSolVal, modSolMaxViol for updating only if
  // modSolVal is smaller than val (or if no solution was recorded previously)
  void update();

  inline int getCardModSol() const {return cardModSol;};
  double *getModSol(const int expectedCard);
  inline double getModSolVal() const {return modSolVal;};
  inline double getModSolMaxViol() const {return modSolMaxViol;};

  // set modSol, modSolVal, and modSolMaxViol to given values; if
  // givenModSol == NULL, only the other two are set
  void setModSol(const double *givenModSol, int givenModCard, 
		 double givenModVal, double givenModMaxViol);

};

#endif
