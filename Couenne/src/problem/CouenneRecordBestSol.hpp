// (C) Copyright Francois Margot and Carnegie Mellon University 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Francois Margot, Tepper School of Business, Carnegie Mellon University,
//
// Date : 3/31/2011

#ifndef RECBESTSOL_H
#define RECBESTSOL_H

#include "CouenneTypes.hpp"


namespace Couenne {

// class to record best found feasible solution
class CouenneRecordBestSol {

public:

  // size of initial domain
  int cardInitDom;
  // vector of length cardInitDom indicating if a variable is integer or not;
  bool *initIsInt;
  // vector of indices of integer variables
  std::vector<int> listInt;
  // copy of initial domain lower bounds
  CouNumber *initDomLb;
  // copy of initial domain upper bounds
  CouNumber *initDomUb;

  // true if a solution value has been recorded, false otherwise
  bool hasSol;
  // size of vector sol
  int cardSol;
  // if not NULL, recorded solution
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
  CouenneRecordBestSol();

  /// Copy constructor
  CouenneRecordBestSol(const CouenneRecordBestSol &other);

  /// Destructor
  ~CouenneRecordBestSol();

  inline int getCardInitDom() const {return cardInitDom;};
  inline bool *getInitIsInt() const {return initIsInt;};
  inline std::vector<int> getListInt() const {return listInt;};
  // set both initIsInt and listInt from the given vector givenIsInt
  void setInitIsInt(const bool *givenIsInt, const int givenCard);
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

  // record givenSol as best solution if givenVal is smaller
  // than val (or if no solution was recorded previously)
  void update(const double *givenSol, const int givenCard, 
	      const double givenVal, const double givenMaxViol);

  // use modSol, modSolVal, modSolMaxViol for updating only if
  // modSolVal is smaller than val (or if no solution was recorded previously)
  void update();

  // compare given two solutions and set sol, solVal, and maxViol to
  // the best of the two with finite value (< 1e49); return -1 if both have
  // infinite value, return 0 if solA is saved, return 1 if solB is saved  
  int compareAndSave(const double *solA, const double solAVal,
		     const double solAMaxViol, 
		     const bool solAIsFeas,
		     const double *solB, const double solBVal,
		     const double solBMaxViol, 
		     const bool solBIsFeas,
		     const int cardSol,
		     const double precision);

  inline int getCardModSol() const {return cardModSol;}
  double *getModSol(const int expectedCard);
  inline double getModSolVal() const {return modSolVal;}
  inline double getModSolMaxViol() const {return modSolMaxViol;}

  // set modSol, modSolVal, and modSolMaxViol to given values; if
  // givenModSol == NULL, only the other two are set
  void setModSol(const double *givenModSol, const int givenModCard, 
		 const double givenModVal, const double givenModMaxViol);

  // print sol, solVal, and maxViol
  void printSol(FILE *fsol) const;
};

}

#endif
