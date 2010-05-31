/* $Id$
 *
 * Name:    CouenneFeasPump.cpp
 * Authors: Pietro Belotti, Lehigh University
 *          Timo Berthold, ZIB Berlin
 * Purpose: Implement the Feasibility Pump heuristic class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneFeasPump.hpp"
#include "BonCouenneInterface.hpp"
#include "CouenneMINLPInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "BonAuxInfos.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

using namespace Couenne;

// Constructor ////////////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump():

  CbcHeuristic (),
  nlp_         (NULL),
  hasCloned_   (false),
  //maxNlpInf_   (maxNlpInf_0),
  numberSolvePerLevel_(-1),
  problem_     (NULL){
  setHeuristicName ("Couenne Feasibility Pump");
}

// Constructor ////////////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (CbcModel &model, 
				  CouenneMINLPInterface &nlp, 
				  bool cloneNlp, 
				  CouenneProblem * couenne):

  CbcHeuristic (model), 
  nlp_         (&nlp), 
  hasCloned_   (cloneNlp), 
  //maxNlpInf_   (maxNlpInf_0),
  numberSolvePerLevel_(-1),
  problem_     (couenne) {

  setHeuristicName("Couenne Feasibility Pump");

  if (cloneNlp)
    nlp_ = dynamic_cast <CouenneMINLPInterface *> (nlp. clone());
  }
  
// Copy constructor ///////////////////////////////////////////// 
CouenneFeasPump::CouenneFeasPump (const CouenneFeasPump &other):

  CbcHeuristic (other), 
  nlp_         (other.nlp_), 
  hasCloned_   (other.hasCloned_),
  //maxNlpInf_   (other.maxNlpInf_),
  numberSolvePerLevel_ (other.numberSolvePerLevel_),
  problem_     (other.problem_) {

  if (hasCloned_ && nlp_ != NULL)
    nlp_ = dynamic_cast <CouenneMINLPInterface *> (other.nlp_ -> clone ());
}

// Clone //////////////////////////////////////////////////////// 
CbcHeuristic *CouenneFeasPump::clone () const 
{return new CouenneFeasPump (*this);}

// Assignment operator ////////////////////////////////////////// 
CouenneFeasPump &CouenneFeasPump::operator= (const CouenneFeasPump & rhs) {

  if (this != &rhs) {

    CbcHeuristic::operator= (rhs);

    if (hasCloned_ && nlp_)
      delete nlp_;
      
    hasCloned_ = rhs.hasCloned_;

    if (nlp_ != NULL){

      if (hasCloned_) nlp_ = dynamic_cast <CouenneMINLPInterface *> (rhs.nlp_ -> clone());
      else            nlp_ = rhs.nlp_;
    }
  }

  //maxNlpInf_           = rhs.maxNlpInf_;
  numberSolvePerLevel_ = rhs.numberSolvePerLevel_;
  problem_             = rhs.problem_;

  return *this;
}

// Destructor /////////////////////////////////////////////////// 
CouenneFeasPump::~CouenneFeasPump () {
  if (hasCloned_)
    delete nlp_;
  nlp_ = NULL;
}

// Pass pointer to NLP solver /////////////////////////////////// 
void CouenneFeasPump::setNlp (CouenneInterface &nlp, bool cloneNlp) {

  // FIXME: need a CouenneMINLPInterface, not a CouenneInterface

  // if (hasCloned_ && nlp_ != NULL)
  //   delete nlp_;

  // hasCloned_ = cloneNlp;

  // if (cloneNlp) nlp_ = dynamic_cast <CouenneMINLPInterface *> (nlp. clone ());
  // else          nlp_ = &nlp;
}

// Pass pointer to problem description ////////////////////////// 
void CouenneFeasPump::setCouenneProblem (CouenneProblem * couenne) {
  problem_ = couenne;
}

/// set new expression as the NLP objective function using
/// argument as point to minimize distance from. Return new
/// objective function
expression *CouenneFeasPump::updateNLPObj (double *) {

  return NULL;
}

/// admits a (possibly fractional) solution and fixes the integer
/// components in the nonlinear problem for later re-solve
void CouenneFeasPump::fixIntVariables (double *sol) {

}
