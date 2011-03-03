/* $Id$
 *
 * Name:    CouenneSolverInterface.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 *          Andreas Waechter, IBM Corp.
 * Purpose: Implementation of the OsiSolverInterface::resolve () method 
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

/// constructor
template <class T> 
CouenneSolverInterface<T>::CouenneSolverInterface 
(CouenneCutGenerator *cg /*= NULL*/):

  T (),
  cutgen_ (cg),
  knowInfeasible_(false),
  knowOptimal_(false),
  knowDualInfeasible_(false) {}
//  doingResolve_ (true) {}


/// copy constructor
template <class T> 
CouenneSolverInterface<T>::CouenneSolverInterface 
(const CouenneSolverInterface &src):

  OsiSolverInterface    (src),
  T                     (src),
  cutgen_               (src.cutgen_),
  knowInfeasible_       (src.knowInfeasible_),
  knowOptimal_          (src.knowOptimal_),
  knowDualInfeasible_   (src.knowDualInfeasible_) {}
//doingResolve_         (src.doingResolve_) {}

/// Destructor
template <class T> 
CouenneSolverInterface<T>::~CouenneSolverInterface () {
  //  if (cutgen_)
  //    delete cutgen_;
}


/// Solve initial LP relaxation 
template <class T> 
void CouenneSolverInterface<T>::initialSolve () {

  knowInfeasible_     = 
  knowOptimal_        = 
  knowDualInfeasible_ = false;

  T::initialSolve ();

  if (T::getObjValue () <= - Couenne_large_bound)
    knowDualInfeasible_ = true;

  // some originals may be unused due to their zero multiplicity (that
  // happens when they are duplicates), restore their value
  if (cutgen_ -> Problem () -> nUnusedOriginals () > 0) {
    CouNumber *x = new CouNumber [T::getNumCols ()];
    CoinCopyN (T::getColSolution (), T::getNumCols (), x);
    cutgen_ -> Problem () -> restoreUnusedOriginals (x);
    T::setColSolution (x);
    delete [] x;
  }
}

template <class T>
bool CouenneSolverInterface<T>::isProvenPrimalInfeasible() const {
  return knowInfeasible_ || T::isProvenPrimalInfeasible();
}

template <class T> 
bool CouenneSolverInterface<T>::isProvenOptimal() const {
  return knowOptimal_ || T::isProvenOptimal();
}

template <class T> 
bool CouenneSolverInterface<T>::isProvenDualInfeasible() const {
  return knowDualInfeasible_ || T::isProvenDualInfeasible();
}

/// Defined in Couenne/src/convex/generateCuts.cpp
void sparse2dense (int, t_chg_bounds *, int *&, int &);


/// Resolve an LP relaxation after problem modification
template <class T> 
void CouenneSolverInterface<T>::resolve () {

  static int count = -1;
  char filename [30];

  // save problem to be loaded later
  if (cutgen_ && (cutgen_ -> check_lp ())) {
    count++;
    sprintf (filename, "resolve_%d", count);
    T::writeMps (filename);
  }

  knowInfeasible_     = 
  knowOptimal_        = 
  knowDualInfeasible_ = false;

  const CoinWarmStart *ws = NULL;

  if (cutgen_ && (cutgen_ -> check_lp ()))
    ws = T::getWarmStart ();

  //deleteScaleFactors ();

  // re-solve problem
  T::resolve ();

  if (T::getObjValue () <= - Couenne_large_bound)
    knowDualInfeasible_ = true;

  CouNumber 
    //objval     = T::getObjValue (),
    curCutoff  = cutgen_ -> Problem () -> getCutOff (),
    objvalGlob = T::getColSolution () [cutgen_ -> Problem () -> Obj (0) -> Body () -> Index ()];  

  // check if resolve found new integer solution
  if (//doingResolve () &&                 // this is not called from strong branching
      isProvenOptimal () &&
      (objvalGlob < curCutoff - COUENNE_EPS) &&
      (cutgen_ -> Problem () -> checkNLP (T::getColSolution (), objvalGlob, true)) &&
      //      (objvalGlo < curCutoff - COUENNE_EPS) && // check again as it may have changed
      (objvalGlob > -COUENNE_INFINITY/2)) {    // check if it makes sense

    // also save the solution so that cbcModel::setBestSolution saves it too

    //printf ("new cutoff from CSI: %g\n", objval);
    cutgen_ -> Problem () -> setCutOff (objvalGlob);
  }

  // some originals may be unused due to their zero multiplicity (that
  // happens when they are duplicates), restore their value
  if (cutgen_ -> Problem () -> nUnusedOriginals () > 0) {
    CouNumber *x = new CouNumber [T::getNumCols ()];
    CoinCopyN (T::getColSolution (), T::getNumCols (), x);
    cutgen_ -> Problem () -> restoreUnusedOriginals (x);
    T::setColSolution (x);
    delete [] x;
  }

  // check LP independently
  if (cutgen_ && (cutgen_ -> check_lp ())) {

    OsiSolverInterface
      *nsi = new T,
      *csi = clone ();

    sprintf (filename, "resolve_%d.mps", count);
    nsi -> readMps (filename);

    nsi -> messageHandler () -> setLogLevel (0);
    nsi -> setWarmStart (ws);

    nsi -> initialSolve ();

    if ((nsi -> isProvenOptimal () && isProvenOptimal ()) ||
	!(nsi -> isProvenOptimal ()) && !isProvenOptimal ()) {

      if (nsi -> isProvenOptimal () &&
	  (fabs (nsi -> getObjValue () - T::getObjValue ()) / 
	   (1. + fabs (nsi -> getObjValue ()) + fabs (T::getObjValue ())) > 1e-2))

	printf ("Warning: discrepancy between saved %g and current %g [%g], file %s\n", 
		nsi -> getObjValue (),  T::getObjValue (),
		nsi -> getObjValue () - T::getObjValue (),
		filename);
    }

    csi -> messageHandler () -> setLogLevel (0);
    csi -> setWarmStart (ws);

    csi -> initialSolve ();

    if ((csi -> isProvenOptimal () && isProvenOptimal ()) ||
	!(csi -> isProvenOptimal ()) && !isProvenOptimal ()) {

      if (csi -> isProvenOptimal () &&
	  (fabs (csi -> getObjValue () - T::getObjValue ()) / 
	   (1. + fabs (csi -> getObjValue ()) + fabs (T::getObjValue ())) > 1e-2))

	printf ("Warning: discrepancy between cloned %g and current %g [%g]\n", 
		csi -> getObjValue (),  T::getObjValue (),
		csi -> getObjValue () - T::getObjValue ());
    }

    delete nsi;
    delete csi;
    
    delete ws;

    //else printf ("Warning: discrepancy between statuses %s -- %s feasible\n", 
    //filename, isProvenOptimal () ? "current" : "saved");
  }
}


/// Create a hot start snapshot of the optimization process.
template <class T> 
void CouenneSolverInterface<T>::markHotStart () 
{OsiSolverInterface::markHotStart ();} // OsiClpSolverInterface::markHotStart() seems not to work


/// Delete the hot start snapshot.
template <class T> 
void CouenneSolverInterface<T>::unmarkHotStart () 
{OsiSolverInterface::unmarkHotStart();}


/// Optimize starting from the hot start snapshot.
template <class T> 
void CouenneSolverInterface<T>::solveFromHotStart () {

  knowInfeasible_     = 
  knowOptimal_        = 
  knowDualInfeasible_ = false;

  resolve();

  if (T::getObjValue () <= - Couenne_large_bound)
    knowDualInfeasible_ = true;

  // some originals may be unused due to their zero multiplicity (that
  // happens when they are duplicates), restore their value
  if (cutgen_ -> Problem () -> nUnusedOriginals () > 0) {
    CouNumber *x = new CouNumber [T::getNumCols ()];
    CoinCopyN (T::getColSolution (), T::getNumCols (), x);
    cutgen_ -> Problem () -> restoreUnusedOriginals (x);
    T::setColSolution (x);
    delete [] x;
  }

  if (isProvenPrimalInfeasible ()) knowInfeasible_     = true;
  if (isProvenOptimal          ()) knowOptimal_        = true;
  if (isProvenDualInfeasible   ()) knowDualInfeasible_ = true;
}

//class CouenneSolverInterface <OsiClpSolverInterface>;
