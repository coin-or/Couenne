/* $Id$
 *
 * Name:    CouenneSolverInterface.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 *          Andreas Waechter, IBM Corp.
 * Purpose: Implementation of the OsiSolverInterface::resolve () method 
 *
 * (C) Carnegie-Mellon University, 2006-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "OsiSolverInterface.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneCutGenerator.hpp"

#include "CouenneRecordBestSol.hpp"                         

//#define FM_CHECK

namespace Couenne {

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

  // printf ("init solution: (");
  // for (int i=0; i< T::getNumCols (); i++)
  //   printf ("%g [%g,%g] ", 
  // 	    T::getColSolution () [i],
  // 	    T::getColLower    () [i],
  // 	    T::getColUpper    () [i]);
  // printf (")\n");

  if (getObjValue () <= - Couenne_large_bound)
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

  // printf ("solution: (");
  // for (int i=0; i< T::getNumCols (); i++)
  //   printf ("%g ", T::getColSolution () [i]);
  // printf (")\n");

  if (getObjValue () <= - Couenne_large_bound)
    knowDualInfeasible_ = true;

  int objind = cutgen_ -> Problem () -> Obj (0) -> Body () -> Index ();

  CouNumber 
    //objval     = T::getObjValue (),
    curCutoff  = cutgen_ -> Problem () -> getCutOff (),
    objvalGlob = objind >= 0 ? T::getColSolution () [objind] : cutgen_ -> Problem () -> Obj (0) -> Body () -> Value ();

  // check if resolve found new integer solution
  bool isChecked = false;  
#ifdef FM_CHECKNLP2
  double curBestVal = 1.e50;
  if(cutgen_->Problem()->getRecordBestSol()->getHasSol()) { 
    curBestVal =  cutgen_->Problem()->getRecordBestSol()->getVal(); 
  }
  curBestVal = (curBestVal < curCutoff ? curBestVal : curCutoff);
  if(isProvenOptimal()) {
    isChecked = cutgen_->Problem()->checkNLP2(T::getColSolution(), 
					      curBestVal, false,
                                              true, // stopAtFirstViol
                                              true, // checkALL
					      cutgen_->Problem()->getFeasTol());
    if(isChecked) {
      objvalGlob = cutgen_->Problem()->getRecordBestSol()->getModSolVal();
      if(!(objvalGlob < curBestVal - COUENNE_EPS)) {
        isChecked = false; 
      }
    }
  }

#ifdef FM_CHECK
  bool ckIsChecked = false;
  double ckObj = objvalGlob;
  if(isProvenOptimal () &&
     (objvalGlob < curCutoff - COUENNE_EPS)) {
    ckIsChecked = cutgen_->Problem()->checkNLP(T::getColSolution (),
					       ckObj, true);
  }
  if(!isChecked && ckIsChecked) {
    printf("CouenneSolverInterface::resolve(): ### ERROR: isChecked: false  ckIsChecked: true\n");
    exit(1);
  }
  else {
    printf("CouenneSolverInterface::resolve(): isChecked == ckIsChecked\n");
  }
#endif

#else /* not FM_CHECKNLP2 */
  if(isProvenOptimal () &&
     (objvalGlob < curCutoff - COUENNE_EPS)) {
    isChecked = cutgen_->Problem()->checkNLP(T::getColSolution (),
                                             objvalGlob, true);
  }
#endif /* not FM_CHECKNLP2 */

  if (//doingResolve () &&    // this is not called from strong branching
      isChecked &&
      (objvalGlob > -COUENNE_INFINITY/2)) {    // check if it makes sense

    // also save the solution so that cbcModel::setBestSolution saves it too

    //printf ("new cutoff from CSI: %g\n", objval);
    cutgen_ -> Problem () -> setCutOff (objvalGlob);

#ifdef FM_TRACE_OPTSOL
#ifdef FM_CHECKNLP2
    cutgen_->Problem()->getRecordBestSol()->update();
#else /* not FM_CHECKNLP2 */

  // some originals may be unused due to their zero multiplicity (that
  // happens when they are duplicates), restore their value
  if (cutgen_ -> Problem () -> nUnusedOriginals () > 0) {
    CouNumber *x = new CouNumber [T::getNumCols ()];
    CoinCopyN (T::getColSolution (), T::getNumCols (), x);
    cutgen_ -> Problem () -> restoreUnusedOriginals (x);
    T::setColSolution (x);
    delete [] x;
  }

  cutgen_->Problem()->getRecordBestSol()->update(T::getColSolution(), 
						 cutgen_->Problem()->nVars(),
						 objvalGlob,
						 cutgen_->Problem()->getFeasTol());
#endif  /* not FM_CHECKNLP2 */
#endif /* FM_TRACE_OPTSOL */

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
	(!(nsi -> isProvenOptimal ()) && !isProvenOptimal ())) {

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
	(!(csi -> isProvenOptimal ()) && !isProvenOptimal ())) {

      if (csi -> isProvenOptimal () &&
	  (fabs (csi -> getObjValue () - getObjValue ()) / 
	   (1. + fabs (csi -> getObjValue ()) + fabs (getObjValue ())) > 1e-2))

	printf ("Warning: discrepancy between cloned %g and current %g [%g]\n", 
		csi -> getObjValue (),  getObjValue (),
		csi -> getObjValue () - getObjValue ());
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

  if (getObjValue () <= - Couenne_large_bound)
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

/// Get the objective function value. Modified due to possible
/// constant objectives passed to Couenne
template <class T> 
inline double CouenneSolverInterface<T>::getObjValue() const {
  return cutgen_ -> Problem () -> constObjVal () + T::getObjValue();
}


}
