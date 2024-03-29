/*
 *
 * Name:    domain.hpp
 * Author:  Pietro Belotti
 * Purpose: class for point and bounding box
 *
 * (C) Carnegie-Mellon University, 2008-09.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_DOMAIN_HPP
#define COUENNE_DOMAIN_HPP

#include <stdlib.h>
#include <stack>

#include "CouenneTypes.hpp"
#include "CouenneConfig.h"

namespace Osi {

  class OsiSolverInterface;
  class OsiCuts;
}


namespace Couenne {

/// Define a point in the solution space and the bounds around it.

class COUENNELIB_EXPORT DomainPoint {

  friend class Domain;

protected:

  int dimension_; ///< dimension of point

  CouNumber *x_;  ///< current value of variables
  CouNumber *lb_; ///< lower bound
  CouNumber *ub_; ///< upper bound

  bool copied_;   ///< true if data has been copied (so we own it, and
		  ///  have to delete it upon destruction)

  bool isNlp_;    ///< true if this point comes from an NLP solver
		  ///  (and is thus nlp feasible)
public:

  /// constructor
  DomainPoint (int dim,
	       CouNumber *x,
	       CouNumber *lb,
	       CouNumber *ub,
	       bool copy = true);

  /// constructor
  DomainPoint (int dim = 0,
	       const CouNumber *x   = NULL,
	       const CouNumber *lb  = NULL,
	       const CouNumber *ub  = NULL,
	       bool copy = true);

  /// destructor
  ~DomainPoint () {
    if (copied_) {
      if (x_)  free (x_);
      if (lb_) free (lb_);
      if (ub_) free (ub_);
    }
  }

  /// copy constructor
  DomainPoint (const DomainPoint &src);

  /// resize domain point (for extending into higher space)
  void resize (int newdim);

  /// return current size
  int size () const {return dimension_;}

  /// return dimension_
  inline int Dimension () {return dimension_;}

  inline CouNumber &x  (int index) {return x_  [index];} ///< return current variable
  inline CouNumber &lb (int index) {return lb_ [index];} ///< return current lower bound
  inline CouNumber &ub (int index) {return ub_ [index];} ///< return current upper bound

  inline CouNumber *x  () {return x_ ;} ///< return current variable vector
  inline CouNumber *lb () {return lb_;} ///< return current lower bound vector
  inline CouNumber *ub () {return ub_;} ///< return current upper bound vector

  /// assignment operator
  DomainPoint &operator= (const DomainPoint &src);

  /// true if this point is the nlp solution
  bool &isNlp ()
  {return isNlp_;}
};


/// Define a dynamic point+bounds, with a way to save and restore
/// previous points+bounds through a LIFO structure

class COUENNELIB_EXPORT Domain {

protected:

  DomainPoint *point_;                  ///< current point
  std::stack <DomainPoint *> domStack_; ///< stack of saved points

public:

  /// basic constructor
  Domain (): point_ (NULL) {}

  /// copy constructor
  Domain (const Domain &src) {
    point_ = new DomainPoint (*(src.point_));
    // TODO -- not important, discard previous points when copying problem
    /*for (std::stack <DomainPoint *>::iterator i = src.domStack_.begin ();
	 i != src.domStack_.end (); ++i)
	 domStack_.push (new DomainPoint (**i));*/
  }

  /// destructor
  ~Domain ();

  /// save current point and start using another
  void push (int dim,
	     CouNumber *x,
	     CouNumber *lb,
	     CouNumber *ub,
	     bool copy = true);

  /// save current point and start using another
  void push (int dim,
	     const CouNumber *x,
	     const CouNumber *lb,
	     const CouNumber *ub,
	     bool copy = true);

  /// save current point and start using another -- retrieve
  /// information from solver interface and from previous column cuts
  void push (const OsiSolverInterface *si,
	     OsiCuts *cs = NULL,
	     bool copy = true);

  /// save current point and start using another
  void push (const DomainPoint &dp, bool copy = true);

  /// restore previous point
  void pop ();

  inline DomainPoint *current ()   {return point_;}                          ///< return current point

  inline CouNumber &x  (int index) {return point_ -> x  (index);}   ///< current variable
  inline CouNumber &lb (int index) {return point_ -> lb (index);}   ///< current lower bound
  inline CouNumber &ub (int index) {return point_ -> ub (index);}   ///< current upper bound

  inline CouNumber *x  () {return point_ -> x  ();}   ///< return current variable vector
  inline CouNumber *lb () {return point_ -> lb ();}   ///< return current lower bound vector
  inline CouNumber *ub () {return point_ -> ub ();}   ///< return current upper bound vector
};

}

#endif
