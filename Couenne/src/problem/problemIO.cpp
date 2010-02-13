/* $Id$
 *
 * Name:    problemIO.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include <fstream>

#include "CoinHelperFunctions.hpp"

#include "expression.hpp"
#include "exprAux.hpp"
#include "CouenneProblem.hpp"


// output content of the problem
void CouenneProblem::print (std::ostream &out) {

  out << "objectives:" << std::endl;
  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i)
    (*i) -> print (out);

  out << "constraints:" << std::endl;
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i)
    (*i) -> print (out);

  out << "variables:" << std::endl;
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)

    if (((*i) -> Type () != AUX) || 
	((*i) -> Multiplicity () > 0)) {

      (*i) -> print (out);

      if (((*i) -> Type () == AUX) && 
	  ((*i) -> Multiplicity () > 0)) {

	out << " (r:" << (*i) -> rank () 
	    << ", m:" << (*i) -> Multiplicity () << ") := ";

	if ((*i) -> Image ())
	  (*i) -> Image () -> print (out, false); 
      }

      CouNumber 
	lb = domain_.lb ((*i) -> Index ()),
	ub = domain_.ub ((*i) -> Index ());

      if ((fabs (lb)     < COUENNE_EPS) &&
	  (fabs (ub - 1) < COUENNE_EPS) &&
	  (*i) -> isInteger ()) out << " binary";
      else {

	out << " [ " << lb << " , " << ub << " ]";
	if ((*i) -> isInteger ()) out << " integer";
      }

      out << std::endl;
    }

  if (commonexprs_.size ()) {
    out << "common expressions:" << std::endl;
    for (std::vector <expression *>::iterator i = commonexprs_.begin ();
	 i != commonexprs_.end (); ++i) {
      out << "v["    << i - commonexprs_.begin () << "] := ";
      (*i) -> print (out);
      out << std::endl;
    }
  }

  if (optimum_) {
    out << "best known solution: (" << *optimum_;
    for (int i=1; i < nVars (); i++)
      out << ' ' << optimum_ [i];
    out << ')' << std::endl;
  }

  if (fabs (bestObj_) < COUENNE_INFINITY)
    out << "best known objective: " << bestObj_ << std::endl;

  out << "end" << std::endl;
} 


/// read optimal solution into member optimum
bool CouenneProblem::readOptimum (std::string *fname) {

  FILE *f;

  if (fname == NULL) {

    fname = &problemName_;

    int base = fname -> rfind ('/'), size;
    if (base < 0) base = 0; else base++;

    size = fname -> find  ('.', base) - base;

    char *filename = new char [size+5];
    CoinFillN (filename, size+5, (char) 0);
    fname -> copy (filename, 1+size, base);
    strcat (filename, "txt");
    f = fopen (filename, "r");
    delete [] filename;
  } else f = fopen (fname -> c_str (), "r");

  if (!f) return false;

  optimum_ = (CouNumber *) realloc (optimum_, nVars () * sizeof (CouNumber));

  CoinFillN (optimum_, nVars (), 0.);

  // read optimal objective function first
  if (fscanf (f, "%lf", &bestObj_) < 1) {
    fclose (f);
    printf ("could not read objective from file \"%s\"\n", fname -> c_str ());
    return false;
  }

  // read optimal values of variables
  for (int i = 0; i < nOrigVars_; i++)
    if (fscanf (f, "%lf", optimum_ + i) < 1) {
      fclose (f);
      printf ("could not read optimal value of x_%d from file \"%s\"\n", i, fname -> c_str ());
      return false;
    }

  if (opt_window_ < 1e50) // restrict solution space around known optimum
    for (int i = 0; i < nOrigVars_; i++) {
      Lb (i) = CoinMax (Lb (i), optimum_ [i] - opt_window_ * (1 + fabs (optimum_ [i])));
      Ub (i) = CoinMin (Ub (i), optimum_ [i] + opt_window_ * (1 + fabs (optimum_ [i])));
    }

  // expand solution to auxiliary space
  getAuxs (optimum_);

  fclose (f);
  return true;
}


/// read cutoff into member optimum
void CouenneProblem::readCutoff (const std::string &fname) {

  CouNumber val;

  FILE *f = fopen (fname.c_str (), "r");
  if (!f) return;

  if (fscanf (f, "%lf", &val) < 1) 
    return;

  setCutOff (val);
}
