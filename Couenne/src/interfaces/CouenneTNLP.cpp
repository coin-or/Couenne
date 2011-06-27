/* $Id$
 *
 * Name:    CouenneTNLP.cpp
 * Authors: Pietro Belotti
 * Purpose: Implementation of an NLP interface with gradient/Jacobian/etc
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "CouenneSparseMatrix.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneExprVar.hpp"
#include "CouenneExprJac.hpp"
#include "CouenneExprHess.hpp"
#include "CouenneTNLP.hpp"

#include <stdio.h>

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

//#define DEBUG

using namespace Ipopt;
using namespace Couenne;

/// Empty constructor
CouenneTNLP::CouenneTNLP ():

  problem_        (NULL),
  sol0_           (NULL),
  sol_            (NULL),
  HLa_            (NULL),

  optHessian_     (NULL),
  saveOptHessian_ (false) {}


/// Empty constructor
CouenneTNLP::~CouenneTNLP () {

  if (sol0_)       delete [] sol0_;
  if (sol_)        delete [] sol_;
  if (HLa_)        delete HLa_;
  if (optHessian_) delete optHessian_;

  for (std::vector <std::pair <int, expression *> >::iterator i = gradient_. begin (); 
       i != gradient_. end (); ++i)
    delete (*i). second;
}


/// Constructor 
CouenneTNLP::CouenneTNLP (CouenneProblem *p):

  problem_        (p),
  sol0_           (NULL),
  sol_            (NULL),
  bestZ_          (COIN_DBL_MAX),
  Jac_            (p),
  HLa_            (new ExprHess (p)),
  optHessian_     (NULL),
  saveOptHessian_ (false) {

  std::set <int> objDep;

  expression *obj = problem_ -> Obj (0) -> Body ();

  // objective of entering problem is reformulated, no need to go
  // further
  obj -> DepList (objDep, STOP_AT_AUX);

  for (std::set <int>::iterator i = objDep.begin (); i != objDep. end (); ++i) {

    expression *gradcomp = obj -> differentiate (*i);
    gradcomp -> realign (problem_);
    gradient_ . push_back (std::pair <int, expression *> (*i, gradcomp));
  }

  // create data structures for nonlinear variables (see
  // get_[number|list]_of_nonlinear_variables () below)

  // constraints

  for (int i = 0; i < problem_ -> nCons (); i++) {

    expression *e = problem_ -> Con (i) -> Body ();

    // if constraint is single variable, don't treat it as constraint
    // but rather as variable bound

    if (e -> Type      () == AUX ||
	e -> Type      () == VAR ||
	e -> Linearity () <= LINEAR)
      continue;

    // constraint is nonlinear, get all variables its left-hand side
    // depends on and make them nonlinear

    e -> DepList (nonLinVars_, STOP_AT_AUX);
  }

  // auxiliaries

  for (int i = 0; i < problem_ -> nVars (); i++) {

    exprVar *e = problem_ -> Var (i);

    if ((e -> Type         () != AUX)    ||
	(e -> Multiplicity () <= 0)      ||
	(e -> Linearity    () <= LINEAR))
      continue;

    e -> Image () -> DepList (nonLinVars_, STOP_AT_AUX);
  }
}


/// Copy constructor 
CouenneTNLP::CouenneTNLP (const CouenneTNLP &rhs) 
{operator= (rhs);}


/// Assignment 
CouenneTNLP &CouenneTNLP::operator= (const CouenneTNLP &rhs) {

  problem_            = rhs.problem_;

  sol0_               = rhs.sol0_ && problem_ ? CoinCopyOfArray (rhs.sol0_, problem_ -> nVars ()) : NULL;
  sol_                = rhs.sol_  && problem_ ? CoinCopyOfArray (rhs.sol_,  problem_ -> nVars ()) : NULL;

  bestZ_              = rhs.bestZ_;
  gradient_           = rhs.gradient_;
  nonLinVars_         = rhs.nonLinVars_;

  Jac_                = rhs.Jac_;
  HLa_                = rhs.HLa_ ? new ExprHess (*(rhs.HLa_)) : NULL;

  optHessian_         = rhs.optHessian_ ? new CouenneSparseMatrix (*(rhs.optHessian_)) : NULL;
  saveOptHessian_     = rhs.saveOptHessian_;

  return *this;
}

/// Clone
CouenneTNLP *CouenneTNLP::clone () 
{return new CouenneTNLP (*this);}


// overload this method to return the number of variables and
// constraints, and the number of non-zeros in the jacobian and the
// hessian. The index_style parameter lets you specify C or Fortran
// style indexing for the sparse matrix iRow and jCol parameters.
// C_STYLE is 0-based, and FORTRAN_STYLE is 1-based.
bool CouenneTNLP::get_nlp_info (Index& n, 
				Index& m, 
				Index& nnz_jac_g,
			        Index& nnz_h_lag, 
				IndexStyleEnum& index_style) {
  n = problem_ -> nVars ();
  m = Jac_. nRows ();

  nnz_jac_g = Jac_ .  nnz ();
  nnz_h_lag = HLa_ -> nnz ();

  index_style = C_STYLE; // what else? ;-)

  return true;
}


/// set initial solution
void CouenneTNLP::setInitSol (double *sol) {

  if (sol) {
    if (!sol0_)
      sol0_ = new CouNumber [problem_ -> nVars ()];
    CoinCopyN (sol, problem_ -> nVars (), sol0_);
  }
}

// overload this method to return the information about the bound on
// the variables and constraints. The value that indicates that a
// bound does not exist is specified in the parameters
// nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
// nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is 1e19. (see
// TNLPAdapter)
bool CouenneTNLP::get_bounds_info (Index n, Number* x_l, Number* x_u,
				   Index m, Number* g_l, Number* g_u) {

  // constraints

#ifdef DEBUG
  printf ("get_bounds_info on %d cons, %d vars\n", m, n);
#endif

  for (int i = 0; i < problem_ -> nCons (); i++) {

    CouenneConstraint *c = problem_ -> Con (i);

    if (c -> Body () -> Type () == AUX ||
	c -> Body () -> Type () == VAR) 
      continue;

    CouNumber
      clb = (*c -> Lb ()) (),
      cub = (*c -> Ub ()) ();

    // prevent ipopt from exiting on inconsistent bounds
    if (clb <= cub) {*g_l++ = clb; *g_u++ = cub;} 
    else            {*g_l++ = cub; *g_u++ = clb;}
  }

  // auxiliaries

  for (int i = 0; i < problem_ -> nVars (); i++) {

    exprVar *e = problem_ -> Var (i);

    CouNumber
      lb = e -> lb (),
      ub = e -> ub ();

    // prevent ipopt from exiting on inconsistent bounds
    if (lb <= ub) {*x_l++ = lb; *x_u++ = ub;} 
    else          {*x_l++ = ub; *x_u++ = lb;}

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    *g_l = (e -> sign () != expression::AUX_GEQ) ? 0. : -COIN_DBL_MAX;
    *g_u = (e -> sign () != expression::AUX_LEQ) ? 0. :  COIN_DBL_MAX;

    ++g_l;
    ++g_u;
  }

  return true;
}


// overload this method to return the variables linearity
// (TNLP::LINEAR or TNLP::NON_LINEAR). The var_types array should be
// allocated with length at least n. (default implementation just
// returns false and does not fill the array).
bool CouenneTNLP::get_variables_linearity (Index n, Ipopt::TNLP::LinearityType* var_types) {

  CoinFillN (var_types, n, Ipopt::TNLP::LINEAR);

  for (std::set <int>:: iterator i = nonLinVars_. begin (); i != nonLinVars_. end (); ++i)
    var_types [*i] = Ipopt::TNLP::NON_LINEAR;

  return true;
}

// overload this method to return the constraint linearity.  array
// should be alocated with length at least n. (default implementation
// just returns false and does not fill the array).
bool CouenneTNLP::get_constraints_linearity (Index m, Ipopt::TNLP::LinearityType* const_types) {

  // constraints

  for (int i = 0; i < problem_ -> nCons (); i++) {

    expression *b = problem_ -> Con (i) -> Body ();

    if (b -> Type () == AUX ||
	b -> Type () == VAR) 
      continue;

    *const_types++ = 
      (b -> Linearity () > LINEAR) ? 
      Ipopt::TNLP::NON_LINEAR : 
      Ipopt::TNLP::LINEAR;
 }

  // auxiliaries

  for (int i = 0; i < problem_ -> nVars (); i++) {

    exprVar *e = problem_ -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    *const_types++ = 
      (e -> Image () -> Linearity () > LINEAR) ? 
      Ipopt::TNLP::NON_LINEAR : 
      Ipopt::TNLP::LINEAR;
  }

  return true;
}


// overload this method to return the starting point. The bool
// variables indicate whether the algorithm wants you to initialize x,
// z_L/z_u, and lambda, respectively.  If, for some reason, the
// algorithm wants you to initialize these and you cannot, return
// false, which will cause Ipopt to stop.  You will have to run Ipopt
// with different options then.
bool CouenneTNLP::get_starting_point (Index n, 
				      bool init_x, Number* x,
				      bool init_z, Number* z_L, Number* z_U,
				      Index m, 
				      bool init_lambda, Number* lambda) {
  if (init_x)
    CoinCopyN (sol0_, n, x);

  assert (!init_z);      // can't initialize bound multipliers
  assert (!init_lambda); // can't initialize Lagrangian multipliers

  return true;
}


// overload this method to return the value of the objective function
bool CouenneTNLP::eval_f (Index n, const Number* x, bool new_x,
			  Number& obj_value) {
  if (new_x)
    CoinCopyN (x, n, problem_ -> X ()); // can't push domain as we
					// don't know when to pop

  obj_value = (*(problem_ -> Obj (0) -> Body ())) ();
  return true;
}


// overload this method to return the vector of the gradient of
// the objective w.r.t. x
bool CouenneTNLP::eval_grad_f (Index n, const Number* x, bool new_x,
			       Number* grad_f) {

#ifdef DEBUG
  printf ("eval_grad_f: [");
  for (int i=0; i<n; i++)
    printf ("%.2g ", x [i]);
  printf ("] --> [");
#endif

  if (new_x)
    CoinCopyN (x, n, problem_ -> X ()); // can't push domain as we
					// don't know when to pop

  CoinFillN (grad_f, n, 0.);

  for (std::vector <std::pair <int, expression *> >::iterator i = gradient_. begin (); 
       i != gradient_. end (); ++i)
    grad_f [i -> first] = (*(i -> second)) ();

#ifdef DEBUG
  for (int i=0; i<n; i++)
    printf ("%.2g ", grad_f [i]);
  printf ("]\n");
#endif

  return true;
}


// overload this method to return the vector of constraint values
bool CouenneTNLP::eval_g (Index n, const Number* x, bool new_x,
			  Index m, Number* g) {

  if (new_x)
    CoinCopyN (x, n, problem_ -> X ()); // can't push domain as we
					// don't know when to pop

#ifdef DEBUG
  if (x) {
    printf ("eval_g: [");
    for (int i=0; i<n; i++)
      printf ("%.2g ", x [i]);
    printf ("] --> [");
  }
#endif

  int nEntries = 0; // FIXME: needs to go

  for (int i = 0; i < problem_ -> nCons (); i++) {

    expression *b = problem_ -> Con (i) -> Body ();

    if (b -> Type () == AUX ||
	b -> Type () == VAR) 
      continue;

    nEntries ++;

    *g++ = (*b) (); // this element of g is the evaluation of the constraint
  }

  // auxiliaries

  assert (n == problem_ -> nVars ());

  for (int i = 0; i < problem_ -> nVars (); i++) {

    exprVar *e = problem_ -> Var (i);

    if ((e -> Type () != AUX) ||
	(e -> Multiplicity () <= 0))
      continue;

    *g++ = (*(e -> Image ())) () - (*e) ();

    nEntries ++;
  }

#ifdef DEBUG
  if (x) {
    for (int i=0; i<nEntries; i++)
      printf ("%.2g ", *(g - nEntries + i));
    printf ("]\n");
  }
#endif

  return true;
}


// overload this method to return the jacobian of the constraints. The
// vectors iRow and jCol only need to be set once. The first call is
// used to set the structure only (iRow and jCol will be non-NULL, and
// values will be NULL) For subsequent calls, iRow and jCol will be
// NULL.
bool CouenneTNLP::eval_jac_g (Index n, const Number* x, bool new_x,
			      Index m, Index nele_jac, Index* iRow,
			      Index *jCol, Number* values) {
  if (new_x)
    CoinCopyN (x, n, problem_ -> X ()); // can't push domain as we
					// don't know when to pop

#ifdef DEBUG
  if (x) {
    printf ("eval_jac_g: ["); fflush (stdout);
    for (int i=0; i<n; i++) 
      {printf ("%.2g ", x [i]); fflush (stdout);}
    printf ("] --> ["); fflush (stdout);
  }
#endif

  if (values == NULL && 
      iRow   != NULL && 
      jCol   != NULL) {

    // initialization of the Jacobian's structure. This has been
    // already prepared by the constructor, so simply copy it

    CoinCopyN (Jac_.iRow (), nele_jac, iRow);
    CoinCopyN (Jac_.jCol (), nele_jac, jCol);

  } else {

    // fill in Jacobian's values. Evaluate each member using the
    // domain modified above by the new value of x

    register expression **e = Jac_. expr ();

    for (register int i=nele_jac; i--;)
      *values++ = (**(e++)) ();
  }

#ifdef DEBUG
  if (values) {
    for (int i=0; i<nele_jac; i++)
      {printf ("%.2g ", *(values - nele_jac + i)); fflush (stdout);}
    printf ("]\n");
  } else printf ("empty\n");
#endif

  return true;
}


// Overload this method to return the hessian of the lagrangian. The
// vectors iRow and jCol only need to be set once (during the first
// call). The first call is used to set the structure only (iRow and
// jCol will be non-NULL, and values will be NULL) For subsequent
// calls, iRow and jCol will be NULL.
//
// This matrix is symmetric - specify the lower diagonal only.
//
// A default implementation is provided, in case the user wants to use
// quasi-Newton approximations to estimate the second derivatives and
// doesn't need to implement this method.
bool CouenneTNLP::eval_h (Index n, const Number* x,      bool new_x,      Number obj_factor, 
			  Index m, const Number* lambda, bool new_lambda, 
			  Index nele_hess,
			  Index* iRow, Index* jCol, Number* values) {

  if (new_x)
    CoinCopyN (x, n, problem_ -> X ()); // can't push domain as we
					// don't know when to pop

#ifdef DEBUG
  if (x) {
    printf ("eval_h: ["); fflush (stdout);
    for (int i=0; i<n; i++)
      {printf ("%.2g ", x [i]); fflush (stdout);}
    printf ("], lambda: ["); fflush (stdout);
    for (int i=0; i<m; i++)
      {printf ("%.2g ", lambda [i]); fflush (stdout);}
    printf ("] --> ["); fflush (stdout);
  }
#endif

  if (values == NULL && 
      iRow   != NULL && 
      jCol   != NULL) {

    /// first call, must determine structure iRow/jCol

    CoinCopyN (HLa_ -> iRow (), nele_hess, iRow);
    CoinCopyN (HLa_ -> jCol (), nele_hess, jCol);

  } else {

    /// generic call, iRow/jCol are known and we should fill in the
    /// values

    CoinZeroN (values, nele_hess);

    for (int i=0; i<nele_hess; i++, values++) {

      int 
	 numL  = HLa_ -> numL () [i],
	*lamI  = HLa_ -> lamI () [i];

      expression
	**expr = HLa_ -> expr () [i];

#ifdef DEBUG
      printf ("[%d %d] %d lambdas: ", HLa_ -> iRow () [i], HLa_ -> jCol () [i], numL); fflush (stdout);
      for (int k=0; k<numL; k++) {
	printf ("%d ", lamI [k]);
	fflush (stdout);
	expr [k] -> print ();
	printf ("\n");
      }
#endif

      // the objective's part of the Hessian can only have level index 0, avoid check

      if (0 == *lamI) {*values += obj_factor           * (*(*expr++)) (); --numL; ++lamI;}
      while (numL--)   *values += lambda [*lamI++ - 1] * (*(*expr++)) ();
    }
  }


#ifdef DEBUG
  if (values) {
    for (int i=0; i<nele_hess; i++)
      {printf ("%.2g ", *(values - nele_hess + i)); fflush (stdout);}
    printf ("]\n");
  } else printf ("empty\n");
#endif

  return true;
}

// Change objective function and modify gradient expressions
// accordingly
void CouenneTNLP::setObjective (expression *newObj) {

  if (HLa_)
    delete HLa_;

  // change the Hessian accordingly

  HLa_ = new ExprHess (problem_);

  std::set <int> objDep;

  // objective of entering problem is reformulated, no need to go
  // further
  newObj -> DepList (objDep, STOP_AT_AUX);

  for (std::vector <std::pair <int, expression *> >::iterator i = gradient_. begin (); 
       i != gradient_. end (); ++i)
    delete (*i). second;

  gradient_ . erase (gradient_ . begin (), gradient_ . end ());

  for (std::set <int>::iterator i = objDep.begin (); i != objDep. end (); ++i) {

    expression
      *gradcomp = newObj -> differentiate (*i),
      *gsimp    = gradcomp -> simplify ();

    if (gsimp) {
      delete gradcomp;
      gradcomp = gsimp;
    }

    gradcomp -> realign (problem_);
    gradient_ . push_back (std::pair <int, expression *> (*i, gradcomp));
  }
}

// This method is called when the algorithm is complete so the TNLP
// can store/write the solution
void CouenneTNLP::finalize_solution (SolverReturn status,
				     Index n, const Number* x, const Number* z_L, const Number* z_U,
				     Index m, const Number* g, const Number* lambda,
				     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq) {

  //printf ("Ipopt[FP] solution (card %d): %12e\n", n, obj_value);

  bestZ_ = obj_value;

  if  (sol_)  CoinCopyN       (x, n, sol_);
  else sol_ = CoinCopyOfArray (x, n);

  // if a save-flag was set, save this solution's lagrangian hessian
  // for later use by the FP

  if (saveOptHessian_) {

    if (!optHessian_)
      optHessian_ = new CouenneSparseMatrix;

    problem_ -> domain () -> push (n, x, NULL, NULL);

    int nnz = HLa_ -> nnz ();

    // resize them to full size (and realloc them to optHessianNum_ later)

    double *&optHessianVal = optHessian_ -> val ();
    int    *&optHessianRow = optHessian_ -> row ();
    int    *&optHessianCol = optHessian_ -> col ();

    int     &optHessianNum = optHessian_ -> num ();    

    optHessianVal = (double *) realloc (optHessianVal, nnz * sizeof (double));
    optHessianRow = (int    *) realloc (optHessianRow, nnz * sizeof (int));
    optHessianCol = (int    *) realloc (optHessianCol, nnz * sizeof (int));

    optHessianNum = 0;    

    for (int i=0; i < HLa_ -> nnz (); ++i) {

      double hessMember = 0.;
      expression **elist = HLa_ -> expr () [i];

      for (int j=0; j < HLa_ -> numL () [i]; ++j) {

	int indLam = HLa_ -> lamI () [i][j];

	hessMember += (indLam == 0) ? 
	  (*(elist [j])) () :                  // this is the objective
	  (*(elist [j])) () * lambda [indLam-1]; // this is a constraint
      }

      if (hessMember != 0.) {

	// printf ("saving: %d, %d --> %g\n", 
	// 	HLa_ -> iRow () [i],
	// 	HLa_ -> jCol () [i], hessMember);

	optHessianVal [optHessianNum]   = hessMember;
	optHessianRow [optHessianNum]   = HLa_ -> iRow () [i];
	optHessianCol [optHessianNum++] = HLa_ -> jCol () [i];
      }
    }

    optHessianVal = (double *) realloc (optHessianVal, optHessianNum * sizeof (double));
    optHessianRow = (int    *) realloc (optHessianRow, optHessianNum * sizeof (int));
    optHessianCol = (int    *) realloc (optHessianCol, optHessianNum * sizeof (int));

    problem_ -> domain () -> pop ();
  }
}


// Intermediate Callback method for the user.  Providing dummy default
// implementation.  For details see IntermediateCallBack in IpNLP.hpp.
bool CouenneTNLP::intermediate_callback (AlgorithmMode mode,
					 Index iter, Number obj_value,
					 Number inf_pr, Number inf_du,
					 Number mu, Number d_norm,
					 Number regularization_size,
					 Number alpha_du, Number alpha_pr,
					 Index ls_trials,
					 const IpoptData* ip_data,
					 IpoptCalculatedQuantities* ip_cq) {

  //printf ("Ipopt FP: iter %4d obj %12e %12e %12e\n", iter, obj_value, inf_pr, inf_du);
  return true;
}


// Methods for quasi-Newton approximation.  If the second derivatives
// are approximated by Ipopt, it is better to do this only in the
// space of nonlinear variables.  The following methods are called by
// Ipopt if the quasi-Newton approximation is selected.  If -1 is
// returned as number of nonlinear variables, Ipopt assumes that all
// variables are nonlinear.
Index CouenneTNLP::get_number_of_nonlinear_variables () 
{return nonLinVars_. size ();}


// Otherwise, it calls get_list_of_nonlinear_variables with an array
// into which the indices of the nonlinear variables should be written
// - the array has the lengths num_nonlin_vars, which is identical
// with the return value of get_number_of_nonlinear_variables().  It
// is assumed that the indices are counted starting with 1 in the
// FORTRAN_STYLE, and 0 for the C_STYLE.
bool CouenneTNLP::get_list_of_nonlinear_variables (Index  num_nonlin_vars,
						   Index* pos_nonlin_vars) {

  for (std::set <int>:: iterator i = nonLinVars_. begin (); i != nonLinVars_. end (); ++i)
    *pos_nonlin_vars++ = *i;

  return true;
}
