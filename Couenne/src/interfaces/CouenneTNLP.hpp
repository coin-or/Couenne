/* $Id$
 *
 * Name:    CouenneTNLP.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Definition of an NLP interface with gradient/Jacobian/etc
 * 
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNETNLP_HPP
#define COUENNETNLP_HPP

#include "IpTNLP.hpp"
#include "CouenneExprJac.hpp"
#include "CouenneExprHess.hpp"

#include <vector>
#include <set>

using namespace Ipopt;

namespace Couenne {

  class CouenneProblem;

  /// Class for handling NLPs using CouenneProblem
  class CouenneTNLP: public Ipopt::TNLP {

  public:

    /// Empty constructor
    CouenneTNLP ();

    /// Constructor 
    CouenneTNLP (CouenneProblem *);

    /// Destructor
    virtual ~CouenneTNLP ();

    /// set initial solution
    void setInitSol (double *sol);

    /// returns best solution (if it exists)
    CouNumber *getSolution ()
    {return sol_;}

    /// returns value of the best solution
    CouNumber getSolValue ()
    {return bestZ_;}

    /// return the number of variables and constraints, and the number
    /// of non-zeros in the jacobian and the hessian. The index_style
    /// parameter lets you specify C or Fortran style indexing for the
    /// sparse matrix iRow and jCol parameters.  C_STYLE is 0-based,
    /// and FORTRAN_STYLE is 1-based.
    virtual bool get_nlp_info (Ipopt::Index& n, 
			       Ipopt::Index& m, 
			       Ipopt::Index& nnz_jac_g,
			       Ipopt::Index& nnz_h_lag, 
			       enum Ipopt::TNLP::IndexStyleEnum& index_style);

    /// return the information about the bound on the variables and
    /// constraints. The value that indicates that a bound does not
    /// exist is specified in the parameters nlp_lower_bound_inf and
    /// nlp_upper_bound_inf.  By default, nlp_lower_bound_inf is -1e19
    /// and nlp_upper_bound_inf is 1e19. (see TNLPAdapter)
    virtual bool get_bounds_info (Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
				  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

    /// return the variables linearity (TNLP::Linear or
    /// TNLP::NonLinear). The var_types array should be allocated with
    /// length at least n. (default implementation just return false
    /// and does not fill the array).
    virtual bool get_variables_linearity (Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types);

    /// return the constraint linearity.  array should be alocated
    /// with length at least n. (default implementation just return
    /// false and does not fill the array).
    virtual bool get_constraints_linearity (Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types);

    /// return the starting point. The bool variables indicate whether
    /// the algorithm wants you to initialize x, z_L/z_u, and lambda,
    /// respectively.  If, for some reason, the algorithm wants you to
    /// initialize these and you cannot, return false, which will
    /// cause Ipopt to stop.  You will have to run Ipopt with
    /// different options then.
    virtual bool get_starting_point (Ipopt::Index n, 
				     bool init_x, Ipopt::Number* x,
				     bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
				     Ipopt::Index m, 
				     bool init_lambda, Ipopt::Number* lambda);

    /// return the value of the objective function
    virtual bool eval_f (Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                         Ipopt::Number& obj_value);

    /// return the vector of the gradient of the objective w.r.t. x
    virtual bool eval_grad_f (Ipopt::Index n, const Ipopt::Number* x, 
			      bool new_x,
			      Ipopt::Number* grad_f);

    /// return the vector of constraint values
    virtual bool eval_g (Ipopt::Index n, const Ipopt::Number* x, bool new_x,
			 Ipopt::Index m, Number* g);

    /// return the jacobian of the constraints. The vectors iRow and
    /// jCol only need to be set once. The first call is used to set
    /// the structure only (iRow and jCol will be non-NULL, and values
    /// will be NULL) For subsequent calls, iRow and jCol will be
    /// NULL.
    virtual bool eval_jac_g (Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                             Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
			     Ipopt::Index *jCol, Ipopt::Number* values);

    /// return the hessian of the lagrangian. The vectors iRow and
    /// jCol only need to be set once (during the first call). The
    /// first call is used to set the structure only (iRow and jCol
    /// will be non-NULL, and values will be NULL) For subsequent
    /// calls, iRow and jCol will be NULL. This matrix is symmetric -
    /// specify the lower diagonal only.  A default implementation is
    /// provided, in case the user wants to se quasi-Newton
    /// approximations to estimate the second derivatives and doesn't
    /// not neet to implement this method.
    virtual bool eval_h (Ipopt::Index n, const Ipopt::Number* x, bool new_x,
			 Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
			 bool new_lambda, Ipopt::Index nele_hess,
			 Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);

    /// This method is called when the algorithm is complete so the TNLP can store/write the solution
    virtual void finalize_solution (SolverReturn status,
				    Index n, const Number* x, const Number* z_L, const Number* z_U,
				    Index m, const Number* g, const Number* lambda,
				    Number obj_value,
				    const IpoptData* ip_data,
				    IpoptCalculatedQuantities* ip_cq);

    /// Intermediate Callback method for the user.  Providing dummy
    /// default implementation.  For details see IntermediateCallBack
    /// in IpNLP.hpp.
    virtual bool intermediate_callback (AlgorithmMode mode,
					Index iter, Number obj_value,
					Number inf_pr, Number inf_du,
					Number mu, Number d_norm,
					Number regularization_size,
					Number alpha_du, Number alpha_pr,
					Index ls_trials,
					const IpoptData* ip_data,
					IpoptCalculatedQuantities* ip_cq);

    /// @name Methods for quasi-Newton approximation.  If the second
    /// derivatives are approximated by Ipopt, it is better to do this
    /// only in the space of nonlinear variables.  The following
    /// methods are call by Ipopt if the quasi-Newton approximation is
    /// selected.  If -1 is returned as number of nonlinear variables,
    /// Ipopt assumes that all variables are nonlinear.  Otherwise, it
    /// calls get_list_of_nonlinear_variables with an array into which
    /// the indices of the nonlinear variables should be written - the
    /// array has the lengths num_nonlin_vars, which is identical with
    /// the return value of get_number_of_nonlinear_variables ().  It
    /// is assumed that the indices are counted starting with 1 in the
    /// FORTRAN_STYLE, and 0 for the C_STYLE.
    virtual Index get_number_of_nonlinear_variables ();

    /// get real list
    virtual bool get_list_of_nonlinear_variables (Index num_nonlin_vars,
						  Index* pos_nonlin_vars);

  private:

    /// Pointer to the object containing all info
    CouenneProblem *problem_;

    /// Initial solution
    CouNumber *sol0_;

    /// Optimal solution
    CouNumber *sol_;

    /// Value of the optimal solution
    CouNumber bestZ_;

    /// expression gradient (packed sparse vector)
    std::vector <std::pair <int, expression *> > gradient_;

    /// list of nonlinear variables
    std::set <int> nonLinVars_;

    /// Jacobian
    ExprJac Jac_;

    /// Hessian --- there are 1+m of them, but all are squeezed in a
    /// single object
    ExprHess HLa_;
  };
}

#endif
