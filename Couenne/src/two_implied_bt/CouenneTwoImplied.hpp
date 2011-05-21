/* $Id$
 *
 * Name:    CouenneTwoImplied.hpp
 * Author:  Pietro Belotti
 * Purpose: Bound Tightening using pairs of linear inequalities or equations
 *
 * (C) Pietro Belotti, 2010.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNETWOIMPLIED_HPP
#define COUENNETWOIMPLIED_HPP

#include "BonRegisteredOptions.hpp"

#include "CglCutGenerator.hpp"
#include "OsiRowCut.hpp"
#include "CouenneJournalist.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
}

namespace Couenne {

  class CouenneProblem;

  /**
     Cut Generator for implied bounds derived from pairs of linear (in)equalities

     Implied bounds usually work on a SINGLE inequality of the form
  
      \f$ \ell_j \le   \sum_{i \in N_+} a_{ji} x_i 
                   + \sum_{i \in N_-} a_{ji} x_i \le u_j \f$ 
  
     where  \f$ a_{ji} > 0 \f$  for  \f$ i \in N_+ \f$  and  \f$ a_{ji} < 0 \f$  for  \f$ i \in N_- \f$ , and allow
     one to infer better bounds  \f$ [x^L_i, x^U_i] \f$  on all variables with
     nonzero coefficients:
  
     (1)  \f$  x^L_i \ge (\ell_j - \sum_{i \in N_+} a_{ji} x^U_i 
                                 - \sum_{i \in N_-} a_{ji} x^L_i ) / a_{ji} \qquad \forall i \in N_+   \f$ 
   
     (2)  \f$  x^U_i \le (u_j - \sum_{i \in N_+} a_{ji} x^L_i 
                              - \sum_{i \in N_-} a_{ji} x^U_i ) / a_{ji} \qquad \forall i \in N_+   \f$ 
   
   
     (3)  \f$  x^L_i \ge (u_j - \sum_{i \in N_+} a_{ji} x^L_i 
                              - \sum_{i \in N_-} a_{ji} x^U_i ) / a_{ji} \qquad \forall i \in N_-   \f$ 
   
     (4)  \f$  x^U_i \le (\ell_j - \sum_{i \in N_+} a_{ji} x^U_i 
                                 - \sum_{i \in N_-} a_{ji} x^L_i ) / a_{ji} \qquad \forall i \in N_+  \f$ 
   
     Consider now two inequalities:
  
      \f$ \ell_h \le   \sum_{i \in N^1_+} a_{hi} x_i  + \sum_{i \in N^1_-} a_{hi} x_i \le u_h \f$ 
  
      \f$ \ell_k \le   \sum_{i \in N^2_+} a_{ki} x_i  + \sum_{i \in N^2_-} a_{ki} x_i \le u_k \f$ 
   
     and their CONVEX combination using \f$ \alpha \f$ and \f$ 1 -
     \alpha \f$ , where \f$ \alpha \in [0,1] \f$ :
  
      \f$ \ell' \le \sum_{i \in N} b_i x_i \le u' \f$ 
  
     with  \f$ N = N^1_+\cup N^1_-\cup N^2_+\cup N^2_- \f$ ,  \f$ \ell' =
     \alpha \ell_h + (1-\alpha) \ell_k \f$ , and  \f$ u' = \alpha u_h + (1-\alpha) u_k \f$ . As
     an example where this might be useful, consider
  
      \f$ x + y \ge 2 \f$ 

      \f$ x - y \ge 1 \f$ 
  
     with  \f$ x \in [0,4] \f$  and  \f$ y \in [0,1] \f$ . (This is similar to an example
     given in Tawarmalani and Sahinidis to explain FBBT != OBBT, I
     believe.)  The sum of the two above inequalities gives  \f$ x \ge 1.5 \f$ ,
     while using only the implied bounds on the single inequalities
     gives  \f$ x \ge 1 \f$ .
   
     The key consideration here is that the  \f$  b_i \f$  coefficients,  \f$ \ell' \f$ , and
      \f$ u' \f$  are functions of  \f$ \alpha \f$ , which determines which, among (1)-(4),
     to apply. In general,
  
     if  \f$ b_i > 0 \f$  then

      \f$ x^L_i \ge (l' - \sum_{j \in N_+'} b_j x^U_j 
     - \sum_{j \in N_-'} b_j x^L_j) / b_i \f$ ,
   						      	    
      \f$ x^U_i \le (u' - \sum_{j \in N_+'} b_j x^L_j
     - \sum_{j \in N_-'} b_j x^U_j) / b_i \f$ ;
   						      	    
     if  \f$ b_i < 0 \f$  then

      \f$ x^L_i \ge (l' - \sum_{j \in N_+'} b_j x^U_j
     - \sum_{j \in N_-'} b_j x^L_j) / b_i \f$ ,
   						      	    
      \f$ x^U_i \le (u' - \sum_{j \in N_+'} b_j x^L_j
     - \sum_{j \in N_-'} b_j x^U_j) / b_i \f$ .
   
     Each lower/upper bound is therefore a piecewise rational function
     of \f$ \alpha \f$ , given that \f$ b_i \f$ and the content of \f$
     N_+' \f$ and \f$ N_-' \f$ depend on \f$ \alpha \f$ . These
     functions are continuous (easy to prove) but not differentiable
     at some points of \f$ [0,1] \f$ .
  
     The purpose of this procedure is to find the maximum of the lower
     bounding function and the minimum of the upper bounding function.
  
     Divide the interval \f$ [0,1] \f$ into at most \f$ m+1 \f$
     intervals (where \f$ m \f$ is the number of coefficients not
     identically zero, or the number of \f$ b_i \f$ that are nonzero
     for at least one value of \f$ \alpha \f$ ). The limits \f$ c_i
     \f$ of the subintervals are the zeros of each coefficient,
     i.e. the values of \f$ \alpha \f$ such that \f$ \alpha a_{ki} +
     (1-\alpha) a_{hi} = 0 \f$ , or \f$ c_i = \frac{-a_{hi}}{a_{ki} -
     a_{hi}} \f$ .
  
     Sorting these values gives us something to do on every interval
      \f$ [c_j, c_{j+1}] \f$ when computing a new value of \f$ x^L_i
      \f$ and \f$ x^U_i \f$ , which I'll denote \f$ L_i \f$ and \f$
      U_i \f$ in the following.
  
     0) if  \f$ c_j = c_i \f$  then 
     - compute \f$ VL =  \lim_{c_j \to \alpha} L_i (\alpha) \f$ 
     - if  \f$ =  +\infty  \f$ , infeasible 
     else compute derivative DL (should be   \f$ +\infty \f$  )
  
     1) else 

     - compute \f$ VL = \lim_{\alpha \to c_j} L_i (\alpha) \f$ (can be
     retrieved from previous interval as \f$ L_i (\alpha) \f$ is
     continuous)

     - compute  \f$ DL = \lim_{\alpha \to c_j} dL_i (\alpha)  \f$ 
  
     update \f$ x^L \f$  with VL if necessary.
  
     2) if  \f$ c_{j+1} = c_i \f$  then 

     - compute \f$ VR = \lim_{\alpha \to c_{j+1}} L_i (\alpha) \f$ 

     - if =  \f$ +\infty  \f$  , infeasible else compute derivative DR (should be
       \f$ -\infty \f$  )
  
     3) else 
     - compute  \f$ VR = \lim_{\alpha \to c_{j+1}} L_i (\alpha)  \f$ 
     - compute  \f$ DR = \lim_{\alpha \to c_{j+1}} dL_i (\alpha)  \f$ 
  
     update  \f$  x^L \f$  with VR if necessary.
  
     if \f$ DL > 0 \f$ and \f$ DR < 0 \f$ , there might be a maximum
     in between, otherwise continue to next interval
  
     compute internal maximum VI, update \f$ x^L \f$ with VI if
     necessary.
  

     Apply a similar procedure for the upper bound

     This should be applied for any \f$ h,k,i \f$ , therefore we might
     have a lot to do. First, select possible pairs \f$ (h,k) \f$
     among those for which there exists at least one variable that
     satisfies neither of the following conditions:
  
     a) same sign coefficient, constraints \f$ (h,k) \f$ both \f$ \ge
     \f$ or both \f$ \le \f$
  
     b) opposite sign coefficient, constraints \f$ (h,k) \f$ \f$
     (\le,\ge) \f$ or \f$ (\ge,\le) \f$
  
     as in those cases, no \f$ c_i \f$ would be in \f$ [0,1] \f$
  */

  class CouenneTwoImplied: public CglCutGenerator {

  public:

    /// constructor
    CouenneTwoImplied (CouenneProblem *,
		       JnlstPtr,
		       const Ipopt::SmartPtr <Ipopt::OptionsList>);

    /// copy constructor
    CouenneTwoImplied  (const CouenneTwoImplied &);

    /// destructor
    ~CouenneTwoImplied ();

    /// clone method (necessary for the abstract CglCutGenerator class)
    CouenneTwoImplied *clone () const
    {return new CouenneTwoImplied (*this);}

    /// the main CglCutGenerator
    void generateCuts (const OsiSolverInterface &, 
		       OsiCuts &, 
		       const CglTreeInfo = CglTreeInfo ()) const;

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  protected:

    /// pointer to problem data structure (used for post-BT)
    CouenneProblem *problem_;

    /// Journalist
    JnlstPtr jnlst_;

    /// maximum number of trials in every call
    int nMaxTrials_;

    /// Total CPU time spent separating cuts
    mutable double totalTime_;

    /// CPU time spent columning the row formulation
    mutable double totalInitTime_;

    /// first call indicator
    mutable bool firstCall_;

    /// Depth of the BB tree where to start decreasing chance of running this
    int depthLevelling_;

    /// Depth of the BB tree where stop separation
    int depthStopSeparate_;
  };
}

#endif
