/* $Id$
 *
 * Name:    CouenneProblem.hpp
 * Author:  Pietro Belotti, Lehigh University
 *          Andreas Waechter, IBM
 * Purpose: define the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-10.
 * This file is licensed under the Eclipse Public License (EPL)
 */

#ifndef COUENNE_PROBLEM_HPP
#define COUENNE_PROBLEM_HPP

#define FM_TRACE_OPTSOL
#define FM_CHECKNLP2

#include <vector>
#include <map>

#include "CouenneConfig.h"

#include "CouenneTypes.hpp"
#include "CouenneExpression.hpp"

#include "CouenneJournalist.hpp"
#include "CouenneDomain.hpp"

namespace Ipopt {
  template <class T> class SmartPtr;
  class OptionsList;
  class Journalist;
}

namespace Bonmin {
  class RegisteredOptions;
  class BabInfo;
  class OsiTMINLPInterface;
  class BabSetupBase;
}

struct ASL;
struct expr;

class CglTreeInfo;
class CbcModel;
class OsiObject;
class CoinWarmStart;

class Nauty;

  class Node{
    int index;
    double coeff;
    double lb;
    double ub;
    int color;
    int code;
    int sign;
  public:
    void node(int, double, double, double, int, int);
    inline void color_vertex (register int k) {color = k;}
    inline int get_index () const {return index;}
    inline double get_coeff () const {return coeff;}
    inline double get_lb () const {return lb;}
    inline double get_ub () const {return ub;}
    inline int get_color () const {return color;}
    inline int get_code () const {return code;}
    inline int get_sign () const {return sign;}
    inline void bounds(register double a, register double b){ lb = a; ub = b;}
  };

#define COUENNE_EPS_SYMM 1e-8

  struct myclass0 {
    inline bool operator() (register const Node &a, register const Node &b) {

      return ((              a.get_code  () <  b.get_code  ())                     ||
	      ((             a.get_code  () == b.get_code  ()                      &&
		((           a.get_coeff () <  b.get_coeff ()  - COUENNE_EPS_SYMM) ||
		 ((fabs     (a.get_coeff () -  b.get_coeff ()) < COUENNE_EPS_SYMM) &&
		  ((         a.get_lb    () <  b.get_lb    ()  - COUENNE_EPS_SYMM) ||
		   ((fabs   (a.get_lb    () -  b.get_lb    ()) < COUENNE_EPS_SYMM) &&
		    ((       a.get_ub    () <  b.get_ub    ()  - COUENNE_EPS_SYMM) ||
		     ((fabs (a.get_ub    () -  b.get_ub    ()) < COUENNE_EPS_SYMM) &&
		      ((     a.get_index () <  b.get_index ())))))))))));

    //   bool is_less = 0;
    //   if(a.get_code() < b.get_code() )
    // 	is_less = 1;
    //   else {
    // 	if(a.get_code() == b.get_code() ) {
    // 	  if(a.get_coeff() < b.get_coeff() )
    // 	    is_less = 1;
    // 	  else{
    // 	    if(a.get_coeff() ==  b.get_coeff() ) {
    // 	      if(a.get_lb() < b.get_lb())
    // 		is_less = 1;
    // 	      else{
    // 		if(a.get_lb() == b.get_lb()) {
    // 		  if(a.get_ub() < b.get_ub())
    // 		    is_less = 1;
    // 		  else{
    // 		    if(a.get_ub() == b.get_ub()) {		    
    // 		      if(a.get_index() < b.get_index())
    // 			is_less = 1;
    // 		    }
    // 		  }
    // 		}
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // return is_less;

    }
  } ;
    
      
  struct myclass {
    inline bool operator() (register const  Node &a, register const Node &b) {
      return (a.get_index() < b.get_index() );
    }
  };


namespace Couenne {

  enum TrilinDecompType {rAI, treeDecomp, bi_tri, tri_bi};

  class exprVar;
  class exprAux;
  class DepGraph;
  class CouenneObject;
  class CouenneCutGenerator;
  class quadElem;
  class LinMap;
  class QuadMap;
  class CouenneConstraint;
  class CouenneObjective;
  class GlobalCutOff;
  //  class JnlstPtr;
  //  class ConstJnlstPtr;
  class CouenneRecordBestSol;

  typedef Ipopt::SmartPtr<Ipopt::Journalist> JnlstPtr;
  typedef Ipopt::SmartPtr<const Ipopt::Journalist> ConstJnlstPtr;

  struct compExpr;

// default tolerance for checking feasibility (and integrality) of NLP solutions
const CouNumber feas_tolerance_default = 1e-5;

/** Class for MINLP problems with symbolic information
 *
 *  It is read from an AMPL .nl file and contains variables, AMPL's
 *  "defined variables" (aka common expressions), objective(s), and
 *  constraints in the form of expression's. Changes throughout the
 *  program occur in standardization.
 */

class CouenneProblem {

  friend class exprMul;

  /// structure to record fixed, non-fixed, and continuous variables
  enum fixType {UNFIXED, FIXED, CONTINUOUS};

 public:

  /// Type of multilinear separation
  enum multiSep {MulSepNone, MulSepSimple, MulSepTight};

  // min depth for strong branching output
  int minDepthPrint_;

  // min number of nodes for strong branching output
  int minNodePrint_;

  // indicate if strong branching output should be printed
  bool doPrint_;

 protected:

  /// problem name
  std::string problemName_;

  std::vector <exprVar           *> variables_;   ///< Variables (original, auxiliary, and defined)
  std::vector <CouenneObjective  *> objectives_;  ///< Objectives
  std::vector <CouenneConstraint *> constraints_; ///< Constraints

  /// AMPL's common expressions (read from AMPL through structures cexps and cexps1)
  std::vector <expression *> commonexprs_; 

  mutable Domain domain_; ///< current point and bounds;

  /// Expression map for comparison in standardization and to count
  /// occurrences of an auxiliary
  std::set <exprAux *, compExpr> *auxSet_;

  /// Number of elements in the x_, lb_, ub_ arrays
  mutable int curnvars_;

  /// Number of discrete variables
  int nIntVars_;

  /// Best solution known to be loaded from file -- for testing purposes
  mutable CouNumber *optimum_;

  /// Best known objective function
  CouNumber bestObj_;

  /// Indices of variables appearing in products (used for SDP cuts)
  int *quadIndex_;

  /// Variables that have commuted to auxiliary
  bool *commuted_;

  /// numbering of variables. No variable xi with associated pi(i)
  /// greater than pi(j) should be evaluated before variable xj
  int *numbering_;

  /// Number of "defined variables" (aka "common expressions")
  int ndefined_;

  /// Dependence (acyclic) graph: shows dependence of all auxiliary
  /// variables on one another and on original variables. Used to
  /// create a numbering of all variables for evaluation and bound
  /// tightening (reverse order for implied bounds)
  DepGraph *graph_;

  /// Number of original variables
  int nOrigVars_;

  /// Number of original constraints (disregarding those that turned
  /// into auxiliary variable definition)
  int nOrigCons_;

  /// Number of original integer variables
  int nOrigIntVars_;

  /// Pointer to a global cutoff object
  mutable GlobalCutOff* pcutoff_;

  /// flag indicating if this class is creator of global cutoff object
  mutable bool created_pcutoff_;

  bool doFBBT_;  ///< do Feasibility-based bound tightening
  bool doRCBT_;  ///< do reduced cost      bound tightening
  bool doOBBT_;  ///< do Optimality-based  bound tightening
  bool doABT_;   ///< do Aggressive        bound tightening

  int logObbtLev_;   ///< frequency of Optimality-based bound tightening
  int logAbtLev_;    ///< frequency of Aggressive       bound tightening

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// window around known optimum (for testing purposes)
  CouNumber opt_window_;

  /// Use quadratic expressions?
  bool useQuadratic_;

  /// feasibility tolerance (to be used in checkNLP)
  CouNumber feas_tolerance_;

  /// inverse dependence structure: for each variable x give set of
  /// auxiliary variables (or better, their indices) whose expression
  /// depends on x
  std::vector <std::set <int> > dependence_;

  /// vector of pointer to CouenneObjects. Used by CouenneVarObjects
  /// when finding all objects related to (having as argument) a
  /// single variable
  std::vector <CouenneObject *> objects_;

  /// each element is true if variable is integer and, if auxiliary,
  /// depends on no integer
  mutable int *integerRank_;

  /// numberInRank_ [i] is the number of integer variables in rank i
  mutable std::vector <int> numberInRank_;

  /// maximum cpu time
  double maxCpuTime_;

  /// options
  Bonmin::BabSetupBase *bonBase_;

  /// AMPL structure pointer (temporary --- looking forward to embedding into OS...)
  ASL *asl_;

  /// some originals may be unused due to their zero multiplicity
  /// (that happens when they are duplicates). This array keeps track
  /// of their indices and is sorted by evaluation order
  int *unusedOriginalsIndices_;

  /// number of unused originals
  int nUnusedOriginals_;

  // to speedup sorting operations in strong branching
  int lastPrioSort_;

  // to record best solution found
  CouenneRecordBestSol *recBSol;

  /// Type of Multilinear separation
  enum multiSep multilinSep_;

  /// number of FBBT iterations
  int max_fbbt_iter_;

  /// true if FBBT exited for iteration limits as opposed to inability
  /// to further tighten bounds
  mutable bool fbbtReachedIterLimit_;

  /// use orbital branching?
  bool orbitalBranching_;

  /// check bounds on auxiliary variables when verifying MINLP
  /// feasibility of a solution. Usually this is not needed, unless
  /// some manipulation on auxiliary variables is done before
  /// Branch-and-Bound
  bool checkAuxBounds_;

  /// return type of decomposition of quadrilinear terms    
  enum TrilinDecompType trilinDecompType_;

  /// constant value of the objective if no variable is declared in it
  double constObjVal_;

 public:

  CouenneProblem  (ASL * = NULL,
		   Bonmin::BabSetupBase *base = NULL,
		   JnlstPtr jnlst = NULL);  ///< Constructor
  CouenneProblem  (const CouenneProblem &); ///< Copy constructor
  ~CouenneProblem ();                       ///< Destructor

  /// initializes parameters like doOBBT
  void initOptions (Ipopt::SmartPtr <Ipopt::OptionsList> options);

  /// Clone method (for use within CouenneCutGenerator::clone)
  CouenneProblem *clone () const
  {return new CouenneProblem (*this);}

  int nObjs     () const {return (int) objectives_.   size ();} ///< Get number of objectives
  int nCons     () const {return (int) constraints_.  size ();} ///< Get number of constraints
  int nOrigCons () const {return nOrigCons_;}                   ///< Get number of original constraints

  inline int nOrigVars    () const {return nOrigVars_;}                ///< Number of orig. variables
  inline int nDefVars     () const {return ndefined_;}                 ///< Number of def'd variables
  inline int nOrigIntVars () const {return nOrigIntVars_;}             ///< Number of original integers
  inline int nIntVars     () const {return nIntVars_;}                 ///< Number of integer variables
  inline int nVars        () const {return (int) variables_. size ();} ///< Total number of variables
  
  void setNDefVars (int ndefined__) {ndefined_ = ndefined__;}

  // Symmetry Info

  std::vector<int>  *Find_Orbit(int) const;
  mutable std::vector<Node> node_info;
  mutable Nauty *nauty_info;

  myclass0  node_sort; 
  myclass index_sort;

  void sym_setup();
  void Compute_Symmetry() const;
  void Print_Orbits() const;
  void ChangeBounds (const double * , const double *, int ) const;
  inline bool compare (register Node &a, register Node &b) const;
  Nauty *getNtyInfo () {return nauty_info;}

  // bool node_sort (  Node  a, Node  b);
  // bool index_sort (  Node  a, Node  b);

  /// empty if no NTY, symmetry data structure setup otherwise
  void setupSymmetry ();
  
  /// get evaluation order index 
  inline int evalOrder (int i) const
  {return numbering_ [i];}

  /// get evaluation order vector (numbering_)
  inline int *evalVector ()
  {return numbering_;}

  // get elements from vectors
  inline CouenneConstraint *Con (int i) const {return constraints_ [i];} ///< i-th constraint
  inline CouenneObjective  *Obj (int i) const {return objectives_  [i];} ///< i-th objective

  /// Return pointer to i-th variable
  inline exprVar *Var   (int i) const 
  {return variables_ [i];}

  /// Return vector of variables (symbolic representation)
  inline std::vector <exprVar *> &Variables () 
  {return variables_;}

  /// Return pointer to set for comparisons
  inline std::set <exprAux *, compExpr> *& AuxSet () 
  {return auxSet_;}

  /// Return pointer to dependence graph
  inline DepGraph *getDepGraph () 
  {return graph_;}

  /// return current point & bounds
  inline Domain *domain () const
  {return &domain_;}
  
  inline std::vector <expression *>& commonExprs() { return commonexprs_; }

  // Get and set current variable and bounds
  inline CouNumber   &X     (int i) const {return domain_.x   (i);} ///< \f$x_i\f$
  inline CouNumber   &Lb    (int i) const {return domain_.lb  (i);} ///< lower bound on \f$x_i\f$
  inline CouNumber   &Ub    (int i) const {return domain_.ub  (i);} ///< upper bound on \f$x_i\f$

  // get and set current variable and bounds
  inline CouNumber  *X     () const {return domain_.x  ();} ///< Return vector of variables
  inline CouNumber  *Lb    () const {return domain_.lb ();} ///< Return vector of lower bounds
  inline CouNumber  *Ub    () const {return domain_.ub ();} ///< Return vector of upper bounds

  // get optimal solution and objective value
  CouNumber  *&bestSol () const {return optimum_;} ///< Best known solution (read from file)
  CouNumber    bestObj () const {return bestObj_;} ///< Objective of best known solution

  /// Get vector of commuted variables
  bool *&Commuted () 
  {return commuted_;}

  /// Add (non linear) objective function
  void addObjective     (expression *, const std::string & = "min");

  // Add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression * = NULL); ///< Add equality constraint \f$ h(x) = b\f$
  void addGEConstraint  (expression *, expression * = NULL); ///< Add \f$\ge\f$ constraint, \f$h(x)\ge b\f$
  void addLEConstraint  (expression *, expression * = NULL); ///< Add \f$\le\f$ constraint, \f$h(x)\le b\f$
  void addRNGConstraint (expression *, expression * = NULL, 
			               expression * = NULL); ///< Add range constraint, \f$a\le h(x)\le b\f$

  /// Add (non linear) objective function
  void setObjective (int indObj = 0, expression * = NULL, const std::string & = "min");

  /// Add original variable.
  ///
  /// @param isint if true, this variable is integer, otherwise it is
  /// continuous
  expression *addVariable (bool isint = false, Domain *d = NULL);

  /// Add auxiliary variable and associate it with expression given as
  /// argument (used in standardization)
  exprAux *addAuxiliary (expression *);

  /// preprocess problem in order to extract linear relaxations etc.
  void reformulate (CouenneCutGenerator * = NULL);

  /// Break problem's nonlinear constraints in simple expressions to
  /// be convexified later. Return true if problem looks feasible,
  /// false if proven infeasible.
  bool standardize ();

  /// Display current representation of problem: objective, linear and
  /// nonlinear constraints, and auxiliary variables.
  void print (std::ostream & = std::cout);

#ifdef COIN_HAS_ASL
  /// Read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL *);

  /// Generate a Couenne expression from an ASL expression
  expression *nl2e (struct expr *, const ASL *asl);
#endif

  // bound tightening parameters
  bool doFBBT () const {return doFBBT_;} ///< shall we do Feasibility Based Bound Tightening?
  bool doRCBT () const {return doRCBT_;} ///< shall we do reduced cost      Bound Tightening?
  bool doOBBT () const {return doOBBT_;} ///< shall we do Optimality  Based Bound Tightening?
  bool doABT  () const {return doABT_;}  ///< shall we do Aggressive        Bound Tightening?

  int  logObbtLev () const {return logObbtLev_;} ///< How often shall we do OBBT?
  int  logAbtLev  () const {return logAbtLev_;}  ///< How often shall we do ABT?

  /// Write nonlinear problem to a .mod file (with lots of defined
  /// variables)
  /// 
  /// @param fname Name of the .mod file to be written
  ///
  /// @param aux controls the use of auxiliaries. If true, a problem
  /// is written with auxiliary variables written with their
  /// associated expression, i.e. \f$w_i = h_i(x,y,w)\f$ and bounds
  /// \f$l_i \le w_i \le u_i\f$, while if false these constraints are
  /// written in the form \f$l_i \le h_i (x,y) \le u_i\f$.
  ///
  /// Note: if used before standardization, writes original AMPL formulation
  void writeAMPL (const std::string &fname, bool aux);

  /// Write nonlinear problem to a .gms file
  /// 
  /// @param fname Name of the .gams file to be written.
  void writeGAMS (const std::string &fname);

  /// Initialize auxiliary variables and their bounds from original
  /// variables
  //void initAuxs (const CouNumber *, const CouNumber *, const CouNumber *);
  void initAuxs () const;

  /// Get auxiliary variables from original variables
  void getAuxs (CouNumber *) const;

  /// tighten bounds using propagation, implied bounds and reduced costs
  bool boundTightening (t_chg_bounds *, 
			Bonmin::BabInfo * = NULL) const;

  /// core of the bound tightening procedure
  bool btCore (t_chg_bounds *chg_bds) const;

  /// Optimality Based Bound Tightening
  int obbt (const CouenneCutGenerator *cg,
	    const OsiSolverInterface &csi,
	    OsiCuts &cs,
	    const CglTreeInfo &info,
	    Bonmin::BabInfo * babInfo,
	    t_chg_bounds *chg_bds);

  /// aggressive bound tightening. Fake bounds in order to cut
  /// portions of the solution space by fathoming on
  /// bounds/infeasibility
  bool aggressiveBT (Bonmin::OsiTMINLPInterface *nlp,
		     t_chg_bounds *, 
		     const CglTreeInfo &info,
		     Bonmin::BabInfo * = NULL) const;

  /// procedure to strengthen variable bounds. Return false if problem
  /// turns out to be infeasible with given bounds, true otherwise.
  int redCostBT (const OsiSolverInterface *psi,
		 t_chg_bounds *chg_bds) const;

  /// "Forward" bound tightening, that is, propagate bound of variable
  /// \f$x\f$ in an expression \f$w = f(x)\f$ to the bounds of \f$w\f$.
  int tightenBounds (t_chg_bounds *) const;

  /// "Backward" bound tightening, aka implied bounds. 
  int impliedBounds (t_chg_bounds *) const;

  /// Look for quadratic terms to be used with SDP cuts
  void fillQuadIndices ();

  /// Fill vector with coefficients of objective function
  void fillObjCoeff (double *&);

  /// Replace all occurrences of original variable with new aux given
  /// as argument
  void auxiliarize (exprVar *, exprVar * = NULL);

  /// Set cutoff 
  void setCutOff (CouNumber cutoff, const CouNumber *sol = NULL) const;

  /// Reset cutoff
  void resetCutOff (CouNumber value = COUENNE_INFINITY) const;

  /// Get cutoff
  CouNumber getCutOff () const;

  /// Get cutoff solution
  CouNumber *getCutOffSol () const;

  /// Make cutoff known to the problem
  void installCutOff () const;

  /// Provide Journalist
  ConstJnlstPtr Jnlst () const;

  /// Check if solution is MINLP feasible
  bool checkNLP (const double *solution, double &obj, bool recompute = false) const;

  /// generate integer NLP point Y starting from fractional solution
  /// using bound tightening
  int getIntegerCandidate (const double *xFrac, double *xInt, double *lb, double *ub) const;

  /// Read best known solution from file given in argument
  bool readOptimum (std::string *fname = NULL);

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  /// standardization of linear exprOp's
  exprAux *linStandardize (bool addAux, 
			   CouNumber c0, 
			   LinMap  &lmap,
			   QuadMap &qmap);

  /// split a constraint w - f(x) = c into w's index (it is returned)
  /// and rest = f(x) + c
  int splitAux (CouNumber, expression *, expression *&, bool *, enum expression::auxSign &);

  /// translates pair (indices, coefficients) into vector with pointers to variables
  void indcoe2vector (int *indexL,
		      CouNumber *coeff,
		      std::vector <std::pair <exprVar *, CouNumber> > &lcoeff);

  /// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
  void indcoe2vector (int *indexI,
		      int *indexJ,
		      CouNumber *coeff,
		      std::vector <quadElem> &qcoeff);

  /// given (expression *) element of sum, returns (coe,ind0,ind1)
  /// depending on element:
  ///
  /// 1) a * x_i ^ 2   ---> (a,i,?)   return COU_EXPRPOW
  /// 2) a * x_i       ---> (a,i,?)   return COU_EXPRVAR
  /// 3) a * x_i * x_j ---> (a,i,j)   return COU_EXPRMUL
  /// 4) a             ---> (a,?,?)   return COU_EXPRCONST
  ///
  /// x_i and/or x_j may come from standardizing other (linear or
  /// quadratic operator) sub-expressions
  void decomposeTerm (expression *term,
		      CouNumber initCoe,
		      CouNumber &c0,
		      LinMap  &lmap,
		      QuadMap &qmap);

  /// return problem name
  const std::string &problemName () const
  {return problemName_;}
  
  void setProblemName(std::string& problemName__)
  { problemName_ = problemName__; }

  /// return inverse dependence structure
  const std::vector <std::set <int> > &Dependence () const
  {return dependence_;}

  /// return object vector
  const std::vector <CouenneObject *> &Objects () const
  {return objects_;}

  /// find SOS constraints in problem
  int findSOS (CbcModel *CbcModelPtr,
	       OsiSolverInterface *solver, 
	       OsiObject ** objects);

  /// set maximum CPU time
  inline void setMaxCpuTime (double time)
  {maxCpuTime_ = time;}

  /// return maximum CPU time
  inline double getMaxCpuTime () const
  {return maxCpuTime_;}

  /// save CouenneBase
  void setBase (Bonmin::BabSetupBase *base);

  /// Some originals may be unused due to their zero multiplicity
  /// (that happens when they are duplicates). This procedure creates
  /// a structure for quickly checking and restoring their value after
  /// solving.
  void createUnusedOriginals ();

  /// Some originals may be unused due to their zero multiplicity (that
  /// happens when they are duplicates). This procedure restores their
  /// value after solving
  void restoreUnusedOriginals (CouNumber * = NULL) const;

  /// return indices of neglected redundant variables
  int *unusedOriginalsIndices () 
  {return unusedOriginalsIndices_;}

  /// number of unused originals
  int nUnusedOriginals ()
  {return nUnusedOriginals_;}

  /// return type of separator for multilinear terms
  enum multiSep MultilinSep () const
  {return multilinSep_;}

  /// true if latest call to FBBT terminated due to iteration limit reached
  bool fbbtReachedIterLimit () const
  {return fbbtReachedIterLimit_;}

  /// return true if orbital branching activated
  bool orbitalBranching () const
  {return orbitalBranching_;}

  /// set the value for checkAuxBounds. When true, all MINLP feasible
  /// solutions will additionally be tested for feasibility with
  /// respect to auxiliary variable bounds. This is normally not needed.
  void setCheckAuxBounds (bool value) 
  {checkAuxBounds_ = value;}

  /// return true if bounds of auxiliary variables have to be satisfied
  /// whenever a solution is tested for MINLP feasibiliry
  bool checkAuxBounds () const
  {return checkAuxBounds_;}

  /// return type of decomposition of quadrilinear terms    
  enum TrilinDecompType getTrilinDecompType ()
  {return trilinDecompType_;}

  /// options
  Bonmin::BabSetupBase *bonBase () const {return bonBase_;}

  /// returns constant objective value if it contains no variables
  inline double constObjVal () const {return constObjVal_;}

protected:

  /// single fake tightening. Return
  ///
  /// -1   if infeasible
  ///  0   if no improvement
  /// +1   if improved
  int fake_tighten (char direction,  ///< 0: left, 1: right
		    int index,       ///< index of the variable tested
		    const double *X, ///< point round which tightening is done
		    CouNumber *olb,  ///< cur. lower bound
		    CouNumber *oub,  ///< cur. upper bound
		    t_chg_bounds *chg_bds,
		    t_chg_bounds *f_chg) const;

  /// Optimality Based Bound Tightening -- inner loop
  int obbtInner (OsiSolverInterface *, 
		 OsiCuts &,
		 t_chg_bounds *, 
		 Bonmin::BabInfo *) const;

  int obbt_iter (OsiSolverInterface *csi, 
		 t_chg_bounds *chg_bds, 
		 const CoinWarmStart *warmstart, 
		 Bonmin::BabInfo *babInfo,
		 double *objcoe,
		 int sense, 
		 int index) const;

  int call_iter (OsiSolverInterface *csi, 
		 t_chg_bounds *chg_bds, 
		 const CoinWarmStart *warmstart, 
		 Bonmin::BabInfo *babInfo,
		 double *objcoe,
		 enum nodeType type,
		 int sense) const;

  /// analyze sparsity of potential exprQuad/exprGroup and change
  /// linear/quadratic maps accordingly, if necessary by adding new
  /// auxiliary variables and including them in the linear map
  void analyzeSparsity (CouNumber, 
			LinMap &,
			QuadMap &);

  /// re-organizes multiplication and stores indices (and exponents) of
  /// its variables
  void flattenMul (expression *mul, 
		   CouNumber &coe, 
		   std::map <int, CouNumber> &indices);

  /// clear all spurious variables pointers not referring to the variables_ vector
  void realign ();

  /// fill dependence_ structure
  void fillDependence (Bonmin::BabSetupBase *base, CouenneCutGenerator * = NULL);

  /// fill freeIntegers_ array
  void fillIntegerRank () const;

  /// Test fixing of an integer variable (used in getIntegerCandidate())
  int testIntFix (int index, 
		  CouNumber xFrac, 
		  enum fixType *fixed,
		  CouNumber *xInt,
		  CouNumber *dualL, CouNumber *dualR,
		  CouNumber *olb,   CouNumber *oub,
		  bool patient) const;

public:

  /// 
  inline int getLastPrioSort() const 
  {return lastPrioSort_;}

  ///
  void setLastPrioSort (int givenLastPS);

  /// returns recorded best solution
  inline CouenneRecordBestSol *getRecordBestSol() const 
  {return recBSol;}

  /// returns feasibility tolerance
  double getFeasTol() {return feas_tolerance_;}

  /// Recompute objective value for sol
  double checkObj(const CouNumber *sol, const double &precision) const;

  /// check integrality of vars in sol with index between from and upto 
  /// (original vars only if origVarOnly == true); 
  /// return true if all integer vars are within precision of an integer value
  bool checkInt(const CouNumber *sol,
		const int from, const int upto, 
		const std::vector<int> listInt,
		const bool origVarOnly, 
		const bool stopAtFirstViol,
		const double precision, double &maxViol) const;

  /// Check bounds; returns true iff feasible for given precision
  bool checkBounds(const CouNumber *sol,
		   const bool stopAtFirstViol,
		   const double precision, double &maxViol) const;

  /// returns true iff value of all auxilliaries are within bounds
  bool checkAux(const CouNumber *sol,
		const bool stopAtFirstViol,
		const double precision, double &maxViol) const;

  /// returns true iff value of all auxilliaries are within bounds
  bool checkCons(const CouNumber *sol,
		 const bool stopAtFirstViol,
		 const double precision, double &maxViol) const;

  /// Return true if either solution or recomputed_solution obtained
  /// using getAuxs() from the original variables in solution is feasible
  /// within precision (the solution with minimum violation is then stored
  /// in recBSol->modSol, as well as its value and violation); 
  /// return false otherwise.
  /// If stopAtFirstViol == true, recBSol->modSol is meaningless upon return.
  /// If stopAtFirstViol == false, recBSol->modSol contains the solution
  /// with minimum violation, although this violation might be larger than 
  /// precision.
  /// This is useful for cases where the current solution must be considered
  /// valid (e.g., because Cbc is going to accept it anyway), although it 
  /// violates precision requirements.

  /// Value of obj matters only if careAboutObj == true;
  /// the code then tries to balance violation of constraints and
  /// value of objective.

  /// if checkAll = false, check only integrality/bounds for 
  /// original vars and constraints; consider only recomputed_sol
  /// if checkAll == true, check also integrality/bounds on auxs;
  /// consider both recomputed_sol and solution

  /// if careAboutObj is set to true, then stopAtFirstViol must be set to 
  /// false too.
  bool checkNLP2(const double *solution,
		 const double obj, const bool careAboutObj,
		 const bool stopAtFirstViol,
		 const bool checkAll,
		 const double precision) const;
};

}

#endif
