/* $Id$
 *
 * Name:    CouenneOrbitObj.hpp
 * Authors: Jim Ostrowski, University of Waterloo
 *          Pietro Belotti, Lehigh University
 * Purpose: Object for orbital branching
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */
/*
#ifndef COUENNEORBITOBJ_HPP
#define COUENNEORBITOBJ_HPP

#include "BonBabSetupBase.hpp"

#include "CouenneExprVar.hpp"
#include "CouenneJournalist.hpp"
#include "OsiBranchingObject.hpp"
#include "CouenneObject.hpp"
#include "Nauty.h"

namespace Couenne {

class Node{
  int index;
  double coeff;
  double lb;
  double ub;
  int color;
  int code;
public:
  void node(int, double, double, double, int);
  void color_vertex(int);
  int get_index () {return index; 
  };
  double get_coeff () {return coeff;      
  };
  double get_lb () {return lb;      
  };
  double get_ub () {return ub ;     
  };
  int get_color () {return color;     
  };
  int get_code () {return code;     
  };
  void bounds( double a, double b){ lb = a; ub = b;
  };
};

  bool compare (  Node a, Node b);
  bool node_sort (  Node a, Node b);
  bool index_sort (  Node a, Node b);


/// OsiObject for Orbital Branching
class CouenneOrbitObj: public CouenneObject {

public:

  /// empty constructor (for unused objects)
  CouenneOrbitObj ();

  /// Constructor with information for branching point selection strategy
  CouenneOrbitObj (CouenneCutGenerator *cutgen,
		   CouenneProblem *p, 
		   exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);


		   

  /// Constructor with lesser information, used for infeasibility only
  CouenneOrbitObj (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Destructor
  ~CouenneOrbitObj () {}

  /// Copy constructor
  CouenneOrbitObj (const CouenneOrbitObj &src);

  /// Cloning method
  virtual CouenneObject * clone () const
  {return new CouenneOrbitObj (*this);}


  
  /// set object parameters by reading from command line
  void setParameters (Bonmin::BabSetupBase *base);
  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// compute infeasibility of this variable, |w - f(x)|, where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double checkInfeasibility (const OsiBranchingInformation * info) const;

  /// fix (one of the) arguments of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject *createBranch (OsiSolverInterface*,const OsiBranchingInformation*, int) const;


  
  void Compute_Symmetry();
  void Print_Orbits();
  void ChangeBounds (CouenneProblem *p);
    

  std::vector<Node> node_info;
  
  Nauty *nauty_info;
  
protected:

};

}

#endif
*/
