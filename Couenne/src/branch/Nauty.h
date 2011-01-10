/* $Id$ 
 *
 * Name:    Nauty.cpp
 * Authors: Jim Ostrowski
 * Purpose: Branching with symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef NAUTY_H
#define NAUTY_H

extern "C" {
#include "nauty.h"
}

#include <cstdio>
#include <map>
#include <vector>

class Nauty
{

public:
  enum VarStatus { FIX_AT_ZERO, FIX_AT_ONE, FREE };
  
  Nauty(int n_);
  ~Nauty();

  void addElement(int ix, int jx);
  void clearPartitions();
  void computeAuto();
  void deleteElement(int ix, int jx);
  void color_node(int ix, int color) { vstat_[ix] = color; }
  void insertRHS(int rhs , int cons) {constr_rhs.insert( std::pair<int,int>(rhs,cons));}
  
  double getGroupSize() const;
  int getNautyCalls() const { return nautyCalls_; }
  double getNautyTime() const { return nautyTime_; }

  int getN() const { return n_; }
  
  int getNumGenerators() const;
  int getNumOrbits() const;

  /// Returns the orbits in a "convenient" form
  std::vector<std::vector<int> > *getOrbits() const;

  void getVstat(double *v, int nv);

  /**
   * Methods to classify orbits.  Not horribly efficient, but gets the job done
   */
  //  bool isAllFixOneOrbit(const std::vector<int> &orbit) const;
  // bool isAllFreeOrbit(const std::vector<int> &orbit) const;
  //bool isAutoComputed() const { return autoComputed_; }
  //bool isConstraintOrbit(const std::vector<int> &orbit) const;
  //bool isMixedFreeZeroOrbit(const std::vector<int> &orbit) const;
  //void makeFree(int ix) { vstat_[ix] = FREE; }  
  //  void setWriteAutoms(const std::string &afilename);
  //void unsetWriteAutoms();

private:

  Nauty ();

  // The base nauty stuff
  graph *G_;
  int *lab_;
  int *ptn_;
  set *active_;
  int *orbits_;
  optionblk *options_;
  statsblk *stats_;
  setword *workspace_;
  int worksize_;
  int m_;
  int n_;
  graph *canonG_;
  
  bool autoComputed_;

  int *vstat_;

  static int nautyCalls_;
  static double nautyTime_;

  std::multimap<int,int> constr_rhs;
  std::multimap<int,int>::iterator it;

  std::pair<std::multimap<int,int>::iterator,
            std::multimap<int,int>::iterator> ret;

  // File pointer for automorphism group
  FILE *afp_;

};

#endif
