/* $Id$ 
 *
 * Name:    Nauty.cpp
 * Authors: Jim Ostrowski
 * Purpose: Branching with symmetry -- implementation of the Nauty object
 * Date:    October 13, 2010
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <cassert>
#include <cmath>
#include "Nauty.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include "CoinTime.hpp"
//#include "OrbitalOptions.h"
//extern OrbitalOptions *options;

int Nauty::nautyCalls_ = 0;
double Nauty::nautyTime_ = 0.0;

Nauty::Nauty(int vertices)
{

  n_ = vertices;
  m_ = (n_ + WORDSIZE - 1)/WORDSIZE;

  //printf ("size of long = %d (%d)\nwordsize = %d\nn,m = %d,%d\n", 
  //          SIZEOF_LONG, sizeof (long), WORDSIZE, n_, m_);

  nauty_check (WORDSIZE, m_, n_, NAUTYVERSIONID);

  /// Apparently sizes are skewed on 64bit machines

#define MULTIPLIER 2

  G_ = (graph *) malloc(MULTIPLIER * m_ * n_ * sizeof(int));
  lab_ = (int *) malloc(MULTIPLIER * n_ * sizeof(int));  
  ptn_ = (int *) malloc(MULTIPLIER * n_ * sizeof(int));
  active_ = NULL;
  orbits_ = (int *) malloc(MULTIPLIER * n_ * sizeof(int));
  options_ = (optionblk *) malloc(MULTIPLIER * sizeof(optionblk));
  stats_ = (statsblk *) malloc(MULTIPLIER * sizeof(statsblk));
  worksize_ = 100*m_;
  workspace_ = (setword *) malloc(MULTIPLIER * worksize_*sizeof(setword));
  canonG_ = NULL;
  if (G_ == 0 || lab_ == 0 || ptn_ == 0 || 
      orbits_ == 0 || options_ == 0 || stats_ == 0 ||
      workspace_ == 0) assert(0);

  // Zero allocated memory
  memset(G_, 0, m_*n_*sizeof(int));
  memset(lab_, 0, n_*sizeof(int));
  memset(ptn_, 0, n_*sizeof(int));
  memset(orbits_, 0, n_*sizeof(int));
  memset(workspace_, 0, worksize_*sizeof(setword));

  // Set the options you want
  options_->getcanon = FALSE;
  options_->digraph = FALSE;
  options_->writeautoms = FALSE;
  options_->writemarkers = FALSE;
  options_->defaultptn = TRUE;
  options_->cartesian = FALSE;
  options_->linelength = 78;
  options_->outfile = NULL;
  options_->userrefproc = NULL;
  options_->userautomproc = NULL;
  options_->userlevelproc = NULL;
  options_->usernodeproc = NULL;
  //  options_->usertcellproc = NULL;
  options_->invarproc = NULL;
  options_->tc_level = 100;
  options_->mininvarlevel = 0;
  options_->maxinvarlevel = 1;
  options_->invararg = 0;
  options_->dispatch = &dispatch_graph;
   // Make an empty graph
  for (int j = 0; j < n_; j++) {
    set *gv = GRAPHROW(G_, j, m_);
    EMPTYSET(gv, m_);
  }

  vstat_ = new int[n_];
   clearPartitions();
   afp_ = NULL;
 }

Nauty::~Nauty()
{
  if (G_) free(G_);
  if (lab_) free(lab_);
  if (ptn_) free(ptn_);
  if (active_) free(active_);
  if (orbits_) free(orbits_);
  if (options_) free(options_);
  if (stats_) free(stats_);
  if (workspace_) free(workspace_);
  if (canonG_) free(canonG_);

  if (vstat_) delete [] vstat_;
}

void
Nauty::addElement(int ix, int jx)
{
  // Right now die if bad index.  Can throw exception later
  //printf("addelement %d %d \n", ix, jx);
  assert(ix < n_ and jx < n_);
  if(ix != jx){  //No Loops
    set *gv = GRAPHROW(G_, ix, m_);
    ADDELEMENT(gv, jx);
    set *gv2 = GRAPHROW(G_, jx, m_);
    ADDELEMENT(gv2, ix);
    autoComputed_ = false;
  }
}

void 
Nauty::clearPartitions()
{
  for (int j = 0; j < n_; j++) {
    vstat_[j] = 1;
    //printf("vstat %d = %d", j, vstat_[j]);
  }

  autoComputed_ = false;
}

void
Nauty::computeAuto()
{

  //  if (autoComputed_) return;

  double startCPU = CoinCpuTime ();

  options_->defaultptn = FALSE;

  // Here we only implement the partitions
  // [ fix1 | fix0 (union) free | constraints ]
  int ix = 0;
  
  for( int color = 1; color <= n_; color++){
    for (int j = 0; j < n_; j++) {
      if (vstat_[j] == color) {
        lab_[ix] = j;
        ptn_[ix] = color;
        ix++;
      }
    }
     if (ix > 0) ptn_[ix-1] = 0;
  }
  
  /*
  for (int j = 0; j < n_; j++)
    printf("ptn %d = %d      lab = %d \n", j, ptn_[j], lab_[j]);
  */
  

  // Should be number of columns
  assert(ix == n_);
  // Now the constraints if needed


    // Compute Partition
    
  nauty(G_, lab_, ptn_, active_, orbits_, options_, 
        stats_, workspace_, worksize_, m_, n_, canonG_);
  autoComputed_ = true;

  double endCPU = CoinCpuTime ();

  nautyCalls_++;
  nautyTime_ += endCPU - startCPU;
  // Need to make sure all generators are written
  if (afp_) fflush(afp_);
   
}

void
Nauty::deleteElement(int ix, int jx)
{
  // Right now die if bad index.  Can throw exception later
  assert(ix < n_ and jx < n_);
  set *gv = GRAPHROW(G_, ix, m_);
  if (ISELEMENT(gv, jx)) {
    DELELEMENT(gv, jx);
  } 
  autoComputed_ = false;
}

double
Nauty::getGroupSize() const
{
  if (!autoComputed_) return -1.0;
  return( stats_->grpsize1 * pow(10.0, (double) stats_->grpsize2) );
}

int
Nauty::getNumGenerators() const
{
  if (!autoComputed_) return -1;
  return(stats_->numgenerators);
}

int
Nauty::getNumOrbits() const
{
  if (!autoComputed_) return -1;
  return(stats_->numorbits);
}

std::vector<std::vector<int> >
*Nauty::getOrbits() const
{
  std::vector<std::vector<int> > *orb = new std::vector<std::vector<int> >;
  if (!autoComputed_) return orb;
  orb -> resize(getNumOrbits());
  std::multimap<int, int> orbmap;
  std::set<int> orbkeys;
  for (int j = 0; j < n_; j++) {
    orbkeys.insert(orbits_[j]);
    orbmap.insert(std::make_pair(orbits_[j], j));
  }

  int orbix = 0;
  for (std::set<int>::iterator it = orbkeys.begin();
       it != orbkeys.end(); ++it) {
    std::multimap<int, int>::iterator pos;
    for (pos = orbmap.lower_bound(*it);
         pos != orbmap.upper_bound(*it); ++pos) {
      (*orb)[orbix].push_back(pos->second);
    }
    orbix++;
  }

  assert(orbix == getNumOrbits());
  return orb;
}

void
Nauty::getVstat(double *v, int nv)
{
  assert(nv == n_);
  memcpy(v, vstat_, nv * sizeof(VarStatus));
}

/*
bool
Nauty::isAllFixOneOrbit(const std::vector<int> &orbit) const
{

  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;
    if (vstat_[*it] != FIX_AT_ONE) return false;
  }
  return true;
}

bool
Nauty::isAllFreeOrbit(const std::vector<int> &orbit) const
{
  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;
    if (vstat_[*it] != FREE) return false;
  }
  return true;
}

bool
Nauty::isConstraintOrbit(const std::vector<int> &orbit) const
{
  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return true;
  }
  return false;
  
}

bool
Nauty::isMixedFreeZeroOrbit(const std::vector<int> &orbit) const
{
  bool containsFree = false;
  bool containsZero = false;

  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;    
    if (vstat_[*it] == FREE) containsFree = true;
    if (vstat_[*it] == FIX_AT_ZERO) containsZero = true;    
    if (containsFree && containsZero) break;
  }  
  return (containsFree && containsZero);
}
*/
/*
void 
Nauty::setWriteAutoms(const std::string &fname)
{
  afp_ = fopen(fname.c_str(), "w");
  options_->writeautoms = TRUE;
  options_->writemarkers = FALSE;
  options_->outfile = afp_;

}

void 
Nauty::unsetWriteAutoms()
{
  fclose(afp_);
  options_->writeautoms = FALSE;
}
*/
