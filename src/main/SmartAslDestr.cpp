// $Id$
//
// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License (EPL).
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
// Pietro Belotti, Lehigh University
//
// Date : 05/31/2010

// Ampl includes

#include "BonCouenneSetup.hpp"

#ifdef COIN_HAS_ASL
#include "asl.h"
#include "getstub.h"
#endif

using namespace Couenne;
  
SmartAsl::~SmartAsl(){
#ifdef COIN_HAS_ASL
  //Code from Ipopt::AmplTNLP to free asl
  if(asl != NULL){
    if (X0) {
      delete [] X0;
      X0 = NULL;
    }
    if (havex0) {
      delete [] havex0;
      havex0 = NULL;
    }
    if (pi0) {
      delete [] pi0;
      pi0 = NULL;
    }
    if (havepi0) {
      delete [] havepi0;
      havepi0 = NULL;
    }
    ASL* asl_to_free = (ASL*)asl;
    ASL_free(&asl_to_free);
    asl = NULL;
  }
  ASL_free(&asl);
#endif
}
  
