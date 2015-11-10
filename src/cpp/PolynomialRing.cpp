//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialRing implementation
//
//----------------------------------------------------------------------

//project headers
#include "PolynomialRing.h"

PolynomialRing::PolynomialRing(void)
  : numVars_(1) {
}

PolynomialRing::PolynomialRing(const Index numVars)
  : numVars_(numVars) {
  if (!(numVars > 0))
    throw PolynomialRingIndexError();
}

PolynomialRing::PolynomialRing(const PolynomialRing& that)
  : numVars_(that.numVars_) {
}

PolynomialRing::~PolynomialRing(void) {
  //nothing
}
