#ifndef POLYNOMIAL_RING_H
#define POLYNOMIAL_RING_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialRing
//
// An implementation of PolynomialRingInterface with explicit numVars.
//
//----------------------------------------------------------------------

//project headers
#include "Polynomial.h"
#include "PolynomialRingInterface.h"

class PolynomialRingError {
};
class PolynomialRingIndexError : public PolynomialRingError {
};

class PolynomialRing : public PolynomialRingInterface {
 public:
  //the big four
  PolynomialRing(void);
  explicit PolynomialRing(const PolynomialRing& ring);
  ~PolynomialRing(void);
  //constructor from number of variables
  explicit PolynomialRing(const Index numVars);

 protected:
  //inspectors
  inline Index getNumVars_(void) const;
  //forbid assignment
  PolynomialRing& operator=(const PolynomialRing& ring);

 private:
  const Index numVars_;
};

Index PolynomialRing::getNumVars_(void) const {
  assert(numVars_ > 0);
  return numVars_;
}

#endif //POLYNOMIAL_RING_H
