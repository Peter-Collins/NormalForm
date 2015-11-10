#ifndef POLYNOMIAL_RING_INTERFACE_H
#define POLYNOMIAL_RING_INTERFACE_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialRingInterface
//
// PURPOSE:
//
// Abstract interface for all polynomial rings, including those that
// will be implemented by graded lie algebras.
//
// NOTES:
//
// To allow such flexibility, we defer knowledge of how many variables
// the ring has to a pure virtual method getNumVars_(), which we will
// override in derived classes.
//
//----------------------------------------------------------------------

#include "Types.h"
#include "Polynomial.h"

class PolynomialRingInterface {
 public:
  //delegated to pure virtual methods
  inline Index getNumVars(void) const;

  //implemented
  inline bool hasElt(const Polynomial& poly) const;
  inline Polynomial one(void) const;
  inline Polynomial zero(void) const;
  inline Polynomial monomial(const Powers& pows, const Coeff& coeff) const;
  inline Polynomial coordinateMonomial(const Index i, const Power p = 1) const;

 protected:
  //three of the big four not needed since we have a pure virtual method
  //PolynomialRingInterface(void);
  //PolynomialRingInterface(const PolynomialRingInterface& ring);
  //PolynomialRingInterface& operator=(const PolynomialRingInterface& ring);

  //virtual destructor; we must implement this
  virtual ~PolynomialRingInterface(void);

  //inspectors
  virtual Index getNumVars_(void) const = 0;

 private:
};

inline
Index PolynomialRingInterface::getNumVars(void) const {
  return getNumVars_();
}

inline
bool PolynomialRingInterface::hasElt(const Polynomial& poly) const {
  return (poly.getNumVars() == getNumVars_());
}

inline
Polynomial PolynomialRingInterface::one(void) const {
  return Polynomial::One(getNumVars_());
}

inline
Polynomial PolynomialRingInterface::zero(void) const {
  return Polynomial::Zero(getNumVars_());
}

inline
Polynomial PolynomialRingInterface::coordinateMonomial(const Index i,
                                                       const Power p) const {
  return Polynomial::CoordinateMonomial(getNumVars_(), i, p);
}

#endif //POLYNOMIAL_RING_INTERFACE_H

