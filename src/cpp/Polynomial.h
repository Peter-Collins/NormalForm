#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Polynomial
//
// PURPOSE:
//
// Implement multivariate polynomials as maps from multi-indices
// (Powers objects) to coefficients.
//
// NOTES:
//
// Of interest is the associated PolynomialRingInterface,
// PolynomialRing, and LieAlgebra modules, which provide more
// convenient ways to express some polynomial operations, for example,
// the factories below for coordinate monomials are more naturally
// invoked as ordinary (non-static) methods from the corresponding
// polynomial ring.
//
// Notice that coefficients are only stored if the boolean
// isValidCoeff is true.  By changing isValidCoeff to measure the size
// of a coefficient, rather than simply returning whether it is
// non-zero, we can make polynomials that truncate small terms
// automically.
//
//----------------------------------------------------------------------

//system headers
#include <assert.h>
#include <map>

//library headers

//project headers
#include "Types.h"
#include "Representations.h"

//type definitions
typedef std::map< Powers, Coeff > PowersToCoeffMap;
typedef bool (*PowersAndCoeffToBool)(const Powers& p, const Coeff& c);

//
// philosophy: throw exceptions for errors that most likely arose in
// client code; use assertions for internal errors
//
class PolynomialError {
};
class PolynomialIndexError : public PolynomialError {
};
class PolynomialSizeMismatchError : public PolynomialError {
};
class PolynomialValueError : public PolynomialError {
};

class Polynomial {

 public:

  //the big four (default, copy, assign, destruct)
  Polynomial(void);
  Polynomial(const Polynomial& that);
  Polynomial& operator=(const Polynomial& that);
  ~Polynomial(void);

  //constructors
  explicit Polynomial(const Index numVars);
  explicit Polynomial(const Powers& powers);
  Polynomial(const Powers& powers, const Coeff& coeff);

  //factories
  static Polynomial Zero(const Index numVars);
  static Polynomial One(const Index numVars);
  static Polynomial CoordinateMonomial(const Index numVars,
                                       const Index i,
                                       const Power p = 1);

  //basic inspectors and getters
  inline Index getNumVars(void) const;
  inline size_t getNumTerms(void) const;
  inline const PowersToCoeffMap& getPowersAndCoeffs(void) const;
  const Coeff& getCoeff(const Powers& powers) const;
  inline const Coeff& operator[](const Powers& powers) const;

  //basic Setters
  Polynomial& setMonomial(const Powers& powers, const Coeff& coeff);
  Polynomial& addMonomial(const Powers& powers, const Coeff& coeff);

  //boolean properties
  inline bool isZero(void) const;
  bool hasPowers(const Powers& powers) const;
  bool isConstant(void) const;
  bool isHomogeneous(void) const;
  bool isHomogeneous(const Power p) const;

  //comparators
  inline bool operator!=(const Polynomial& that) const;
  inline bool operator==(const Polynomial& that) const;

  //unary negation
  Polynomial operator-(void) const;

  //arithmetic assignments
  Polynomial& operator+=(const Polynomial& that);
  Polynomial& operator-=(const Polynomial& that);
  Polynomial& operator*=(const Polynomial& that);
  Polynomial& operator*=(const Coeff& c);

  //binary arithmetic operators (some are derived from arithmetic assigns)
  inline friend Polynomial operator+(const Polynomial& a, const Polynomial& b);
  inline friend Polynomial operator-(const Polynomial& a, const Polynomial& b);
  inline friend Polynomial operator*(const Polynomial& p, const Coeff& c);
  inline friend Polynomial operator*(const Coeff& c, const Polynomial& p);
  friend Polynomial operator*(const Polynomial& a, const Polynomial& b);

  //arithmetic methods
  Polynomial pow(const Power p) const;

  //evaluate and differentiate
  Coeff operator()(const std::vector< Coeff >& v) const;
  Polynomial operator()(const std::vector< Polynomial >& v) const;
  Polynomial diff(const Index i) const;

  //homogeneous structure and partitions
  Power degree(void) const;
  Polynomial homogeneousPart(const Power deg, const Power upper = 0) const;
  std::pair< Polynomial, Polynomial > partition(const PowersAndCoeffToBool& pred) const;

  //norms
  Real lInfinityNorm(void) const;
  Real lOneNorm(void) const;

  //checking
  inline bool isValidIndex(const Index i) const;
  inline bool hasSameNumVars(const Powers& powers) const;
  inline bool hasSameNumVars(const Polynomial& that) const;

  //pretty printing
  friend std::ostream& operator<<(std::ostream& ostr, const Polynomial& poly);

 protected:

  Index numVars_;
  PowersToCoeffMap terms_;

}; //Polynomial

inline Index Polynomial::getNumVars(void) const {
  return numVars_;
}
inline size_t Polynomial::getNumTerms(void) const {
  return terms_.size();
}
inline bool Polynomial::isValidIndex(const Index i) const {
  return ((i >= 0) && (i < numVars_));
}
inline bool Polynomial::hasSameNumVars(const Powers& powers) const {
  return (powers.numVars_ == numVars_);
}
inline bool Polynomial::hasSameNumVars(const Polynomial& that) const {
  return (that.numVars_ == numVars_);
}
inline const Coeff& Polynomial::operator[](const Powers& powers) const {
  return getCoeff(powers);
}
inline const PowersToCoeffMap& Polynomial::getPowersAndCoeffs(void) const {
  return terms_;
}
inline bool Polynomial::operator!=(const Polynomial& that) const {
  if (!(numVars_ == that.numVars_))
    throw PolynomialSizeMismatchError();
  if (that.terms_ != terms_)
    return true;
  return false;
}
inline bool Polynomial::operator==(const Polynomial& that) const {
  return !(operator!=(that));
}
inline bool Polynomial::isZero(void) const {
  return terms_.empty();
}
inline Polynomial operator+(const Polynomial& a, const Polynomial& b) {
  Polynomial result(a);
  return (result += b);
}
inline Polynomial operator-(const Polynomial& a, const Polynomial& b) {
  Polynomial result(a);
  return (result -= b);
}
inline Polynomial operator*(const Polynomial& p, const Coeff& c) {
  Polynomial result(p);
  return (result *= c);
}
inline Polynomial operator*(const Coeff& c, const Polynomial& p) {
  Polynomial result(p);
  return (result *= c);
}

#endif //POLYNOMIAL_H
