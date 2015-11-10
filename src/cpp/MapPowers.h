#ifndef MAP_POWERS_H
#define MAP_POWERS_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: MapPowers
//
// PURPOSE:
//
// A product of pure (non-negative) powers, represented efficiently by
// using a map from coordinate indices to non-zero powers; especially
// suitable for large numbers of variables.  Multiplying two such
// objects (having the same number of variables) results in adding
// together the powers component-wise, for example.
//
// NOTES:
//
// On one level, these are multi-indices.
//
// We do NOT overload ::operator[] for lvalues, as this leads to some
// over-complication.  Instead, we provide a ::set(index, power)
// method.  Let's not worry too much about this lack of syntactic
// sugar; it doesn't look worth the effort to provide it.
//
//----------------------------------------------------------------------

//system headers
#include <vector>
#include <iostream>

//library headers

//project headers
#include "Types.h"

//philosophy: throw exceptions for errors due to client code; use
//assertions for internal errors
class MapPowersError {
};
class MapPowersIndexError : public MapPowersError {
};
class MapPowersSizeMismatchError : public MapPowersError {
};
class MapPowersValueError : public MapPowersError {
};
// gcc 4 needs this here for forward declaration of polynomial below
class Polynomial;

class MapPowers {
 public:
  //the big four
  MapPowers(void);
  MapPowers(const MapPowers& that);
  MapPowers& operator=(const MapPowers& that);
  ~MapPowers(void);

  //other constructors and factories
  explicit MapPowers(const Index len);
  //the following left implicit for convenience
  MapPowers(const std::vector< Power >& powers);
  MapPowers(const Index len, const IndexToPowerMap& powers);

  //inspectors
  inline Index getNumVars(void) const;
  inline Index getNumPowersStored(void) const;
  Power getPower(const Index i) const;
  inline Power operator[](const Index i) const;

  //item assignment
  MapPowers& setPower(const Index i, const Power p);

  //checking
  inline bool isValidIndex(const Index i) const;
  inline bool hasSameNumVars(const MapPowers& that) const;

  //basic comparators
  friend bool operator!=(const MapPowers& a, const MapPowers& b);
  friend bool isLessThanLexicographically(const MapPowers& a, const MapPowers& b);

  //arithmetic assignments
  MapPowers& operator*=(MapPowers& that); //that could be self
  MapPowers& operator/=(MapPowers& that); //that could be self
  MapPowers& operator*=(const MapPowers& that); //that not self
  MapPowers& operator/=(const MapPowers& that); //that not self

  //utility members
  Power degree(void) const;
  std::pair<Coeff, MapPowers> diff(const Index i) const;
  MapPowers pow(const Power p) const;

  //evaluation and substitution
  Coeff operator()(const std::vector< Coeff >& v) const;

  //factory
  inline static MapPowers CoordinatePower(const Index numVars,
                                          const Index i,
                                          const Power p = 1);

  //forward declaration of polynomial for use in evaluation at polynomial
  friend class Polynomial;
  Polynomial operator()(const std::vector< Polynomial >& v) const;
  friend std::ostream& operator<<(std::ostream& ostr, const MapPowers& powers);

 protected:
  Index numVars_;
  IndexToPowerMap powers_;

 private:
}; //MapPowers

//all binary operators (except friends) must be outside the class [Style]:-
inline bool operator==(const MapPowers& a, const MapPowers& b);
inline bool operator<(const MapPowers& a, const MapPowers& b);

// Derived arithmetic operators
inline MapPowers operator*(const MapPowers& a, const MapPowers& b);
inline MapPowers operator/(const MapPowers& a, const MapPowers& b);

//---inline implementations follow

inline bool MapPowers::isValidIndex(const Index i) const {
  return (i >= 0 && i < numVars_);
}

inline Index MapPowers::getNumVars(void) const {
  return numVars_;
}

inline Index MapPowers::getNumPowersStored(void) const {
  return powers_.size();
}

inline bool MapPowers::hasSameNumVars(const MapPowers& that) const {
  return (numVars_ == that.numVars_);
}

inline Power MapPowers::operator[](const Index i) const {
  return getPower(i);
}

inline bool operator==(const MapPowers& a, const MapPowers& b) {
  return !(a != b);
}

inline bool operator<(const MapPowers& p, const MapPowers& q) {
  //we reverse order here
  return isLessThanLexicographically(q, p);
}

inline bool operator>(const MapPowers& a, const MapPowers& b) {
  //symmetry
  return operator<(b, a);
}

inline bool operator<=(const MapPowers& a, const MapPowers& b) {
  //since < is a total ordering
  return !operator<(b, a);
}

inline bool operator>=(const MapPowers& a, const MapPowers& b) {
  //since < is a total ordering
  return !operator<(a, b);
}

inline MapPowers operator*(const MapPowers& a, const MapPowers& b) {
  MapPowers result(a);
  return (result*=b);
}

inline MapPowers operator/(const MapPowers& a, const MapPowers& b) {
  MapPowers result(a);
  return (result/=b);
}

inline
MapPowers MapPowers::CoordinatePower(const Index numVars,
                                     const Index i,
                                     const Power p) {
  MapPowers xI(numVars);
  xI.setPower(i, p);
  return xI;
}

#endif //MAP_POWERS_H

