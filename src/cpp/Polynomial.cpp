//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Polynomial implementation
//
//----------------------------------------------------------------------

//system
#include <assert.h>
#include <iostream>

//project
#include "Polynomial.h"

//----------------------------------------------------------------------
//the big four
//----------------------------------------------------------------------

Polynomial::Polynomial(void)
  : numVars_(1) {
}

Polynomial::Polynomial(const Polynomial& that)
  : numVars_(that.numVars_) {
  PowersToCoeffMap::const_iterator iter;
  for (iter = that.terms_.begin(); iter != that.terms_.end(); ++iter) {
    assert(hasSameNumVars(iter->first));
    assert(isValidCoeff(iter->second));
    terms_.insert(PowersToCoeffMap::value_type(iter->first, iter->second));
  }
  assert(getNumTerms() == that.getNumTerms());
}

Polynomial& Polynomial::operator=(const Polynomial& that) {
  if (&that != this) {
    if (!hasSameNumVars(that))
      throw PolynomialSizeMismatchError();
    terms_.clear();
    assert(terms_.empty());
    PowersToCoeffMap::const_iterator iter;
    for (iter = that.terms_.begin(); iter != that.terms_.end(); ++iter) {
      assert(hasSameNumVars(iter->first));
      assert(isValidCoeff(iter->second));
      terms_.insert(PowersToCoeffMap::value_type(iter->first, iter->second));
    }
  }
  assert(getNumTerms() == that.getNumTerms());
  return *this;
}

Polynomial::~Polynomial(void) {
  //do nothing
}

//----------------------------------------------------------------------
//constructors and factories:
//----------------------------------------------------------------------

//zero polynomial
Polynomial::Polynomial(const Index numVars)
  : numVars_(numVars) {
  if (!(numVars_ > 0))
    throw PolynomialIndexError();
}

//unit monomial
Polynomial::Polynomial(const Powers& powers)
  : numVars_(powers.numVars_) {
  terms_.insert(PowersToCoeffMap::value_type(powers, CoeffOne));
}

//monomial
Polynomial::Polynomial(const Powers& powers, const Coeff& coeff)
  : numVars_(powers.numVars_) {
  if (!hasSameNumVars(powers))
    throw PolynomialSizeMismatchError();
  if (isValidCoeff(coeff)) {
    terms_.insert(PowersToCoeffMap::value_type(powers, coeff));
  }
}

Polynomial& Polynomial::setMonomial(const Powers& powers, const Coeff& coeff) {
  if (!hasSameNumVars(powers))
    throw PolynomialSizeMismatchError();
  PowersToCoeffMap::iterator pc = terms_.find(powers);
  if (pc == terms_.end()) {
    //not found
    if (isValidCoeff(coeff))
      terms_.insert(PowersToCoeffMap::value_type(powers, coeff));
  }
  else {
    //found
    if (isValidCoeff(coeff))
      pc->second = coeff;
    else
      terms_.erase(pc->first);
  }
  return *this;
}

Polynomial& Polynomial::addMonomial(const Powers& powers, const Coeff& coeff) {
  if (!hasSameNumVars(powers))
    throw PolynomialSizeMismatchError();
  PowersToCoeffMap::iterator pc = terms_.find(powers);
  if (pc == terms_.end()) {
    //not found
    if (isValidCoeff(coeff))
      terms_.insert(PowersToCoeffMap::value_type(powers, coeff));
  }
  else {
    //found
    pc->second += coeff;
    if (!isValidCoeff(pc->second))
      terms_.erase(pc->first);
  }
  return *this;
}

Polynomial Polynomial::Zero(const Index numVars) {
  return Polynomial(numVars);
}

Polynomial Polynomial::One(const Index numVars) {
  MapPowers zeroPowers(numVars);
  return Polynomial(zeroPowers);
}

Polynomial Polynomial::CoordinateMonomial(const Index numVars,
					  const Index i,
					  const Power p) {
  if (!(numVars > 0))
    throw PolynomialIndexError();
  if (!isNonnegativePower(p))
    throw PolynomialValueError();
  MapPowers powers(numVars);
  powers.setPower(i, p);
  return Polynomial(powers);
}

//----------------------------------------------------------------------
//inspectors
//----------------------------------------------------------------------

const Coeff& Polynomial::getCoeff(const Powers& powers) const {
  if (!hasSameNumVars(powers))
    throw PolynomialSizeMismatchError();
  PowersToCoeffMap::const_iterator iter = terms_.find(powers);
  if (iter == terms_.end())
    return CoeffZero;
  assert(isValidCoeff(iter->second));
  return iter->second;
}

bool Polynomial::hasPowers(const Powers& powers) const {
  if (!hasSameNumVars(powers))
    throw PolynomialSizeMismatchError();
  PowersToCoeffMap::const_iterator iter = terms_.find(powers);
  return (iter != terms_.end());
}

//----------------------------------------------------------------------
//call; evaluation at a vector of scalars or polynomials
//----------------------------------------------------------------------

Coeff Polynomial::operator()(const std::vector< Coeff >& v) const {
  if (!(v.size() == numVars_))
    throw PolynomialSizeMismatchError();
  Coeff result(CoeffZero);
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    assert(hasSameNumVars(pc->first));
    result += (pc->second) * ((pc->first)(v));
  }
  return result;
}

Polynomial Polynomial::operator()(const std::vector< Polynomial >& v) const {
  if (!(v.size() == numVars_))
    throw PolynomialSizeMismatchError();
  //important: the number of variables depends on the elements of v!
  Polynomial result(v[0].getNumVars());
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    result += (pc->second) * ((pc->first)(v));
  }
  return result;
}

//----------------------------------------------------------------------
//unary operators
//----------------------------------------------------------------------

Polynomial Polynomial::operator-(void) const {
  Polynomial res(numVars_);
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    res.terms_[pc->first] = -(pc->second);
  }
  assert(hasSameNumVars(res));
  assert(res.terms_.size() == terms_.size());
  return res;
}

//----------------------------------------------------------------------
//arithmetic assigns
//----------------------------------------------------------------------

Polynomial& Polynomial::operator+=(const Polynomial& that) {
  if (&that == this) {
    PowersToCoeffMap::iterator pc;
    for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
      assert(isValidCoeff(pc->second));
      (pc->second) *= 2.0;
    }
  }
  else {
    if (!hasSameNumVars(that))
      throw PolynomialSizeMismatchError();
    PowersToCoeffMap::const_iterator thatPc;
    for (thatPc = that.terms_.begin(); thatPc != that.terms_.end(); ++thatPc) {
      assert(isValidCoeff(thatPc->second));
      addMonomial(thatPc->first, thatPc->second);
    }
  }
  return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& that) {
  //does not matter if &that == this
  this->operator=((*this) * that);
  return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& that) {
  if (&that == this) {
    terms_.clear();
    assert(terms_.empty());
  }
  else {
    if (!hasSameNumVars(that))
      throw PolynomialSizeMismatchError();
    PowersToCoeffMap::const_iterator thatPc;
    for (thatPc = that.terms_.begin(); thatPc != that.terms_.end(); ++thatPc) {
      assert(isValidCoeff(thatPc->second));
      addMonomial(thatPc->first, -thatPc->second);
    }
  }
  return *this;
}

Polynomial& Polynomial::operator*=(const Coeff& c) {
  if (!isValidCoeff(c)) {
    terms_.clear();
    assert(terms_.empty());
  }
  else {
    PowersToCoeffMap::iterator pc;
    PowersToCoeffMap::iterator pcTemp;
    for (pc = terms_.begin(); pc != terms_.end();) {
      assert(isValidCoeff(pc->second));
      //GOTCHA: careful not to introduce a tiny coeff!
      const Coeff newCoeff((pc->second) * c);
      if (isValidCoeff(newCoeff)) {
        (pc->second) = newCoeff;
        ++pc;
      }
      else {
        pcTemp = pc;
        ++pcTemp;
        //delete current element; invalidates the old iterator
        terms_.erase(pc->first);
        pc = pcTemp;
      }
    }
  }
  return *this;
}

//----------------------------------------------------------------------
//graded (degree/homogeneous) structure
//----------------------------------------------------------------------

bool Polynomial::isConstant(void) const {
  if (terms_.empty())
    return true;
  if (terms_.size() > 1)
    return false;
  PowersToCoeffMap::const_iterator pc = terms_.begin();
  assert(hasSameNumVars(pc->first));
  assert(isValidCoeff(pc->second));
  for (Index i = 0; i < numVars_; ++i)
    if ((pc->first)[i] > 0)
      return false;
  return true;
}

bool Polynomial::isHomogeneous(void) const {
  if (terms_.size() <= 1)
    return true;
  PowersToCoeffMap::const_iterator pc = terms_.begin();
  assert(hasSameNumVars(pc->first));
  assert(isValidCoeff(pc->second));
  Power deg_found = (pc->first).degree();
  Power deg;
  for (++pc; pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    deg = (pc->first).degree();
    if (deg != deg_found)
      return false;
  }
  return true;
}

bool Polynomial::isHomogeneous(const Power p) const {
  if (!isNonnegativePower(p))
    throw PolynomialValueError();
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    if ((pc->first).degree() != p)
      return false;
  }
  return true;
}

Power Polynomial::degree(void) const {
  Power maxDegree(0);
  PowersToCoeffMap::const_iterator pc;
  Power degree;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    degree = (pc->first).degree();
    if (degree > maxDegree) {
      maxDegree = degree;
    }
  }
  return maxDegree;
}

Polynomial Polynomial::homogeneousPart(const Power deg,
				       const Power lessThan) const {
  if (!isNonnegativePower(deg))
    throw PolynomialValueError();
  if (!isNonnegativePower(lessThan))
    throw PolynomialValueError();
  Power upper = (deg < lessThan) ? lessThan : deg+1;
  Polynomial res(numVars_);
  Power degree;
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    degree = (pc->first).degree();
    if ((deg <= degree) && (degree < upper))
      res.terms_.insert(PowersToCoeffMap::value_type(pc->first, pc->second));
  }
  return res;
}

//----------------------------------------------------------------------
//binary arithmetic operators
//----------------------------------------------------------------------

Polynomial operator*(const Polynomial& a, const Polynomial& b) {
  if (!a.hasSameNumVars(b))
    throw PolynomialSizeMismatchError();
  Polynomial res(a.numVars_);
  PowersToCoeffMap::const_iterator apc;
  PowersToCoeffMap::const_iterator bpc;
  for (apc = a.terms_.begin(); apc != a.terms_.end(); ++apc) {
    assert(isValidCoeff(apc->second));
    for (bpc = b.terms_.begin(); bpc != b.terms_.end(); ++bpc) {
      assert(isValidCoeff(bpc->second));
      res.addMonomial((apc->first) * (bpc->first),
                      (apc->second) * (bpc->second));
    }
  }
  return res;
}

//----------------------------------------------------------------------
//norms
//----------------------------------------------------------------------

Real Polynomial::lOneNorm(void) const {
  Real sum(RealZero);
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    sum += fabs(pc->second);
  }
  return sum;
}

Real Polynomial::lInfinityNorm(void) const {
  Real largest(RealZero);
  Real current;
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    current = fabs(pc->second); //bad!
    if (current > largest) {
      largest = current;
    }
  }
  return largest;
}

//----------------------------------------------------------------------
//i/o
//----------------------------------------------------------------------

std::ostream& operator<<(std::ostream& ostr, const Polynomial& poly) {
  bool isFirst = true;
  PowersToCoeffMap::const_iterator pc = poly.terms_.begin();
  ostr << "= ";
  if (pc == poly.terms_.end())
    return (ostr << "0\n");
  for (; pc != poly.terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    if (isFirst)
      isFirst = false;
    else
      ostr << "+ ";
    ostr << (pc->second) << " " << (pc->first) << "\n";
  }
  return ostr;
}

//----------------------------------------------------------------------
//miscellany
//----------------------------------------------------------------------

Polynomial Polynomial::diff(const Index i) const {
  if (!isValidIndex(i))
    throw PolynomialIndexError();
  Polynomial res(numVars_);
  std::pair< Coeff, MapPowers > diffCp;
  PowersToCoeffMap::const_iterator pc;
  PowersToCoeffMap::iterator currentPc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(hasSameNumVars(pc->first));
    assert(isValidCoeff(pc->second));
    //coeff and diff powers
    diffCp = (pc->first).diff(i);
    if (diffCp.first != CoeffZero) {
      res.addMonomial(diffCp.second, (pc->second) * (diffCp.first));
    }
  }
  return res;
}

std::pair< Polynomial, Polynomial >
Polynomial::partition(const PowersAndCoeffToBool& pred) const {
  Polynomial polyTrue(numVars_);
  Polynomial polyFalse(numVars_);
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms_.begin(); pc != terms_.end(); ++pc) {
    assert(isValidCoeff(pc->second));
    if ((*pred)(pc->first, pc->second))
      polyTrue.setMonomial(pc->first, pc->second);
    else
      polyFalse.setMonomial(pc->first, pc->second);
  }
  assert(polyTrue.getNumTerms() + polyFalse.getNumTerms() == getNumTerms());
  return std::pair< Polynomial, Polynomial >(polyTrue, polyFalse);
}

Polynomial Polynomial::pow(const Power p) const {
  if (!isNonnegativePower(p))
    throw PolynomialValueError();
  if (p == 1)
    return Polynomial(*this);
  Polynomial res = Polynomial::One(numVars_);
  if (p > 0) {
    Power mask(1);
    //position of p's most-significant bit
    while (mask <= p)
      mask <<= 1;
    mask >>= 1;
    //use repeated squaring for binary expansion power
    while (mask > 0) {
      res *= res;
      if (p & mask)
	res *= (*this);
      mask >>= 1;
    }
  }
  return res;
}

