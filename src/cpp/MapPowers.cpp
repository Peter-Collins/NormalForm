//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: MapPowers implementation
//
//----------------------------------------------------------------------

//system
#include <iostream>
#include <vector>
#include <assert.h>

//project
#include "Types.h"
#include "Utility.h"
#include "MapPowers.h"
#include "Polynomial.h"

// [4] Note that, for the default-value logic to work correctly, we
// must use find(key) on the map, and must _not_ rely on a lookup
// map[key], as the latter would erroneously inset a new zero-valued
// key into the map, defeating our desire for a sparse representation.

//construct default
MapPowers::MapPowers(void)
  : numVars_(1) {
  //nothing; thus the default is $x_{0}^{0}$, i.e., 1-variable version of 1.
}

//construct copy
MapPowers::MapPowers(const MapPowers& that)
  : numVars_(that.numVars_) {
  IndexToPowerMap::const_iterator iter;
  for (iter = that.powers_.begin(); iter != that.powers_.end(); ++iter) {
    assert(isPositivePower(iter->second));
    powers_[iter->first] = iter->second;
  }
}

//assign
MapPowers& MapPowers::operator=(const MapPowers& that) {
  if (&that != this) {
    numVars_ = that.numVars_;
    powers_.clear();
    IndexToPowerMap::const_iterator iter;
    for (iter = that.powers_.begin(); iter != that.powers_.end(); ++iter) {
      assert(isPositivePower(iter->second));
      powers_[iter->first] = iter->second;
    }
  }
  return *this;
}

//destruct
MapPowers::~MapPowers(void) {
  //nothing
}

//----------------------------------------------------------------
// Other constructors and factories
//----------------------------------------------------------------

MapPowers::MapPowers(const Index numVars)
  : numVars_(numVars) {
  if (!(numVars_ > 0))
    throw MapPowersIndexError();
}

MapPowers::MapPowers(const std::vector< Power >& powers)
  : numVars_(powers.size()) {
  if (!(numVars_ > 0))
    throw MapPowersIndexError();
  for (Index i = 0; i < numVars_; ++i) {
    if (!isNonnegativePower(powers[i]))
      throw MapPowersValueError();
    if (powers[i] > 0)
      powers_[i] = powers[i];
  }
}

MapPowers::MapPowers(const Index numVars, const IndexToPowerMap& powers)
  : numVars_(numVars) {
  IndexToPowerMap::const_iterator iter;
  for (iter = powers.begin(); iter != powers.end(); ++iter) {
    if (!isValidIndex(iter->first))
      throw MapPowersIndexError();
    if (!isNonnegativePower(iter->second))
      throw MapPowersValueError();
    if (iter->second > 0)
      powers_[iter->first] = iter->second;
  }
}

//----------------------------------------------------------------
// Inspectors
//----------------------------------------------------------------

Power MapPowers::getPower(const Index i) const {
  if (!isValidIndex(i))
    throw MapPowersIndexError();
  IndexToPowerMap::const_iterator iter = powers_.find(i);
  if (iter == powers_.end())
    return 0;
  else {
    assert(isPositivePower(iter->second));
    return iter->second;
  }
}

MapPowers& MapPowers::setPower(const Index i, const Power p) {
  if (!isValidIndex(i))
    throw MapPowersIndexError();
  if (!isNonnegativePower(p))
    throw MapPowersValueError();
  IndexToPowerMap::iterator iter = powers_.find(i);
  if (iter == powers_.end()) {
    if (p > 0)
      powers_[i] = p;
  }
  else {
    if (p > 0) {
      assert(isPositivePower(p));
      iter->second = p;
    }
    else
      powers_.erase(i);
  }
  return *this;
}

//----------------------------------------------------------------
// Comparators
//----------------------------------------------------------------

bool operator!=(const MapPowers& p, const MapPowers& q) {
  if (!p.hasSameNumVars(q))
    throw MapPowersSizeMismatchError();
  if (p.powers_ != q.powers_)
    return true;
  return false;
}

bool isLessThanLexicographically(const MapPowers& p, const MapPowers& q) {
  if (!p.hasSameNumVars(q))
    throw MapPowersSizeMismatchError();
  IndexToPowerMap::const_iterator a;
  IndexToPowerMap::const_iterator b;
  const IndexToPowerMap::const_iterator a_end = p.powers_.end();
  const IndexToPowerMap::const_iterator b_end = q.powers_.end();
  for (a = p.powers_.begin(), b = q.powers_.begin();
       (a != a_end) && (b != b_end);
       ++a, ++b) {
    assert(isPositivePower(a->second));
    assert(isPositivePower(b->second));
    //whichever has the later-indexed first non-zero power is the lesser:
    if (a->first > b->first)
      return true;
    if (b->first > a->first)
      return false;
    //the first non-zero indices are equal; compare the powers:
    if (a->second < b->second)
      return true;
    if (b->second < a->second)
      return false;
  }
  //at this point, a or b or both have ended
  if (b != b_end) {
    //b has more powers, so a < b.
    return true;
  }
  if (a != a_end) {
    //a has more powers, so b < a.
    return false;
  }
  //both have run out, with no result.
  return false;
}

//----------------------------------------------------------------
// Arithmetic assignments
//----------------------------------------------------------------

MapPowers& MapPowers::operator*=(MapPowers& that) {
  if (&that == this) {
    IndexToPowerMap::iterator iter;
    for (iter = powers_.begin(); iter != powers_.end(); ++iter) {
      iter->second *= 2;
      assert(isPositivePower(iter->second));
    }
  }
  else {
    if (!hasSameNumVars(that))
      throw MapPowersSizeMismatchError();
    IndexToPowerMap::const_iterator q;
    IndexToPowerMap::iterator p;
    for (q = that.powers_.begin(); q != that.powers_.end(); ++q) {
      p = powers_.find(q->first);
      if (p == powers_.end()) {
	assert(isPositivePower(q->second));
	powers_[q->first] = q->second;
      }
      else {
	p->second += q->second;
	assert(isPositivePower(p->second));
      }
    }
  }
  return *this;
}

MapPowers& MapPowers::operator*=(const MapPowers& that) {
  assert(!(&that == this));
  if (!hasSameNumVars(that))
    throw MapPowersSizeMismatchError();
  IndexToPowerMap::const_iterator q;
  IndexToPowerMap::iterator p;
  for (q = that.powers_.begin(); q != that.powers_.end(); ++q) {
    p = powers_.find(q->first);
    if (p == powers_.end()) {
      assert(isPositivePower(q->second));
      powers_[q->first] = q->second;
    }
    else {
      p->second += q->second;
      assert(isPositivePower(p->second));
    }
  }
  return *this;
}

MapPowers& MapPowers::operator/=(MapPowers& that) {
  if (&that == this)
    powers_.clear();
  else {
    if (!hasSameNumVars(that))
      throw MapPowersSizeMismatchError();
    IndexToPowerMap::const_iterator q;
    IndexToPowerMap::iterator p;
    for (q = that.powers_.begin(); q != that.powers_.end(); ++q) {
      p = powers_.find(q->first);
      //std::cerr << "Attempt to divide by absent power." << std::endl;
      assert(p != powers_.end());
      //std::cerr << "Attempt to divide by greater power." << std::endl;
      assert(p->second >= q->second);
      p->second -= q->second;
      if (p->second == 0)
	powers_.erase(p->first);
    }
  }
  return *this;
}

MapPowers& MapPowers::operator/=(const MapPowers& that) {
  assert(!(&that == this));
  if (!hasSameNumVars(that))
    throw MapPowersSizeMismatchError();
  IndexToPowerMap::const_iterator q;
  IndexToPowerMap::iterator p;
  for (q = that.powers_.begin(); q != that.powers_.end(); ++q) {
    p = powers_.find(q->first);
    //std::cerr << "Attempt to divide by absent power." << std::endl;
    assert(p != powers_.end());
    //std::cerr << "Attempt to divide by greater power." << std::endl;
    assert(p->second >= q->second);
    p->second -= q->second;
    if (p->second == 0)
      powers_.erase(p->first);
  }
  return *this;
}

//----------------------------------------------------------------
// Utility members
//----------------------------------------------------------------

Power MapPowers::degree(void) const {
  Power sum = Power(0);
  IndexToPowerMap::const_iterator iter;
  for (iter = powers_.begin(); iter != powers_.end(); ++iter) {
    assert(isPositivePower(iter->second));
    sum += iter->second;
  }
  return sum;
}

std::pair<Coeff, MapPowers> MapPowers::diff(const Index i) const {
  //differentiate with respect to the indexed coordinate
  if (!isValidIndex(i))
    throw MapPowersIndexError();
  IndexToPowerMap::const_iterator iter = powers_.find(i);
  if (iter == powers_.end()) {
    MapPowers res(numVars_);
    return std::pair<Coeff, MapPowers>(CoeffZero, res);
  }
  else {
    assert(isPositivePower(iter->second));
    MapPowers res(*this);
    //set copes with zero
    res.setPower(i, (iter->second)-1);
    return std::pair<Coeff, MapPowers>(Coeff(iter->second), res);
  }
}

MapPowers MapPowers::pow(const Power p) const {
  if (!isNonnegativePower(p))
    throw MapPowersValueError();
  MapPowers res(numVars_);
  if (p > 0) {
    IndexToPowerMap::const_iterator iter;
    for (iter = powers_.begin(); iter != powers_.end(); ++iter) {
      assert(isPositivePower(iter->second));
      res.powers_[iter->first] = p*(iter->second);
    }
  }
  return res;
}

Coeff MapPowers::operator()(const std::vector< Coeff >& v) const {
  if (!(v.size() == numVars_))
    throw MapPowersSizeMismatchError();
  Coeff result(CoeffOne);
  IndexToPowerMap::const_iterator ip;
  for (ip = powers_.begin(); ip != powers_.end(); ++ip) {
    assert(isValidIndex(ip->first));
    assert(isPositivePower(ip->second));
    result *= coeffPow(v[ip->first], ip->second);
  }
  return result;
}

Polynomial MapPowers::operator()(const std::vector< Polynomial >& v) const {
  if (!(v.size() == numVars_))
    throw MapPowersSizeMismatchError();
  //important: the number of variables depends on the elements of v!
  Polynomial result = Polynomial::One(v[0].getNumVars());
  IndexToPowerMap::const_iterator ip;
  for (ip = powers_.begin(); ip != powers_.end(); ++ip) {
    assert(isValidIndex(ip->first));
    assert(isPositivePower(ip->second));
    result *= v[ip->first].pow(ip->second);
  }
  return result;
}

std::ostream& operator<<(std::ostream& ostr, const MapPowers& powers) {
  bool isFirst = true;
  IndexToPowerMap::const_iterator ip;
  for (ip = powers.powers_.begin(); ip != powers.powers_.end(); ++ip) {
    if (isFirst)
      isFirst = false;
    else
      ostr << " ";
    if ((ip->second) > 1)
      ostr << "x_{" << (ip->first) << "}^{" << (ip->second) << "}";
    else
      ostr << "x_{" << (ip->first) << "}";
  }
  return ostr;
}

