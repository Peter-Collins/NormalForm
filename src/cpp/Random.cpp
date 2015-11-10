//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Random implementation
//
//----------------------------------------------------------------------

//standard headers
#include <vector>
#include <iostream>

//library headers

//project headers
#include "Types.h"
#include "Random.h"
#include "Polynomial.h"

std::vector< Power > randomVectorPowers(const Index length,
					const Power maxDegree) {
  //
  // Note that this does not choose monomials with equal probability,
  // but does ensure that the total degree of the output monomial is
  // in the specified bounds.
  //
  assert(length > 0);
  assert(maxDegree >= 0);
  std::vector< Power > v(length);
  for (Index i = 0; i < length; ++i)
    v[i] = 0;
  const Power deg = randomPower(maxDegree);
  Index j;
  for (Power p = 0; p < deg; ++p) {
    j = randomIndex(length);
    v[j] = (v[j] + 1);
  }
  assert(v.size() == length);
  return v;
}

MapPowers randomPowers(const Index length, const Power maxDegree) {
  assert(length > 0);
  assert(maxDegree >= 0);
  MapPowers powers(randomVectorPowers(length, maxDegree));
  assert(powers.getNumVars() == length);
  assert(powers.degree() < maxDegree);
  return powers;
}

Polynomial randomPoly(const Index numVars,
		      const Index maxNumVars,
		      const Power maxDegree,
		      const size_t maxNumTerms,
		      const Real maxAbsCoeff) {
  assert(numVars >= 0);
  assert(maxNumVars >= 0);
  assert(maxDegree >= 0);
  assert(maxNumTerms >= 0);
  Index numVars_;
  if (numVars == 0) {
    numVars_ = randomIndex(maxNumVars);
  } else {
    numVars_ = numVars;
  }
  Polynomial poly(numVars_);
  size_t numTerms_ = randomSize(maxNumTerms);
  for (size_t term = 0; term < numTerms_; ++term) {
    MapPowers powers = randomPowers(numVars_, maxDegree);
    poly.addMonomial(powers, randomCoeff(maxAbsCoeff));
  }
  assert(poly.getNumVars() == numVars_);
  assert(poly.getNumTerms() <= numTerms_);
  return poly;
}
