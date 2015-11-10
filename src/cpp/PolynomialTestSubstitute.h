//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef POLYNOMIAL_TEST_SUBSTITUTE_H
#define POLYNOMIAL_TEST_SUBSTITUTE_H

//system
#include <utility> //for pair
#include <set>

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"

static const size_t poly_substitute_n_cases = 8;
static const Index poly_substitute_max_vars = 10;
static const Real poly_substitute_tol = 1.0e-14;

class PolynomialTestSubstitute : public CppUnit::TestFixture {

 public:

  //The following create static CppUnit::TestSuite *suite():
  CPPUNIT_TEST_SUITE(PolynomialTestSubstitute);

  //Substitution / evaluation at a vector of Polynomials:
  CPPUNIT_TEST(test_substitute_identity);
  CPPUNIT_TEST(test_substitute_squares);
  CPPUNIT_TEST(test_substitute_example);
  CPPUNIT_TEST(test_substitute_example2);
  CPPUNIT_TEST(test_substitute_sum);
  CPPUNIT_TEST(test_substitute_neg);

  //End the test suite:
  CPPUNIT_TEST_SUITE_END();

  void test_substitute_identity(void) {
    for (size_t n=0; n<poly_substitute_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_substitute_max_vars);
      const Polynomial a = randomPoly(vars);
      std::vector< Polynomial > id;
      for (Index i=0; i<vars; ++i) {
	id.push_back(Polynomial::CoordinateMonomial(vars, i));
	CPPUNIT_ASSERT(id[i].getNumVars() == vars);
      }
      CPPUNIT_ASSERT(id.size() == vars);
      Polynomial b = a(id);
      CPPUNIT_ASSERT(b.getNumVars() == a.getNumVars());
      CPPUNIT_ASSERT(b.getNumTerms() == a.getNumTerms());
      CPPUNIT_ASSERT((b - a).lInfinityNorm() < poly_substitute_tol);
    }
  }

  void test_substitute_squares(void) {
    for (size_t n=0; n<poly_substitute_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_substitute_max_vars);
      const Polynomial a = randomPoly(vars);
      std::vector< Polynomial > sq;
      for (Index i=0; i<vars; ++i) {
	sq.push_back(Polynomial::CoordinateMonomial(vars, i, 2));
      }
      CPPUNIT_ASSERT(sq.size() == vars);
      Polynomial b = a(sq);
      CPPUNIT_ASSERT(b.getNumTerms() == a.getNumTerms());
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(fabs(b[(pc->first).pow(2)] - pc->second) < poly_substitute_tol);
      }
    }
  }

  void test_substitute_example(void) {
    const Polynomial x = Polynomial::CoordinateMonomial(2, 0);
    const Polynomial y = Polynomial::CoordinateMonomial(2, 1);
    Polynomial p = x+y;
    Polynomial s = Polynomial::CoordinateMonomial(1, 0, 2);
    Polynomial q = x*x + Coeff(2.0)*x*y + y*y;
    Polynomial c(p);
    c *= p;
    CPPUNIT_ASSERT(q == c);
    CPPUNIT_ASSERT(q == p*p);
    CPPUNIT_ASSERT(q == p.pow(2));
    CPPUNIT_ASSERT(s(std::vector< Polynomial >(1, p)) == q);
  }

  void test_substitute_example2(void) {
    const Polynomial x = Polynomial::CoordinateMonomial(2, 0);
    const Polynomial y = Polynomial::CoordinateMonomial(2, 1);
    Polynomial p = x+y;
    Polynomial s = Polynomial::CoordinateMonomial(1, 0, 3);
    Polynomial q = x*x*x + Coeff(3.0)*x*x*y + Coeff(3.0)*x*y*y + y*y*y;
    Polynomial c(p);
    c *= p;
    c *= p;
    CPPUNIT_ASSERT(q == c);
    CPPUNIT_ASSERT(q == p*p*p);
    CPPUNIT_ASSERT(q == p.pow(3));
    CPPUNIT_ASSERT(s(std::vector< Polynomial >(1, p)) == q);
  }

  void test_substitute_sum(void) {
    const Index vars = 3;
    Polynomial p(vars);
    const Coeff a(-0.3);
    const Coeff b(0.995);
    std::vector< Power > powers;
    init(powers) = 2, 1, 3;
    p += a*Polynomial(powers);
    init(powers) = 1, 1, 0;
    p += b*Polynomial(powers);
    std::vector< Polynomial > ids;
    for (Index i=0; i<vars; ++i) {
      init(ids) += Polynomial::CoordinateMonomial(vars, i);
    }
    CPPUNIT_ASSERT(ids.size() == vars);
    ids[0] += Coeff(2.0)*Polynomial::CoordinateMonomial(vars, 2);
    Polynomial q = p(ids);
    CPPUNIT_ASSERT(q.getNumTerms() == 5);
    init(powers) = 2, 1, 3;
    CPPUNIT_ASSERT(q[powers] == a);
    init(powers) = 1, 1, 4;
    CPPUNIT_ASSERT(q[powers] == Coeff(4.0)*a);
    init(powers) = 0, 1, 5;
    CPPUNIT_ASSERT(q[powers] == Coeff(4.0)*a);
    init(powers) = 1, 1, 0;
    CPPUNIT_ASSERT(q[powers] == b);
    init(powers) = 0, 1, 1;
    CPPUNIT_ASSERT(q[powers] == Coeff(2.0)*b);
  }

  void test_substitute_neg(void) {
    const Index vars = 3;
    Polynomial p(vars);
    const Coeff a(-0.3);
    const Coeff b(0.995);
    std::vector< Power > powers;
    init(powers) = 2, 2, 3;
    p += a*Polynomial(powers);
    init(powers) = 1, 1, 0;
    p += b*Polynomial(powers);
    std::vector< Polynomial > neg_ids;
    for (Index i=0; i<vars; ++i) {
      init(neg_ids) += -Polynomial::CoordinateMonomial(vars, i);
    }
    Polynomial q = p(neg_ids);
    CPPUNIT_ASSERT(q.getNumTerms() == p.getNumTerms());
    init(powers) = 2, 2, 3;
    CPPUNIT_ASSERT(q[powers] == -a);
    init(powers) = 1, 1, 0;
    CPPUNIT_ASSERT(q[powers] == b);
  }
        
}; //PolynomialTestSubstitute

#endif //POLYNOMIAL_TEST_SUBSTITUTE_H
