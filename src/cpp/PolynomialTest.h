//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialTest
//
//----------------------------------------------------------------------

#ifndef POLYNOMIAL_TEST_H
#define POLYNOMIAL_TEST_H

//system headers
#include <utility>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"

class PolynomialTest : public CppUnit::TestFixture {
  //The following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(PolynomialTest);
  //Construction:
  CPPUNIT_TEST(testExampleMonomial);
  CPPUNIT_TEST(testExampleMonomialUnit);
  CPPUNIT_TEST(testMakeConstantZero);
  CPPUNIT_TEST(testMakeConstantOne);
  CPPUNIT_TEST_EXCEPTION(testMakeZeroLengthPolynomial,
			 PolynomialIndexError);
  CPPUNIT_TEST(testMonomialCreationUsingStlInit);
  //Basic inspectors:
  CPPUNIT_TEST(testDefaultNumVars);
  CPPUNIT_TEST(testDefaultNumTerms);
  CPPUNIT_TEST(testDefaultPowers);
  CPPUNIT_TEST(testExampleNumTerms);
  CPPUNIT_TEST(testExampleNumVars);
  CPPUNIT_TEST(testPolyWithNoTermsIsConstant);
  CPPUNIT_TEST(testPolyWithNoTermsIsZeroDegree);
  CPPUNIT_TEST(testPolyWithNoTermsIsFalse);
  CPPUNIT_TEST(testPolyWithSomeTermsIsTrue);
  CPPUNIT_TEST(testIsConstantZero);
  CPPUNIT_TEST(testIsConstantOneViaFactory);
  CPPUNIT_TEST(testIsConstantTooManyTerms);
  CPPUNIT_TEST(testIsConstantSingleNonConstTerm);
  //Comparators:
  CPPUNIT_TEST(testDefaultEqualsSelf);
  CPPUNIT_TEST(testDefaultEqualsSame);
  CPPUNIT_TEST(testDefaultNotEqOne);
  //Key existence:
  CPPUNIT_TEST(testMonomialGetItemExists);
  CPPUNIT_TEST(testZeroMonomialGetItemNotExists);
  CPPUNIT_TEST(testExampleGetItemExists);
  CPPUNIT_TEST(testExampleGetOtherItemNotExists);
  //Coefficient access:
  CPPUNIT_TEST(testExampleHasCoeffOne);
  //Call:
  CPPUNIT_TEST(testExampleCall);
  CPPUNIT_TEST(testExampleCallCountNonZeros);
  CPPUNIT_TEST(testExampleCallSameSumMonomials);
  CPPUNIT_TEST(testCallNegativeCount);
  CPPUNIT_TEST(testCallMul);
  //Diff:
  CPPUNIT_TEST(testDiffSum);
  //Degree:
  CPPUNIT_TEST(testExampleDegree);
  CPPUNIT_TEST(testDegreeOne);
  //Norms:
  CPPUNIT_TEST(testLInfinityZero);
  CPPUNIT_TEST(testLOneZero);
  CPPUNIT_TEST(testLInfinityExample);
  CPPUNIT_TEST(testLOneExample);
  CPPUNIT_TEST_SUITE_END();
private:
protected:
  void testDefaultNumTerms(void) {
    Polynomial p;
    CPPUNIT_ASSERT(p.getNumTerms() == 0);
  }
  void testDefaultNumVars(void) {
    Polynomial p;
    CPPUNIT_ASSERT(p.getNumVars() == 1);
  }
  void testDefaultPowers(void) {
    MapPowers m;
    Polynomial p(m);
    CPPUNIT_ASSERT(p.getNumVars() == 1);
    CPPUNIT_ASSERT(p.getNumTerms() == 1);
    MapPowers n;
    Polynomial q(n, CoeffOne);
    CPPUNIT_ASSERT(p == q);
  }
  void testDefaultEqualsSelf(void) {
    Polynomial p;
    CPPUNIT_ASSERT(p == p);
  }
  void testDefaultEqualsSame(void) {
    Polynomial p;
    Polynomial q;
    CPPUNIT_ASSERT(p == q);
  }
  void testDefaultNotEqOne(void) {
    Polynomial zero;
    MapPowers zero_power;
    Polynomial one(zero_power);
    CPPUNIT_ASSERT(!(zero == one));
    CPPUNIT_ASSERT(zero != one);
  }
  void testMonomialGetItemExists(void) {
    MapPowers m;
    Polynomial p(m);
    CPPUNIT_ASSERT(p.hasPowers(m));
    CPPUNIT_ASSERT(p[m] == CoeffOne);
  }
  void testZeroMonomialGetItemNotExists(void) {
    MapPowers m;
    Polynomial p(m, CoeffZero);
    CPPUNIT_ASSERT(!(p.hasPowers(m)));
  }
  void testExampleNumTerms(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 3);
    Polynomial p(m);
    CPPUNIT_ASSERT(p.getNumTerms() == 1);
  }
  void testExampleNumVars(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 3);
    Polynomial p(m);
    CPPUNIT_ASSERT(p.getNumVars() == 3);
  }
  void testExampleHasCoeffOne(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 3);
    Polynomial p(m);
    CPPUNIT_ASSERT(p[m] == CoeffOne);
  }
  void testExampleGetItemExists(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 3);
    Polynomial p(m);
    CPPUNIT_ASSERT(p.hasPowers(m));
  }
  void testExampleGetOtherItemNotExists(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 3);
    Polynomial p(m);
    m.setPower(1, 3);
    CPPUNIT_ASSERT(!(p.hasPowers(m)));
    CPPUNIT_ASSERT(p[m] != CoeffOne);
    CPPUNIT_ASSERT(p[m] == CoeffZero);
    CPPUNIT_ASSERT(p.getNumTerms() == 1);
    CPPUNIT_ASSERT(p.getNumVars() == 3);
  }
  void testExampleCall(void) {
    MapPowers pows(3);
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    Polynomial p(3);
    p.setMonomial(pows, Coeff(-1.5));
    std::vector< Coeff > v(3);
    v[0] = CoeffOne;
    v[1] = CoeffOne;
    v[2] = CoeffOne;
    CPPUNIT_ASSERT(v.size() == pows.getNumVars());
    CPPUNIT_ASSERT(p(v) == Coeff(-1.5));
  }
  void testExampleMonomial(void) {
    MapPowers pows(3);
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    Polynomial p(3);
    p.setMonomial(pows, Coeff(-1.5));
    Polynomial q(pows, Coeff(-1.5));
    CPPUNIT_ASSERT(q == p);
  }
  void testExampleMonomialUnit(void) {
    MapPowers pows(3);
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    Polynomial p(3);
    p.setMonomial(pows, CoeffOne);
    Polynomial q(pows);
    CPPUNIT_ASSERT(q == p);
  }
  void testZeroCall(void) {
    Polynomial p(3);
    std::vector< Coeff > v(3);
    v[0] = CoeffOne;
    v[1] = -99.0;
    v[2] = 123.456;
    CPPUNIT_ASSERT(p(v) == CoeffZero);
  }
  void testExampleCallCountNonZeros(void) {
    Polynomial p(3);
    MapPowers pows(3);
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    p.setMonomial(pows, CoeffOne);
    CPPUNIT_ASSERT(p[pows] == CoeffOne);
    pows.setPower(0, 0);
    pows.setPower(1, 7);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(-17.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(-17.0));
    pows.setPower(0, 0);
    pows.setPower(1, 0);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(42.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(42.0));
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 5);
    p.setMonomial(pows, CoeffZero);
    CPPUNIT_ASSERT(p[pows] == CoeffZero);
    CPPUNIT_ASSERT(p.getNumTerms() == 3);
    std::vector< Coeff > v(3);
    for (Index i=0; i<3; ++i) {
      v[i] = CoeffOne;
    }
     
    CPPUNIT_ASSERT(p(v) == CoeffOne+Coeff(-17.0)+Coeff(42.0));
  }
  void testExampleCallSameSumMonomials(void) {
    Polynomial p(3);
    MapPowers pows(3);
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    p.setMonomial(pows, CoeffOne);
    CPPUNIT_ASSERT(p[pows] == CoeffOne);
    Polynomial p0(pows);
    Polynomial p00(pows, CoeffOne);
    pows.setPower(0, 0);
    pows.setPower(1, 7);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(-17.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(-17.0));
    Polynomial p1(pows);
    Polynomial p11(pows, Coeff(-17.0));
    pows.setPower(0, 0);
    pows.setPower(1, 0);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(42.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(42.0));
    Polynomial p2(pows);
    Polynomial p22(pows, Coeff(42.0));
    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 5);
    p.setMonomial(pows, CoeffZero);
    CPPUNIT_ASSERT(p[pows] == CoeffZero);
    Polynomial p3(pows);
    Polynomial p33(pows, CoeffZero);
    CPPUNIT_ASSERT(p.getNumTerms() == 3);
    std::vector< Coeff > v(3);
    for (Index i=0; i<3; ++i) {
      v[i] = i;
    }
    CPPUNIT_ASSERT(p(v) == CoeffOne*p0(v)+Coeff(-17.0)*p1(v)+Coeff(42.0)*p2(v)+CoeffZero*p3(v));
    CPPUNIT_ASSERT(p(v) == p00(v)+p11(v)+p22(v)+p33(v));
  }
  void testMakeConstantZero(void) {
    Polynomial zero = Polynomial::Zero(1);
    CPPUNIT_ASSERT(zero.getNumTerms() == 0);
    std::vector< Coeff > v(1);
    v[0] = Coeff(4.5);
    CPPUNIT_ASSERT(zero(v) == CoeffZero);
  }
  void testMakeConstantOne(void) {
    Polynomial one = Polynomial::One(1);
    CPPUNIT_ASSERT(one.getNumTerms() == 1);
    MapPowers m(1);
    m.setPower(0, 0);
    CPPUNIT_ASSERT(one[m] == CoeffOne);
    std::vector< Coeff > v(1);
    v[0] = Coeff(4.5);
    CPPUNIT_ASSERT(one(v) == CoeffOne);
  }
  void testMakeZeroLengthPolynomial(void) {
    Polynomial p(0); //throws index error
  }
  void testCallNegativeCount(void) {
    for (Index length=1; length<10; ++length) {
      MapPowers ones(length);
      for (Index j=0; j<length; ++j) {
	ones.setPower(j, 1);
      }
      Polynomial p(ones);
      p *= Coeff(2.0);
      std::vector< Coeff > v(length);
      for (Index j=0; j<length; ++j) {
	v[j] = Coeff(-1.0);
      }
      CPPUNIT_ASSERT(p(v) == ((length%2) ? Coeff(-2.0) : Coeff(+2.0)));
    }
  }
  void testPolyWithNoTermsIsConstant(void) {
    Polynomial p(6);
    CPPUNIT_ASSERT(p.isConstant());
  }
  void testPolyWithNoTermsIsZeroDegree(void) {
    Polynomial p(6);
    CPPUNIT_ASSERT(p.degree() == 0);
  }
  void testPolyWithNoTermsIsFalse(void) {
    const Polynomial t(52);
    CPPUNIT_ASSERT(t.isZero());
    Powers powers(53);
    powers.setPower(31, 1);
    Polynomial p(powers);
    CPPUNIT_ASSERT(!(p.isZero()));
    p *= CoeffZero;
    CPPUNIT_ASSERT(p.isZero());
  }
  void testPolyWithSomeTermsIsTrue(void) {
    Powers powers(53);
    powers.setPower(31, 1);
    Polynomial p(powers);
    CPPUNIT_ASSERT(!(p.isZero()));
    p *= CoeffZero;
    CPPUNIT_ASSERT(p.isZero());
  }
  void testMulFloat(void) {
    Powers powers(17);
    powers.setPower(3, 1);
    Polynomial t(powers);
    CPPUNIT_ASSERT(t.getNumTerms() == 1);
    t *= CoeffZero;
    CPPUNIT_ASSERT(t.getNumTerms() == 0);
    CPPUNIT_ASSERT(t.isZero());
  }
  void testMonomialCreationUsingStlInit(void) {
    std::vector< Power > powers;
    init(powers) = 0, 0, 0, 1, 0, 2, 3, 0, 1;
    Polynomial p(powers); //implicit MapPowers creation
    CPPUNIT_ASSERT(p.getNumVars() == 9);
    std::vector< Coeff > v;
    init(v) = CoeffOne, Coeff(2.0), Coeff(3.0), Coeff(4.0), Coeff(5.0), Coeff(6.0), Coeff(7.0), Coeff(8.0), Coeff(9.0);
    CPPUNIT_ASSERT(v.size() == 9);
    Coeff result = p(v);
    CPPUNIT_ASSERT(result == Coeff(4.0*6.0*6.0*7.0*7.0*7.0*9.0));
  }
  void testExampleDegree(void) {
    Polynomial p(5);
    std::vector< Power > powers;
    init(powers) = 1, 0, 1, 1, 0;
    p += Coeff(0.1)*Polynomial(powers);
    init(powers) = 1, 7, 1, 1, 0;
    p -= Coeff(0.2)*Polynomial(powers);
    init(powers) = 0, 0, 1, 1, 2;
    p += Coeff(0.3)*Polynomial(powers);
    init(powers) = 1, 0, 9, 1, 0;
    p -= Coeff(0.4)*Polynomial(powers);
    CPPUNIT_ASSERT(p.degree() == 11);
  }
  void testLInfinityZero(void) {
    Polynomial p(3);
    CPPUNIT_ASSERT(p.lInfinityNorm() == Real(0.0));
  }
  void testLOneZero(void) {
    Polynomial p(3);
    CPPUNIT_ASSERT(p.lOneNorm() == Real(0.0));
  }
  void testLInfinityExample(void) {
    //in this test, we use integer coefficients to ensure exact
    Polynomial p(4);
    std::vector< Power > powers;
    init(powers) = 1, 1, 0, 1;
    p += Coeff(-5)*Polynomial(powers);
    init(powers) = 5, 3, 0, 7;
    p += Coeff(3)*Polynomial(powers);
    CPPUNIT_ASSERT(p.getNumTerms() == 2);
    CPPUNIT_ASSERT(p.lInfinityNorm() == Real(5));
  }
  void testLOneExample(void) {
    //in this test, we use integer coefficients to ensure exact
    Polynomial p(4);
    std::vector< Power > powers;
    init(powers) = 1, 1, 0, 1;
    p += Coeff(-5)*Polynomial(powers);
    init(powers) = 5, 3, 0, 7;
    p += Coeff(3)*Polynomial(powers);
    CPPUNIT_ASSERT(p.getNumTerms() == 2);
    CPPUNIT_ASSERT(p.lOneNorm() == Real(8));
  }
  void testDegreeOne(void) {
    std::vector< Power > powers;
    for (Power j=0; j<5; ++j) {
      for (Index i=0; i<6; ++i) {
	init(powers) = 0, 0, 0, 0, 0, 0;
	powers[i] = j;
	Polynomial p(powers, CoeffOne);
	CPPUNIT_ASSERT(p.degree() == j);
      }
    }
  }
  void testCallMul(void) {
    std::vector< Power > powers;
    for (Power i=0; i<10; ++i) { init(powers) += i; }
    CPPUNIT_ASSERT(powers.size() == 10);
    Polynomial p(powers);
    Polynomial p_accumulate = Polynomial::One(10);
    Polynomial p_accumulate2 = Polynomial::One(10);
    std::vector< Coeff > x;
    init(x) = Coeff(0.5), Coeff(-0.2), Coeff(0.3), Coeff(1.67), Coeff(-2.001), Coeff(0.1), Coeff(0.2), Coeff(0.3), Coeff(0.4), Coeff(-0.9);
    Coeff prod(1.0);
    for (Index i=0; i<10; ++i) {
      std::vector< Power > one;
      init(one) = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      CPPUNIT_ASSERT(one.size() == 10);
      one[i] = powers[i];
      Polynomial p_one(powers);
      prod *= p_one(x);
      p_accumulate = p_accumulate * p_one;
      p_accumulate2 *= p_one;
    }
    const Real tol = 1.0e-14;
    CPPUNIT_ASSERT(fabs(p(x) - prod) < tol);
    CPPUNIT_ASSERT(fabs(p_accumulate(x) - prod) < tol);
    CPPUNIT_ASSERT(fabs(p_accumulate2(x) - prod) < tol);
  }
  void testDiffSum(void) {
    for (Index vars=1; vars<100; ++vars) {
      std::vector< Power > powers;
      for (Index i=0; i<vars; ++i) {
	init(powers) += i+1;
      }
      CPPUNIT_ASSERT(powers.size() == vars);
      MapPowers triangle(powers);
      Polynomial t(triangle, Coeff(-1.0));
      Polynomial u(vars);
      for (Index i=0; i<vars; ++i) {
	u += t.diff(i);
      }
      CPPUNIT_ASSERT(u.getNumVars() == vars);
      CPPUNIT_ASSERT(u.getNumTerms() == vars);
      std::vector< Coeff > ones;
      for (Index i=0; i<vars; ++i) {
	init(ones) += CoeffOne;
      }
      CPPUNIT_ASSERT(ones.size() == vars);
      CPPUNIT_ASSERT(u(ones) == -Coeff(vars*(vars+1)/2));
    }
  }
  void testIsConstantZero(void) {
    Polynomial p(11);
    CPPUNIT_ASSERT(p.isConstant());
  }
  void testIsConstantOneViaFactory(void) {
    Polynomial p = Coeff(-4.2)*Polynomial::One(11);
    CPPUNIT_ASSERT(p.isConstant());
  }
  void testIsConstantTooManyTerms(void) {
    Polynomial pol(11);
    Powers p(11);
    p.setPower(3, 1);
    pol += Polynomial(p);
    Powers q(11);
    q.setPower(7, 1);
    pol -= Coeff(4.2)*Polynomial(q);
    pol = pol + Coeff(3.112)*Polynomial::One(11);
    CPPUNIT_ASSERT(!(pol.isConstant()));
  }
  void testIsConstantSingleNonConstTerm(void) {
    const Index numVars = 11;
    for (Index i = 0; i < numVars; ++i) {
      Polynomial p(numVars);
      CPPUNIT_ASSERT(p.isZero());
      CPPUNIT_ASSERT(p.isConstant());
      Powers pows(11);
      pows.setPower(i, 1);
      p += Polynomial(pows);
      CPPUNIT_ASSERT(!(p.isZero()));
      CPPUNIT_ASSERT(!(p.isConstant()));
    }
  }
}; //PolynomialTest

#endif //POLYNOMIAL_TEST_H
