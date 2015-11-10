//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef POLYNOMIAL_TEST_RING_H
#define POLYNOMIAL_TEST_RING_H

//system
#include <utility> //for pair

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"

static const size_t poly_ring_n_cases = 32;
static const Index poly_ring_max_vars = 50;
static const Real poly_ring_tol(1.0e-12);

class PolynomialTestRing : public CppUnit::TestFixture {

 public:

  //The following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(PolynomialTestRing);

  //Basic assign vs. copy:
  CPPUNIT_TEST(test_copy_assign);

  //Arithmetic:
  CPPUNIT_TEST(test_example_call_same_accumulate_monomials);
  CPPUNIT_TEST(test_example_call_same_accumulate_decumulate);
  CPPUNIT_TEST(test_iadd_self);
  CPPUNIT_TEST(test_isub_self);
  CPPUNIT_TEST(test_mul_zero);
  CPPUNIT_TEST(test_mul_example);
  CPPUNIT_TEST(test_call_add);

  //Ring; equality:
  CPPUNIT_TEST(test_ring_sub_for_use_in_eq);

  //Ring; addition:
  CPPUNIT_TEST(test_ring_commutative_add);
  CPPUNIT_TEST(test_ring_associative_add);
  CPPUNIT_TEST(test_ring_identity_add);
  CPPUNIT_TEST(test_ring_inverse_add);
  CPPUNIT_TEST(test_ring_add_sub);
  CPPUNIT_TEST(test_ring_zero_add_sub);
  CPPUNIT_TEST(test_ring_add_neg_sub);
  CPPUNIT_TEST(test_ring_iadd);
  CPPUNIT_TEST(test_ring_iadd_self);
  CPPUNIT_TEST(test_ring_isub);
  CPPUNIT_TEST(test_ring_isub_self);

  //Ring; multiplication:
  CPPUNIT_TEST(test_ring_commutative_mul);
  CPPUNIT_TEST(test_ring_associative_mul);
  CPPUNIT_TEST(test_ring_identity_mul);
  CPPUNIT_TEST(test_ring_absorb_mul);
  CPPUNIT_TEST(test_ring_imul);
  CPPUNIT_TEST(test_ring_imul_self);

  //Ring; combined:
  CPPUNIT_TEST(test_ring_distributive_left);
  CPPUNIT_TEST(test_ring_distributive_right);
  CPPUNIT_TEST(test_ring_scalar_multiply);
  CPPUNIT_TEST(test_ring_scalar_multiply_twice);
  CPPUNIT_TEST(test_ring_scalar_multiply_vs_unary_minus);
  CPPUNIT_TEST(test_ring_zero_neg);

  CPPUNIT_TEST_SUITE_END();

private:

protected:

  // void setUp(void) { }

  // void tearDown(void) { }

  void test_copy_assign(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1 + randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs(a);
      Polynomial rhs(a.getNumVars());
      rhs = a;
      CPPUNIT_ASSERT((rhs - lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_example_call_same_accumulate_monomials(void) {
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
      v[i] = Coeff(i);
    }
    
    Polynomial qq(p00);
    qq += p11;
    qq += p22;
    qq += p33;
    CPPUNIT_ASSERT(p == qq);
    CPPUNIT_ASSERT(p(v) == qq(v));
  }

  void test_example_call_same_accumulate_decumulate(void) {
    Polynomial p(3);
    MapPowers pows(3);

    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    p.setMonomial(pows, CoeffOne);
    CPPUNIT_ASSERT(p[pows] == CoeffOne);
    Polynomial p00(pows, CoeffOne);

    pows.setPower(0, 0);
    pows.setPower(1, 7);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(-17.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(-17.0));
    Polynomial p11(pows, Coeff(-17.0));

    pows.setPower(0, 0);
    pows.setPower(1, 0);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(42.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(42.0));
    Polynomial p22(pows, Coeff(42.0));

    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 5);
    p.setMonomial(pows, CoeffZero);
    CPPUNIT_ASSERT(p[pows] == CoeffZero);
    Polynomial p33(pows, CoeffZero);

    CPPUNIT_ASSERT(p.getNumTerms() == 3);

    Polynomial qq(p00);
    qq += p11;
    qq += p22;
    qq += p33;
    CPPUNIT_ASSERT(p == qq);

    p -= p00;
    p -= p11;
    p -= p22;
    p -= p33;
    CPPUNIT_ASSERT(p.lInfinityNorm() < 1.0e-14); //zero
    //WARNING: p should have zero terms, but might have small ones?!
    //maybe we should set a global truncation?
  }

  void test_iadd_self(void) {
    Polynomial p(3);
    MapPowers pows(3);

    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    p.setMonomial(pows, CoeffOne);
    CPPUNIT_ASSERT(p[pows] == CoeffOne);
    MapPowers pows0(pows);

    pows.setPower(0, 0);
    pows.setPower(1, 7);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(-17.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(-17.0));
    MapPowers pows1(pows);

    Polynomial q(p);
    q += q;
    CPPUNIT_ASSERT(q[pows0] == Coeff(2.0)*p[pows0]);
    CPPUNIT_ASSERT(q[pows1] == Coeff(2.0)*p[pows1]);
  }

  void test_isub_self(void) {
    Polynomial p(3);
    MapPowers pows(3);

    pows.setPower(0, 1);
    pows.setPower(1, 2);
    pows.setPower(2, 3);
    p.setMonomial(pows, CoeffOne);
    CPPUNIT_ASSERT(p[pows] == CoeffOne);
    MapPowers pows0(pows);

    pows.setPower(0, 0);
    pows.setPower(1, 7);
    pows.setPower(2, 0);
    p.setMonomial(pows, Coeff(-17.0));
    CPPUNIT_ASSERT(p[pows] == Coeff(-17.0));
    MapPowers pows1(pows);

    Polynomial q(p);
    q -= q;
    CPPUNIT_ASSERT(q.getNumTerms() == 0);
  }

  void test_mul_zero(void) {
    std::vector< Power > powers(5);
    init(powers) = 2, 2, 3, 4, 5;
    Polynomial t(powers);
    CPPUNIT_ASSERT(t.getNumTerms() == 1);
    init(powers) = 1, 2, 3, 4, 5;
    Polynomial tt(powers);
    CPPUNIT_ASSERT(tt.getNumTerms() == 1);
    t += tt;
    CPPUNIT_ASSERT(t.getNumTerms() == 2);
    init(powers) = 5, 2, 8, 4, 5;
    Polynomial zero = CoeffZero*Polynomial(powers);
    CPPUNIT_ASSERT(zero.isZero());
    Polynomial t0 = t*zero;
    CPPUNIT_ASSERT(t0.isZero());
    Polynomial r1 = zero*t;
    CPPUNIT_ASSERT(r1.isZero());
  }

  void test_mul_example(void) {
    std::vector< Power > powers(2);
    init(powers) = 1, 0;
    Polynomial x(powers);
    CPPUNIT_ASSERT(x == Polynomial::CoordinateMonomial(2, 0));
    init(powers) = 0, 1;
    Polynomial y(powers);
    CPPUNIT_ASSERT(y == Polynomial::CoordinateMonomial(2, 1));
    Polynomial xy = x*y;
    Polynomial yx = y*x;
    CPPUNIT_ASSERT(xy == yx);
    CPPUNIT_ASSERT(xy.getNumTerms() == 1);
    Polynomial t = xy;
    t += Coeff(2)*x;
    t += y*Coeff(3);
    t += xy;
    t *= (xy+Polynomial::One(2));
    CPPUNIT_ASSERT(t == (Coeff(2)*x*y + Coeff(2)*x + Coeff(3)*y + Coeff(2)*x*x*y*y + Coeff(2)*x*x*y + Coeff(3)*x*y*y));
  }

  void test_call_add(void) {
    std::vector< Power > powers;
    for (Power i = 0; i < 10; ++i) { init(powers) += i; }
    CPPUNIT_ASSERT(powers.size() == 10);
    Polynomial p_as_one(powers);
    std::vector< Coeff > x;
    init(x) = Coeff(0.5), Coeff(-0.2), Coeff(0.3), Coeff(1.67), Coeff(-2.001), Coeff(0.1), Coeff(0.2), Coeff(0.3), Coeff(0.4), Coeff(-0.9);
    CPPUNIT_ASSERT(x.size() == 10);
    Coeff coeff = CoeffOne/Coeff(x.size());
    Coeff total = CoeffZero;
    Polynomial p_as_many(x.size());
    Polynomial p_accumulate(x.size());
    Polynomial p_accumulate2(x.size());
    for (Index i=0; i<x.size(); ++i) {
      MapPowers m_one(powers);
      Polynomial p_one(powers, coeff);
      total += p_one(x);
      p_as_many.setMonomial(m_one, coeff);
      p_accumulate += p_one;
      p_accumulate2 = p_one + p_accumulate2;
    }
    CPPUNIT_ASSERT(fabs(p_as_one(x) - total) < poly_ring_tol);
    CPPUNIT_ASSERT(fabs(p_as_many(x) - total) < poly_ring_tol);
    CPPUNIT_ASSERT(fabs(p_accumulate(x) - total) < poly_ring_tol);
    CPPUNIT_ASSERT(fabs(p_accumulate2(x) - total) < poly_ring_tol);
  }

  void test_ring_sub_for_use_in_eq(void) {
    for (size_t n_case=0; n_case<poly_ring_n_cases; ++n_case) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z = Polynomial(vars);
      Polynomial neg_a = z-a;
      CPPUNIT_ASSERT(a+neg_a == z);
    }
  }

  void test_ring_associative_add(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = (a+b) + c;
      Polynomial rhs = a + (b+c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_associative_mul(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = (a*b) * c;
      Polynomial rhs = a * (b*c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_commutative_add(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial lhs = a + b;
      Polynomial rhs = b + a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_commutative_mul(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial lhs = a * b;
      Polynomial rhs = b * a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_identity_add(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = a + z;
      Polynomial rhs = a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_identity_mul(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial one = Polynomial::One(vars);
      Polynomial lhs = a * one;
      Polynomial rhs = a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_absorb_mul(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = a * z;
      Polynomial rhs = z;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_inverse_add(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = a + (-a);
      Polynomial rhs = z;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_distributive_left(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = a * (b + c);
      Polynomial rhs = (a*b) + (a*c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_distributive_right(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = (a + b) * c;
      Polynomial rhs = (a*c) + (b*c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_scalar_multiply(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = Coeff(0.3)*a + a*Coeff(0.7);
      Polynomial rhs = a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_scalar_multiply_twice(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = Coeff(2.0)*a;
      Polynomial rhs = a + a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_scalar_multiply_vs_unary_minus(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = Coeff(-1.0)*a;
      Polynomial rhs = (-a);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_zero_neg(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = z - a;
      Polynomial rhs = Coeff(-1.0)*a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_add_sub(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = a - (b+c);
      Polynomial rhs = (a-b)-c;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_zero_add_sub(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = (z+a) - a;
      Polynomial rhs = z;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_add_neg_sub(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial c = randomPoly(vars);
      Polynomial lhs = -(a+b+c);
      Polynomial rhs = (-a) + (-b) + (-c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_iadd(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial lhs = a;
      lhs += b;
      Polynomial rhs = a+b;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_iadd_self(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = a;
      lhs += lhs;
      Polynomial rhs = Coeff(2.0)*a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_isub(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial lhs = a;
      lhs -= b;
      Polynomial rhs = a-b;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_isub_self(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial z(vars);
      Polynomial lhs = a;
      lhs -= lhs;
      Polynomial rhs = z;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_imul(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial b = randomPoly(vars);
      Polynomial lhs = a;
      lhs *= b;
      Polynomial rhs = a*b;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void test_ring_imul_self(void) {
    for (size_t n=0; n<poly_ring_n_cases; ++n) {
      Index vars = 1+randomIndex(poly_ring_max_vars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = a;
      lhs *= lhs;
      Polynomial rhs = a*a;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

}; //PolynomialTestRing

#endif //POLYNOMIAL_TEST_RING_H
