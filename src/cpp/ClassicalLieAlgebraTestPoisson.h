//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
//----------------------------------------------------------------------

#ifndef CLASSICAL_LIE_ALGEBRA_TEST_POISSON_H
#define CLASSICAL_LIE_ALGEBRA_TEST_POISSON_H

//system
#include <utility> //for pair

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"
#include "LieAlgebraBase.h"
#include "ClassicalLieAlgebra.h"

static const size_t PolyPoissonNumCases = 8;
static const Index PolyPoissonMaxDof = 4;
static const Real PolyPoissonTolerance = 1.0e-12;

class ClassicalLieAlgebraTestPoisson : public CppUnit::TestFixture {

 public:

  //The following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(ClassicalLieAlgebraTestPoisson);

  //Poisson bracket:
  CPPUNIT_TEST(testPoissonLinearInFirstMul);
  CPPUNIT_TEST(testPoissonLinearInSecondMul);
  CPPUNIT_TEST(testPoissonLinearInFirstAdd);
  CPPUNIT_TEST(testPoissonLinearInSecondAdd);
  CPPUNIT_TEST(testPoissonAntiSymmetric);
  CPPUNIT_TEST(testPoissonJacobiIdentity);
  CPPUNIT_TEST(testPoissonSelfIsZero);

  CPPUNIT_TEST_EXCEPTION(testPoissonChecksEvenNumVars,
			 LieAlgebraSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testPoissonChecksSameNumVars,
                         LieAlgebraSizeMismatchError);

  CPPUNIT_TEST_SUITE_END();

  void testPoissonLinearInFirstMul(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(Coeff(2.0)*a, b);
      const Polynomial rhs = Coeff(2.0)*alg.lieBracket(a, b);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void testPoissonLinearInSecondMul(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, Coeff(-3.0)*b);
      const Polynomial rhs = Coeff(-3.0)*alg.lieBracket(a, b);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void testPoissonLinearInFirstAdd(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial c = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a + b, c);
      const Polynomial rhs = alg.lieBracket(a, c) + alg.lieBracket(b, c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void testPoissonLinearInSecondAdd(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial c = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, b + c);
      const Polynomial rhs = alg.lieBracket(a, b) + alg.lieBracket(a, c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void testPoissonAntiSymmetric(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial c = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, b);
      const Polynomial rhs = -alg.lieBracket(b, a);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }

  void testPoissonJacobiIdentity(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars());
      const Polynomial c = randomPoly(alg.getNumVars());
      Polynomial lhs(alg.zero());
      lhs += alg.lieBracket(a, alg.lieBracket(b, c));
      lhs += alg.lieBracket(b, alg.lieBracket(c, a));
      lhs += alg.lieBracket(c, alg.lieBracket(a, b));
      CPPUNIT_ASSERT(lhs.lInfinityNorm() < 1.0e-8); //matters.
    }
  }

  void testPoissonChecksEvenNumVars(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(25);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars() + 1);
      const Polynomial b = randomPoly(alg.getNumVars() + 1);
      const Polynomial lhs = alg.lieBracket(a, b);
    }
  }

  void testPoissonChecksSameNumVars(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(25);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial b = randomPoly(alg.getNumVars() + 2);
      const Polynomial lhs = alg.lieBracket(a, b);
    }
  }
 
  void testPoissonSelfIsZero(void) {
    for (size_t n = 0; n < PolyPoissonNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyPoissonMaxDof);
      const ClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPoly(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, a);
      CPPUNIT_ASSERT(lhs.lInfinityNorm() < 1.0e-12); //matters?
    }
  }

}; //ClassicalLieAlgebraTestPoisson

#endif //CLASSICAL_LIE_ALGEBRA_TEST_POISSON_H
