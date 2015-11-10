//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: SemiClassicalLieAlgebraTestMoyal
//
// PURPOSE:
//
// Test the Lie (Moyal) bracket for the semi-classical Lie algebra.
//
// NOTES:
//
// Note that the numerical tolerances here need to be somewhat looser
// than for the corresponding tests in the classical version;  this is
// especially true when using double precision arithmetic, rather than
// the GNU Multi-Precision library.
//
// TO DO:
//
// I need to add the commutator relations.
//
// At present, there seems to be a problem with the Jacobi identity.
//
//----------------------------------------------------------------------

#ifndef SEMI_CLASSICAL_LIE_ALGEBRA_TEST_MOYAL_H
#define SEMI_CLASSICAL_LIE_ALGEBRA_TEST_MOYAL_H

//system
#include <utility> //for pair

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "Random.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"
#include "LieAlgebraBase.h"
#include "SemiClassicalLieAlgebra.h"

static const size_t PolyMoyalNumCases = 4;
static const size_t PolyMoyalNumConstantsPerRelation = 3;
static const Index PolyMoyalMaxDof = 3;
static const Power PolyMoyalMaxDegree = 8;
static const Index PolyMoyalMaxNumTerms = 16;

Polynomial stripHBar(const Polynomial& poly) {
  Polynomial result(poly.getNumVars());
  const PowersToCoeffMap& pToC = poly.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  for (pc = pToC.begin(); pc != pToC.end(); ++pc) {
    Powers powers = pc->first;
    if (powers[poly.getNumVars() - 1] > 0)
      powers.setPower(poly.getNumVars() - 1, 0);
    result.setMonomial(powers, pc->second);
  }
  return result;
}

Polynomial randomPolyMoyal(const Index numVars) {
  Polynomial p = randomPoly(numVars, numVars, PolyMoyalMaxDegree, PolyMoyalMaxNumTerms);
  //return stripHBar(p);
  return p;
}

class SemiClassicalLieAlgebraTestMoyal : public CppUnit::TestFixture {

 public:

  //The following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(SemiClassicalLieAlgebraTestMoyal);

  //Moyal bracket:
/*
  CPPUNIT_TEST(testMoyalAnyZeroIsZero);
  CPPUNIT_TEST(testMoyalLinearInFirstMulExample);
  CPPUNIT_TEST(testMoyalLinearInFirstMul);
  CPPUNIT_TEST(testMoyalLinearInSecondMulExample);
  CPPUNIT_TEST(testMoyalLinearInSecondMul);
  CPPUNIT_TEST(testMoyalLinearInFirstAdd);
  CPPUNIT_TEST(testMoyalLinearInSecondAdd);
  CPPUNIT_TEST(testMoyalAntiSymmetric);
  CPPUNIT_TEST(testMoyalSelfIsZero);
*/
  CPPUNIT_TEST(testMoyalJacobiIdentity);

  //range testing
  CPPUNIT_TEST_EXCEPTION(testMoyalChecksEvenNumVars, LieAlgebraSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testMoyalChecksSameNumVars, LieAlgebraSizeMismatchError);

  CPPUNIT_TEST_SUITE_END();

  void testMoyalAnyZeroIsZero(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = alg.zero();
      const Polynomial lhs = alg.lieBracket(a, b);
      const Polynomial rhs = alg.zero();
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalLinearInFirstMulExample(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(Coeff(2.0)*a, b);
      const Polynomial rhs = Coeff(2.0)*alg.lieBracket(a, b);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalLinearInFirstMul(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      for (size_t i = 0; i < PolyMoyalNumConstantsPerRelation; ++i) {
        Coeff c = randomCoeff(5.0);
        const Polynomial lhs = alg.lieBracket(c*a, b);
        const Polynomial rhs = c*alg.lieBracket(a, b);
        CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
      }
    }
  }
  void testMoyalLinearInSecondMulExample(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, Coeff(-3.0)*b);
      const Polynomial rhs = Coeff(-3.0)*alg.lieBracket(a, b);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalLinearInSecondMul(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      for (size_t i = 0; i < PolyMoyalNumConstantsPerRelation; ++i) {
        Coeff c = randomCoeff(5.0);
        const Polynomial lhs = alg.lieBracket(a, c*b);
        const Polynomial rhs = c*alg.lieBracket(a, b);
        CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
      }
    }
  }
  void testMoyalLinearInFirstAdd(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial c = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a + b, c);
      const Polynomial rhs = alg.lieBracket(a, c) + alg.lieBracket(b, c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalLinearInSecondAdd(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial c = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, b + c);
      const Polynomial rhs = alg.lieBracket(a, b) + alg.lieBracket(a, c);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalAntiSymmetric(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial c = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, b);
      const Polynomial rhs = -alg.lieBracket(b, a);
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < poly_ring_tol);
    }
  }
  void testMoyalSelfIsZero(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial lhs = alg.lieBracket(a, a);
      CPPUNIT_ASSERT(lhs.lInfinityNorm() < 1.0e-12); //matters?
    }
  }
  void testMoyalJacobiIdentity(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(PolyMoyalMaxDof);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars());
      const Polynomial c = randomPolyMoyal(alg.getNumVars());
      Polynomial lhs(alg.zero());
      lhs += alg.lieBracket(a, alg.lieBracket(b, c));
      lhs += alg.lieBracket(b, alg.lieBracket(c, a));
      lhs += alg.lieBracket(c, alg.lieBracket(a, b));
      Real lInfinityNorm = lhs.lInfinityNorm();
      if (lInfinityNorm > 1.0e-10) {
        std::cerr << "a" << std::endl << "= " << a;
        std::cerr << "b" << std::endl << "= " << b;
        std::cerr << "c" << std::endl << "= " << c;
        std::cerr << "lhs" << std::endl << "= " << lhs;
        std::cerr << "Jacobi identity l-infinity norm " << lInfinityNorm << std::endl;
      }
      CPPUNIT_ASSERT(lInfinityNorm < 1.0e-10); //matters.
    }
  }
  void testMoyalChecksEvenNumVars(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = randomIndex(25);
      const Index numVars = 1 + (1 + 2*dof);
      const Polynomial a = randomPolyMoyal(numVars);
      const Polynomial b = randomPolyMoyal(numVars);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial lhs = alg.lieBracket(a, b);
    }
  }
  void testMoyalChecksSameNumVars(void) {
    for (size_t n = 0; n < PolyMoyalNumCases; ++n) {
      const Index dof = 1 + randomIndex(25);
      const SemiClassicalLieAlgebra alg(dof);
      const Polynomial a = randomPolyMoyal(alg.getNumVars());
      const Polynomial b = randomPolyMoyal(alg.getNumVars() + 2);
      const Polynomial lhs = alg.lieBracket(a, b);
    }
  }
 
}; //SemiClassicalLieAlgebraTestMoyal

#endif //SEMI_CLASSICAL_LIE_ALGEBRA_TEST_MOYAL_H
