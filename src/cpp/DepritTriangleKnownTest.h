//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangleKnownTest
//
// PURPOSE:
//
// Tests for DepritTriangle with a known generator.
//
// NOTES:
//
// In particular the tests for invariance of the generator must a
// handled carefully: to transform a generator term via a Deprit
// triangle, we require the generator term of grade one higher.  This
// means that one cannot check the highest grade term for invariance?
//
// TODO:
//
// This seems odd and I will have another look at it.
//
// (20050228) We need a cubic and higher generator, where each term of
// (20050228) grade g (>=3) is correctly scaled via factorial(g - 3).
//
//----------------------------------------------------------------------

#ifndef DEPRIT_TRIANGLE_KNOWN_TEST_H
#define DEPRIT_TRIANGLE_KNOWN_TEST_H

//standard headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "DefaultNumericalPrecision.h"
#include "DepritTriangleKnown.h"

const size_t NumDepritTriangleKnownTestCases(16);
const Power NumDepritTriangleKnownTestRows(5);

std::vector< Polynomial > randInnerGenTerms(const LieAlgebraBase& algebra) {
  const Polynomial w(randomPoly(algebra.getNumVars()));
  std::vector< Polynomial > wInnerTerms;
  wInnerTerms.push_back(algebra.zero());
  for (Power j = 1; j < NumDepritTriangleKnownTestRows; ++j) {
    const Coeff preFactor(factorial(j - 1));
    const Polynomial wInnerJ(algebra.isoGradePart(w, j + 2) * preFactor);
    wInnerTerms.push_back(wInnerJ);
  }
  return wInnerTerms;
}

class DepritTriangleKnownTest : public CppUnit::TestFixture {
 public:
  CPPUNIT_TEST_SUITE(DepritTriangleKnownTest);
  //
  //CPPUNIT_TEST(testFixtureIncluded);
  //
  //CPPUNIT_TEST(testForwardGeneratorInvariant);
  //CPPUNIT_TEST(testReverseGeneratorInvariant);
  //
  CPPUNIT_TEST(testReverseTransformZeroIsZeroExample);
  CPPUNIT_TEST(testReverseTransformZeroIsZero);
  CPPUNIT_TEST(testReverseTransformLinearAdd);
  CPPUNIT_TEST(testReverseTransformLinearMul);
  CPPUNIT_TEST(testForwardTransformZeroIsZeroExample);
  CPPUNIT_TEST(testForwardTransformZeroIsZero);
  CPPUNIT_TEST(testForwardTransformLinearAdd);
  CPPUNIT_TEST(testForwardTransformLinearMul);
  //
  CPPUNIT_TEST_SUITE_END();
 protected:
  //void testFixtureIncluded(void) {
  //  CPPUNIT_ASSERT(false);
  //}
  void testReverseTransformZeroIsZeroExample(void) {
    const Index dof(2);
    const ClassicalLieAlgebra algebra(dof);
    //a generator beginning with quadratic terms
    const Polynomial w0(algebra.zero());
    const Polynomial w1(algebra.q(0, 3));
    const Polynomial w2(algebra.p(1, 4));
    std::vector< Polynomial > wInnerTerms;
    wInnerTerms.push_back(w0);
    wInnerTerms.push_back(w1);
    wInnerTerms.push_back(w2);
    //transform zero polynomial
    DepritTriangleKnown transform(algebra, wInnerTerms);
    Polynomial result(algebra.zero());
    for (Power i = 0; i < 3; ++i) {
      transform.reverse(algebra.zero(), result);
      CPPUNIT_ASSERT(result.isZero());
    }
  }
  void testReverseTransformZeroIsZero(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //transform zero
      DepritTriangleKnown transform(algebra, wInnerTerms);
      Polynomial result(algebra.zero());
      for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
        transform.reverse(algebra.zero(), result);
        CPPUNIT_ASSERT(result.isZero());
      }
    }
  }
  void testReverseTransformLinearMul(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //linearity by multiplication
      for (size_t m = 0; m < NumDepritTriangleKnownTestCases; ++m) {
        const Polynomial f(randomPoly(algebra.getNumVars()));
        for (size_t k = 0; k < NumDepritTriangleKnownTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          DepritTriangleKnown transformLhs(algebra, wInnerTerms);
          DepritTriangleKnown transformRhs(algebra, wInnerTerms);
          Polynomial resultLhs(algebra.zero());
          Polynomial resultRhs(algebra.zero());
          for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
            const Polynomial term(algebra.isoGradePart(f, row + 2));
            transformLhs.reverse(term, resultLhs);
            transformRhs.reverse(term * c, resultRhs);
            const Real error((resultLhs * c - resultRhs).lInfinityNorm());
            CPPUNIT_ASSERT(error < 1.0e-12);
          }
        }
      }
    }
  }
  void testReverseTransformLinearAdd(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //linearity by addition
      for (size_t m = 0; m < NumDepritTriangleKnownTestCases; ++m) {
        const Polynomial polA(randomPoly(algebra.getNumVars()));
        const Polynomial polB(randomPoly(algebra.getNumVars()));
        const Polynomial polC(polA + polB);
        for (size_t k = 0; k < NumDepritTriangleKnownTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          DepritTriangleKnown transformA(algebra, wInnerTerms);
          DepritTriangleKnown transformB(algebra, wInnerTerms);
          DepritTriangleKnown transformC(algebra, wInnerTerms);
          Polynomial newA(algebra.zero());
          Polynomial newB(algebra.zero());
          Polynomial newC(algebra.zero());
          for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
            const Polynomial termA(algebra.isoGradePart(polA, row + 2));
            const Polynomial termB(algebra.isoGradePart(polB, row + 2));
            const Polynomial termC(algebra.isoGradePart(polC, row + 2));
            transformA.reverse(termA, newA);
            transformB.reverse(termB, newB);
            transformC.reverse(termC, newC);
            CPPUNIT_ASSERT(((newA + newB) - newC).lInfinityNorm() < 1.0e-12);
          }
        }
      }
    }
  }
  void testForwardTransformZeroIsZeroExample(void) {
    const Index dof(2);
    const ClassicalLieAlgebra algebra(dof);
    //a generator beginning with zero quadratic terms
    const Polynomial w0(algebra.zero());
    const Polynomial w1(algebra.q(0, 3));
    const Polynomial w2(algebra.p(1, 4));
    const Polynomial w3(algebra.q(0, 2) * algebra.p(0, 2) * algebra.p(1, 1));
    std::vector< Polynomial > wInnerTerms;
    wInnerTerms.push_back(w0);
    wInnerTerms.push_back(w1);
    wInnerTerms.push_back(w2);
    wInnerTerms.push_back(w3);
    //transform zero polynomial
    DepritTriangleKnown transform(algebra, wInnerTerms);
    Polynomial result(algebra.zero());
    for (Power i = 0; i < 4; ++i) {
      transform.forward(algebra.zero(), result);
      CPPUNIT_ASSERT(result.isZero());
    }
  }
  void testForwardTransformZeroIsZero(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //transform zero
      DepritTriangleKnown transform(algebra, wInnerTerms);
      Polynomial result(algebra.zero());
      for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
        transform.forward(algebra.zero(), result);
        CPPUNIT_ASSERT(result.isZero());
      }
    }
  }
  void testForwardTransformLinearMul(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //linearity by multiplication
      for (size_t m = 0; m < NumDepritTriangleKnownTestCases; ++m) {
        const Polynomial f(randomPoly(algebra.getNumVars()));
        for (size_t k = 0; k < NumDepritTriangleKnownTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          DepritTriangleKnown transformLhs(algebra, wInnerTerms);
          DepritTriangleKnown transformRhs(algebra, wInnerTerms);
          Polynomial resultLhs(algebra.zero());
          Polynomial resultRhs(algebra.zero());
          for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
            const Polynomial term(algebra.isoGradePart(f, row + 2));
            transformLhs.forward(term, resultLhs);
            transformRhs.forward(term * c, resultRhs);
            const Real error((resultLhs * c - resultRhs).lInfinityNorm());
            CPPUNIT_ASSERT(error < 1.0e-12);
          }
        }
      }
    }
  }
  void testForwardTransformLinearAdd(void) {
    for (size_t n = 0; n < NumDepritTriangleKnownTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randInnerGenTerms(algebra));
      //linearity by addition
      for (size_t m = 0; m < NumDepritTriangleKnownTestCases; ++m) {
        const Polynomial polA(randomPoly(algebra.getNumVars()));
        const Polynomial polB(randomPoly(algebra.getNumVars()));
        const Polynomial polC(polA + polB);
        for (size_t k = 0; k < NumDepritTriangleKnownTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          DepritTriangleKnown transformA(algebra, wInnerTerms);
          DepritTriangleKnown transformB(algebra, wInnerTerms);
          DepritTriangleKnown transformC(algebra, wInnerTerms);
          Polynomial newA(algebra.zero());
          Polynomial newB(algebra.zero());
          Polynomial newC(algebra.zero());
          for (Power row = 0; row < NumDepritTriangleKnownTestRows; ++row) {
            const Polynomial termA(algebra.isoGradePart(polA, row + 2));
            const Polynomial termB(algebra.isoGradePart(polB, row + 2));
            const Polynomial termC(algebra.isoGradePart(polC, row + 2));
            transformA.forward(termA, newA);
            transformB.forward(termB, newB);
            transformC.forward(termC, newC);
            CPPUNIT_ASSERT(((newA + newB) - newC).lInfinityNorm() < 1.0e-12);
          }
        }
      }
    }
  }
//   void testSomeExample(void) {
//     const DefaultNumericalPrecision precision(256);
//     std::cerr << precision << std::endl;

//     //construct example algebra
//     const Index dof(3);
//     const ClassicalLieAlgebra algebra(dof);

//     //construct quadratic part of hamiltonian
//     Polynomial h2(algebra.zero());
//     std::vector< Coeff > freqs;
//     freqs.push_back(Coeff(-1));
//     freqs.push_back(Coeff(0, 3));
//     freqs.push_back(Coeff(0.123456));
//     for (size_t t = 0; t < freqs.size(); ++t)
//       h2 += freqs[t] * (algebra.q(t) * algebra.p(t));

//     //prepare some polynomials to hold results
//     Polynomial kTerm(algebra.zero());
//     Polynomial wTerm(algebra.zero());

//     //compute first deprit triangle row
//     DepritTriangle triangle(algebra, h2);
//     triangle.computeNormalFormAndGenerator(h2, kTerm, wTerm);
//     std::cerr << kTerm << std::endl;
//     std::cerr << wTerm << std::endl;
//     CPPUNIT_ASSERT(kTerm == h2);
//     CPPUNIT_ASSERT(wTerm.isZero());

//     //construct cubic part of hamiltonian
//     Polynomial h3(algebra.zero());
//     h3 += Coeff(4.0) * algebra.q(0, 2) * algebra.p(1);
//     h3 += CoeffJ * algebra.q(1) * algebra.p(1) * algebra.q(0);

//     //compute second deprit triangle row
//     triangle.computeNormalFormAndGenerator(h3, kTerm, wTerm);
//     std::cerr << kTerm << std::endl;
//     std::cerr << wTerm << std::endl;

//     //construct a quartic part
//     Polynomial h4(algebra.zero());
//     h4 += algebra.q(2) * h3;
//     h4 += h2 * h2;
//     h4 += algebra.q(0) * algebra.q(1) * h2;

//     //compute second deprit triangle row
//     triangle.computeNormalFormAndGenerator(h4, kTerm, wTerm);
//     std::cerr << kTerm << std::endl;
//     std::cerr << wTerm << std::endl;

//     //compute some dummy rows
//     for (size_t n = 0; n < 3; ++n) {
//       triangle.computeNormalFormAndGenerator(algebra.zero(), kTerm, wTerm);
//       std::cerr << kTerm << std::endl;
//       std::cerr << wTerm << std::endl;
//     }

//     //try a coordinate transform
//     LieTransform transform0(algebra, triangle.getInnerGeneratorTerms());
//     Polynomial q0Diag(algebra.q(0));
//     Polynomial q0Norm(algebra.zero());
//     std::cerr << q0Diag << std::endl;
//     CPPUNIT_ASSERT(algebra.isIsoGrade(q0Norm));
//     for (size_t n = 0; n < triangle.getCurrentGrade() - 1; ++n) {
//       if (n == 0)
//         transform0.reverseLieTransform(q0Diag, q0Norm);
//       else
//         transform0.reverseLieTransform(algebra.zero(), q0Norm);
//       std::cerr << q0Norm << std::endl;
//       CPPUNIT_ASSERT(algebra.isIsoGrade(q0Norm));
//     }
//   }
 private:
}; //DepritTriangleKnownTest

#endif //DEPRIT_TRIANGLE_KNOWN_TEST_H
