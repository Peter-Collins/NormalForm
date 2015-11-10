#ifndef COORDINATE_CHANGE_TEST_H
#define COORDINATE_CHANGE_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: CoordinateChangeTest
//
// PURPOSE:
//
// Tests for CoordinateChange.
//
// NOTES:
//
// [a] To ensure consistency of the forward and inverse Deprit
// triangles, we perform computations with randomly chosen generators,
// resulting in forward-inverse compositions with thousands of terms
// (for degrees of freedom 1..3, for example) and large infinity
// norms.  Testing the low order portions (up to the last computed
// Deprit triangle row) confirms the extremely high accuracy of these
// computations.
//
// [b] TODO: It is essential also to test consistency of both the
// direct and inverse coordinate changes with the original and
// normalised Hamiltonians from the normalisation Deprit triangle.
//
// [c(i)] Generators must have cubic and higher terms and their terms
// must be correctly factorially weighted.  We being with w_0 = 0
// (quadratic), after which for j >= 1 we have isoGradeParts(w, j +
// 2)*factorial(j - 1), so that this is equivalent to j >= 0
// isoGradeParts(w, j + 3)*factorial(j).  Note that the difference
// between factorial argument and grade is always 3 for the generator!
//
// [c(ii)] The factorial weight always matches the power of epsilon,
// which for w_1 (grade 3) is 0.
//
// [d] The previous incarnation of this file included tests which were
// based on the erroneous assumption that w is the generator of the
// flow in the usual sense, and therefore invariant under the
// coordinate mapping, and that one can compute the generator of the
// inverse by the negative of the generator.  As such, they are
// omitted from compilation.  I leave the source code in the file
// CoordinateChangeTestOld.h, as it might be of use for the testing of
// implementation of other normal form formats, for which those
// assumptions hold.
//
//----------------------------------------------------------------------

//standard headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "LieAlgebraBase.h"
#include "DefaultNumericalPrecision.h"
#include "CoordinateChange.h"

// NumCoordinateChangeTestCases was 8 and took ~8hrs PRC
const size_t NumCoordinateChangeTestCases(2);
const Power NumCoordinateChangeTestRows(6);

std::vector< Polynomial > randomInnerGenTerms(const LieAlgebraBase& algebra) {
  const Polynomial w(randomPoly(algebra.getNumVars()));
  std::vector< Polynomial > wInnerTerms;
  wInnerTerms.push_back(algebra.zero());
  for (Power j = 1; j < NumCoordinateChangeTestRows; ++j) {
    const Coeff preFactor(factorial(j - 1));
    const Polynomial wInnerJ(algebra.isoGradePart(w, j + 2) * preFactor);
    assert(algebra.isIsoGrade(wInnerJ, j + 2));
    wInnerTerms.push_back(wInnerJ);
  }
  return wInnerTerms;
}

class CoordinateChangeTest : public CppUnit::TestFixture {
 public:
  CPPUNIT_TEST_SUITE(CoordinateChangeTest);
  //
  //CPPUNIT_TEST(testFixtureIncluded);
  //
  CPPUNIT_TEST(testForwardReverseComposition);
  //
  CPPUNIT_TEST(testReverseTransformZeroIsZeroExample);
  CPPUNIT_TEST(testReverseTransformZeroIsZero);
  CPPUNIT_TEST(testReverseTransformIsNearIdentity);
  CPPUNIT_TEST(testReverseTransformLinearAdd);
  CPPUNIT_TEST(testReverseTransformLinearMulOnCoordinates);
  CPPUNIT_TEST(testReverseTransformLinearMul);
  CPPUNIT_TEST(testForwardTransformZeroIsZeroExample);
  CPPUNIT_TEST(testForwardTransformZeroIsZero);
  CPPUNIT_TEST(testForwardTransformIsNearIdentity);
  CPPUNIT_TEST(testForwardTransformLinearAdd);
  CPPUNIT_TEST(testForwardTransformLinearMulOnCoordinates);
  CPPUNIT_TEST(testForwardTransformLinearMul);
  //
  CPPUNIT_TEST_SUITE_END();
 protected:
  //void testFixtureIncluded(void) {
  //  CPPUNIT_ASSERT(false);
  //}
  void testForwardReverseComposition(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      //
      // we take a random generator for the deprit coordinate changes
      //
      const Index dof(1 + (n % 2));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //
      // we transform each coordinate, both forward and backward
      //
      std::vector< Polynomial > x;
      std::vector< Polynomial > yFromX;
      std::vector< Polynomial > xFromY;
      for (Index i = 0; i < algebra.getNumVars(); ++i) {
        std::cerr << ".";
        const Polynomial xI(algebra.coordinateMonomial(i));
        x.push_back(xI);
        CoordinateChange forwardI(algebra, wInnerTerms, false);
        CoordinateChange reverseI(algebra, wInnerTerms, true);
        Polynomial xIForward(algebra.zero());
        Polynomial xIReverse(algebra.zero());
        Polynomial result(algebra.zero());
        for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
          if (row == 0)
            forwardI.computeNextTerm(xI, result);
          else
            forwardI.computeNextTerm(algebra.zero(), result);
          xIForward += result;
          if (row == 0)
            reverseI.computeNextTerm(xI, result);
          else
            reverseI.computeNextTerm(algebra.zero(), result);
          xIReverse += result;
        }
        yFromX.push_back(xIForward);
        xFromY.push_back(xIReverse);
      }
      //
      // next, we compose forward and reverse maps in both directions
      // for each coordinate; this is a large computation!  indeed,
      // this is far more work than will ever get done in a real
      // normal form computation!
      //
      for (Index i = 0; i < algebra.getNumVars(); ++i) {
        std::cerr << ".";
        const Polynomial& xI(x.at(i));
        const Polynomial& yIFromX(yFromX.at(i));
        const Polynomial& xIFromY(xFromY.at(i));
        const Polynomial xIFR(xIFromY(yFromX));
        const Polynomial xIRF(yIFromX(xFromY));
        //std::cerr << xIFR.getNumTerms() << " " << xIFR.lInfinityNorm() << std::endl;
        //std::cerr << xIRF.getNumTerms() << " " << xIRF.lInfinityNorm() << std::endl;
        for (Power gra = 0; gra < NumCoordinateChangeTestRows + 1; ++gra) {
          CPPUNIT_ASSERT((algebra.isoGradePart(xIFR, gra) - algebra.isoGradePart(xI, gra)).lInfinityNorm() < 1.0e-12);
          CPPUNIT_ASSERT((algebra.isoGradePart(xIRF, gra) - algebra.isoGradePart(xI, gra)).lInfinityNorm() < 1.0e-12);
        }
      }
    }
  }
  void testReverseTransformZeroIsZeroExample(void) {
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
    CoordinateChange transform(algebra, wInnerTerms, true);
    Polynomial result(algebra.zero());
    for (Power i = 0; i < 4; ++i) {
      transform.computeNextTerm(algebra.zero(), result);
      CPPUNIT_ASSERT(result.isZero());
    }
  }
  void testReverseTransformZeroIsZero(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //transform zero
      CoordinateChange transform(algebra, wInnerTerms, true);
      Polynomial result(algebra.zero());
      for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
        transform.computeNextTerm(algebra.zero(), result);
        CPPUNIT_ASSERT(result.isZero());
      }
    }
  }
  void testReverseTransformIsNearIdentity(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //transform single coordinate one term
      for (Index i = 0; i < algebra.getNumVars(); ++i) {
        CoordinateChange transform(algebra, wInnerTerms, true);
        const Polynomial xI(algebra.coordinateMonomial(i));
        Polynomial result(algebra.zero());
        transform.computeNextTerm(xI, result);
        CPPUNIT_ASSERT(result == xI);
      }
    }
  }
  void testReverseTransformLinearMulOnCoordinates(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by multiplication
      for (Index s = 0; s < algebra.getNumVars(); ++s) {
        const Polynomial xS(algebra.coordinateMonomial(s));
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformLhs(algebra, wInnerTerms, true);
          CoordinateChange transformRhs(algebra, wInnerTerms, true);
          Polynomial termLhs(algebra.zero());
          Polynomial termRhs(algebra.zero());
          for (Power i = 0; i < NumCoordinateChangeTestRows; ++i) {
            if (i == 0) {
              transformLhs.computeNextTerm(xS, termLhs);
              transformRhs.computeNextTerm(xS * c, termRhs);
              CPPUNIT_ASSERT((termRhs - (xS * c)).lInfinityNorm() < 1.0e-12);
            }
            else {
              transformLhs.computeNextTerm(algebra.zero(), termLhs);
              transformRhs.computeNextTerm(algebra.zero(), termRhs);
            }
            CPPUNIT_ASSERT((termRhs - termLhs * c).lInfinityNorm() < 1.0e-12);
          }
        }
      }
    }
  }
  void testReverseTransformLinearMul(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by multiplication
      for (size_t m = 0; m < NumCoordinateChangeTestCases; ++m) {
        const Polynomial f(randomPoly(algebra.getNumVars()));
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformLhs(algebra, wInnerTerms, true);
          CoordinateChange transformRhs(algebra, wInnerTerms, true);
          Polynomial resultLhs(algebra.zero());
          Polynomial resultRhs(algebra.zero());
          for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
            const Polynomial term(algebra.isoGradePart(f, row + 1));
            transformLhs.computeNextTerm(term, resultLhs);
            transformRhs.computeNextTerm(term * c, resultRhs);
            const Real error((resultLhs * c - resultRhs).lInfinityNorm());
            CPPUNIT_ASSERT(error < 1.0e-12);
          }
        }
      }
    }
  }
  void testReverseTransformLinearAdd(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by addition
      for (size_t m = 0; m < NumCoordinateChangeTestCases; ++m) {
        const Polynomial polA(randomPoly(algebra.getNumVars()));
        const Polynomial polB(randomPoly(algebra.getNumVars()));
        const Polynomial polC(polA + polB);
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformA(algebra, wInnerTerms, true);
          CoordinateChange transformB(algebra, wInnerTerms, true);
          CoordinateChange transformC(algebra, wInnerTerms, true);
          Polynomial newA(algebra.zero());
          Polynomial newB(algebra.zero());
          Polynomial newC(algebra.zero());
          for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
            const Polynomial termA(algebra.isoGradePart(polA, row + 1));
            const Polynomial termB(algebra.isoGradePart(polB, row + 1));
            const Polynomial termC(algebra.isoGradePart(polC, row + 1));
            transformA.computeNextTerm(termA, newA);
            transformB.computeNextTerm(termB, newB);
            transformC.computeNextTerm(termC, newC);
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
    CoordinateChange transform(algebra, wInnerTerms, false);
    Polynomial result(algebra.zero());
    for (Power i = 0; i < 4; ++i) {
      transform.computeNextTerm(algebra.zero(), result);
      CPPUNIT_ASSERT(result.isZero());
    }
  }
  void testForwardTransformZeroIsZero(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //transform zero
      CoordinateChange transform(algebra, wInnerTerms, false);
      Polynomial result(algebra.zero());
      for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
        transform.computeNextTerm(algebra.zero(), result);
        CPPUNIT_ASSERT(result.isZero());
      }
    }
  }
  void testForwardTransformIsNearIdentity(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //transform single coordinate one term
      for (Index i = 0; i < algebra.getNumVars(); ++i) {
        CoordinateChange transform(algebra, wInnerTerms, false);
        const Polynomial xI(algebra.coordinateMonomial(i));
        Polynomial result(algebra.zero());
        transform.computeNextTerm(xI, result);
        CPPUNIT_ASSERT(result == xI);
      }
    }
  }
  void testForwardTransformLinearMulOnCoordinates(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by multiplication
      for (Index s = 0; s < algebra.getNumVars(); ++s) {
        const Polynomial xS(algebra.coordinateMonomial(s));
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformLhs(algebra, wInnerTerms, false);
          CoordinateChange transformRhs(algebra, wInnerTerms, false);
          Polynomial termLhs(algebra.zero());
          Polynomial termRhs(algebra.zero());
          for (Power i = 0; i < NumCoordinateChangeTestRows; ++i) {
            if (i == 0) {
              transformLhs.computeNextTerm(xS, termLhs);
              transformRhs.computeNextTerm(xS * c, termRhs);
              CPPUNIT_ASSERT((termRhs - (xS * c)).lInfinityNorm() < 1.0e-12);
            }
            else {
              transformLhs.computeNextTerm(algebra.zero(), termLhs);
              transformRhs.computeNextTerm(algebra.zero(), termRhs);
            }
            CPPUNIT_ASSERT((termRhs - termLhs * c).lInfinityNorm() < 1.0e-12);
          }
        }
      }
    }
  }
  void testForwardTransformLinearMul(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by multiplication
      for (size_t m = 0; m < NumCoordinateChangeTestCases; ++m) {
        const Polynomial f(randomPoly(algebra.getNumVars()));
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformLhs(algebra, wInnerTerms, false);
          CoordinateChange transformRhs(algebra, wInnerTerms, false);
          Polynomial resultLhs(algebra.zero());
          Polynomial resultRhs(algebra.zero());
          for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
            const Polynomial term(algebra.isoGradePart(f, row + 1));
            transformLhs.computeNextTerm(term, resultLhs);
            transformRhs.computeNextTerm(term * c, resultRhs);
            const Real error((resultLhs * c - resultRhs).lInfinityNorm());
            CPPUNIT_ASSERT(error < 1.0e-12);
          }
        }
      }
    }
  }
  void testForwardTransformLinearAdd(void) {
    for (size_t n = 0; n < NumCoordinateChangeTestCases; ++n) {
      const Index dof(1 + randomIndex(5));
      const ClassicalLieAlgebra algebra(dof);
      std::vector< Polynomial > wInnerTerms(randomInnerGenTerms(algebra));
      //linearity by addition
      for (size_t m = 0; m < NumCoordinateChangeTestCases; ++m) {
        const Polynomial polA(randomPoly(algebra.getNumVars()));
        const Polynomial polB(randomPoly(algebra.getNumVars()));
        const Polynomial polC(polA + polB);
        for (size_t k = 0; k < NumCoordinateChangeTestCases; ++k) {
          const Coeff c(randomCoeff(4.0));
          CoordinateChange transformA(algebra, wInnerTerms, false);
          CoordinateChange transformB(algebra, wInnerTerms, false);
          CoordinateChange transformC(algebra, wInnerTerms, false);
          Polynomial newA(algebra.zero());
          Polynomial newB(algebra.zero());
          Polynomial newC(algebra.zero());
          for (Power row = 0; row < NumCoordinateChangeTestRows; ++row) {
            const Polynomial termA(algebra.isoGradePart(polA, row + 1));
            const Polynomial termB(algebra.isoGradePart(polB, row + 1));
            const Polynomial termC(algebra.isoGradePart(polC, row + 1));
            transformA.computeNextTerm(termA, newA);
            transformB.computeNextTerm(termB, newB);
            transformC.computeNextTerm(termC, newC);
            CPPUNIT_ASSERT(((newA + newB) - newC).lInfinityNorm() < 1.0e-12);
          }
        }
      }
    }
  }
 private:
}; //CoordinateChangeTest

#endif //COORDINATE_CHANGE_TEST_H
