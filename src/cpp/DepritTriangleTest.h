//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangleTest
//
// PURPOSE:
//
// Tests for DepritTriangle.
//
// NOTES:
//
// In particular the tests for invariance of the generator must a
// handled carefully: to transform a generator term via a Deprit
// triangle, we require the generator term of grade one higher.  This
// means that one cannot check the highest grade term for invariance?
//
//----------------------------------------------------------------------

#ifndef DEPRIT_TRIANGLE_TEST_H
#define DEPRIT_TRIANGLE_TEST_H

//standard headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "DefaultNumericalPrecision.h"
#include "DepritTriangle.h"

const size_t NumDepritTriangleTestCases(16);

class DepritTriangleTest : public CppUnit::TestFixture {
 public:
  CPPUNIT_TEST_SUITE(DepritTriangleTest);
  //
  //CPPUNIT_TEST(testFixtureIncluded);
  //
  CPPUNIT_TEST_SUITE_END();
 protected:
  //void testFixtureIncluded(void) {
  //  CPPUNIT_ASSERT(false);
  //}
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
}; //DepritTriangleTest

#endif //DEPRIT_TRIANGLE_TEST_H
