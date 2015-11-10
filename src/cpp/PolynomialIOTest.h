#ifndef POLYNOMIAL_IO_TEST_H
#define POLYNOMIAL_IO_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialIOTest
//
//----------------------------------------------------------------------

//system headers
#include <utility>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "Polynomial.h"
#include "PolynomialIO.h"

class PolynomialIOTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(PolynomialIOTest);
  //CPPUNIT_TEST(testFixtureIncluded);
  //
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testX);
  CPPUNIT_TEST(testXVector);
  //
  CPPUNIT_TEST_SUITE_END();
private:
protected:
  //void testFixtureIncluded(void) {
  //  CPPUNIT_ASSERT(false);
  //}
  void testZero(void) {
    Polynomial p(3);
    std::string oStr;
    std::stringstream str(oStr, std::stringstream::out | std::stringstream::in);
    seWritePolynomial(str, p);
    str.flush();
    const Polynomial q(seReadPolynomial(str));
    CPPUNIT_ASSERT(q == p);
  }
  void testX(void) {
    PolynomialRing ring(5);
    for (Index i = 0; i < ring.getNumVars(); ++i) {
      const Polynomial p(ring.coordinateMonomial(i));
      std::string oStr;
      std::stringstream str(oStr, std::stringstream::out | std::stringstream::in);
      seWritePolynomial(str, p);
      str.flush();
      const Polynomial q(seReadPolynomial(str));
      CPPUNIT_ASSERT(q == p);
    }
  }
  void testXVector(void) {
    PolynomialRing ring(5);
    std::vector< Polynomial > polys;
    for (Index i = 0; i < ring.getNumVars(); ++i) {
      const Polynomial p(ring.coordinateMonomial(i));
      polys.push_back(p);
    }
    std::string oStr;
    std::stringstream str(oStr, std::stringstream::out | std::stringstream::in);
    seWriteVectorOfPolynomials(str, polys);
    str.flush();
    std::vector< Polynomial > result;
    seReadVectorOfPolynomials(str, result);
    CPPUNIT_ASSERT(result == polys);
  }
}; //PolynomialIOTest

#endif //POLYNOMIAL_IO_TEST_H
