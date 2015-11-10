#ifndef TYPES_TEST_H
#define TYPES_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
//----------------------------------------------------------------------

//system headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"

class TypesTest : public CppUnit::TestFixture {

  //the following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(TypesTest);

  CPPUNIT_TEST(testRandomSizeZero);
  CPPUNIT_TEST(testRandomSize);

  CPPUNIT_TEST(testRandomIndexZero);
  CPPUNIT_TEST(testRandomIndex);

  CPPUNIT_TEST(testRandomPowerZero);
  CPPUNIT_TEST(testRandomPower);

  CPPUNIT_TEST(testFabsReal);
  CPPUNIT_TEST(testFabsComplex);

  CPPUNIT_TEST_SUITE_END();

protected:

  void testRandomSizeZero(void) {
    for (size_t iCase = 0; iCase < 100; ++iCase)
      CPPUNIT_ASSERT(randomSize(0) == 0);
  }
  void testRandomSize(void) {
    for (size_t maxSize = 1; maxSize < 100; ++maxSize) {
      size_t p = randomSize(maxSize);
      CPPUNIT_ASSERT(0 <= p);
      CPPUNIT_ASSERT(p < maxSize);
    }
  }
  void testRandomIndexZero(void) {
    for (size_t iCase = 0; iCase < 100; ++iCase)
      CPPUNIT_ASSERT(randomIndex(0) == 0);
  }
  void testRandomIndex(void) {
    for (Index maxIndex = 1; maxIndex < 100; ++maxIndex) {
      Index p = randomIndex(maxIndex);
      CPPUNIT_ASSERT(0 <= p);
      CPPUNIT_ASSERT(p < maxIndex);
    }
  }
  void testRandomPowerZero(void) {
    for (size_t iCase = 0; iCase < 100; ++iCase)
      CPPUNIT_ASSERT(randomPower(0) == 0);
  }
  void testRandomPower(void) {
    for (Power maxPower = 1; maxPower < 100; ++maxPower) {
      Power p = randomPower(maxPower);
      CPPUNIT_ASSERT(0 <= p);
      CPPUNIT_ASSERT(p < maxPower);
    }
  }
  void testFabsReal(void) {
    for (size_t iCase = 0; iCase < 100; ++iCase) {
      const Real r(randomReal(10.0));
      if (r < 0)
        CPPUNIT_ASSERT(fabs(r) == (-r));
      else
        CPPUNIT_ASSERT(fabs(r) == r);
    }
  }
  void testFabsComplex(void) {
    CPPUNIT_ASSERT(fabs(Complex(Real(0), Real(0))) == Real(0));
    //
    CPPUNIT_ASSERT(fabs(Complex(Real(-3), Real(-4))) == Real(5));
    CPPUNIT_ASSERT(fabs(Complex(Real(-3), Real(4))) == Real(5));
    CPPUNIT_ASSERT(fabs(Complex(Real(3), Real(-4))) == Real(5));
    CPPUNIT_ASSERT(fabs(Complex(Real(3), Real(4))) == Real(5));
    //
    for (size_t i = 0; i < 100; ++i ) {
      const Real r(randomReal(1000.0));
      CPPUNIT_ASSERT(fabs(Complex(r, RealZero)) == fabs(r));
      CPPUNIT_ASSERT(fabs(Complex(RealZero, r)) == fabs(r));
    }
  }

private:

}; //TypesTest

#endif //TYPES_TEST_H
