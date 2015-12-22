#ifndef DEFAULT_NUMERICAL_PRECISION_TEST_H
#define DEFAULT_NUMERICAL_PRECISION_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DefaultNumericalPrecisionTest
//
// PURPOSE:
//
// Tests for DefaultNumericalPrecision.
//
// NOTES:
//
//----------------------------------------------------------------------

//standard headers

//library headers
#include <cppunit/extensions/HelperMacros.h>
#include <gmpxx.h>

//project headers
#include "DefaultNumericalPrecision.h"

class DefaultNumericalPrecisionTest : public CppUnit::TestFixture {
 public:
  CPPUNIT_TEST_SUITE(DefaultNumericalPrecisionTest);
  //CPPUNIT_TEST(testFixtureIncluded);
  CPPUNIT_TEST(testSetPrecision);
  CPPUNIT_TEST(testGetNumBits);
  CPPUNIT_TEST(testGetNumDigits);
  CPPUNIT_TEST(testScoping);
  CPPUNIT_TEST_SUITE_END();
 protected:
  //void testFixtureIncluded(void) {
  //  CPPUNIT_ASSERT(false);
  //}
  void testSetPrecision(void) {
    DefaultNumericalPrecision prec0(64);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 64);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 128);
    DefaultNumericalPrecision prec1(128);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 128);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 256);
    DefaultNumericalPrecision prec2(256);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 256);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 512);
    DefaultNumericalPrecision prec3(128);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 128);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 256);
    DefaultNumericalPrecision prec4(64);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 64);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 128);
  }
  void testGetNumBits(void) {
    DefaultNumericalPrecision prec(64);
    const size_t p = mpf_get_default_prec();
    CPPUNIT_ASSERT(prec.getNumBits() == p);
  }
  void testGetNumDigits(void) {
    DefaultNumericalPrecision prec(53);
    CPPUNIT_ASSERT(prec.getNumDigits() == 19);
  }
  void testScoping(void) {
    DefaultNumericalPrecision prec(64);
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 64);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 128);
    {
      DefaultNumericalPrecision prec1(128);
      CPPUNIT_ASSERT(mpf_get_default_prec() >= 128);
      CPPUNIT_ASSERT(mpf_get_default_prec() < 256);
      {
        DefaultNumericalPrecision prec2(72);
        CPPUNIT_ASSERT(mpf_get_default_prec() >= 72);
        CPPUNIT_ASSERT(mpf_get_default_prec() <= 128);
      }
      CPPUNIT_ASSERT(mpf_get_default_prec() >= 128);
      CPPUNIT_ASSERT(mpf_get_default_prec() < 256);
    }
    CPPUNIT_ASSERT(mpf_get_default_prec() >= 64);
    CPPUNIT_ASSERT(mpf_get_default_prec() < 128);
  }
 private:
}; //DefaultNumericalPrecisionTest

#endif //DEFAULT_NUMERICAL_PRECISION_TEST_H
