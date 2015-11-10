#ifndef GMP_TEST_H
#define GMP_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: GmpTest
//
// PURPOSE:
//
// Tests for GNU Multi-Precision support.
//
// NOTES:
//
//----------------------------------------------------------------------

//system
#include <iostream>

//library
#include <cppunit/extensions/HelperMacros.h>
#include <gmpxx.h>

class GmpBasicTest : public CppUnit::TestFixture {

  //the following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(GmpBasicTest);

  //Make 1.0:
  CPPUNIT_TEST(testMakeOne);
  CPPUNIT_TEST(testMakeOneBits1);
  CPPUNIT_TEST(testMakeOneBits129);
  CPPUNIT_TEST(testMakeOnePlusEpsilon);

  //Default precision:
  CPPUNIT_TEST(testGetDefaultPrecision);

  CPPUNIT_TEST_SUITE_END();

private:

protected:
  // void setUp() { }
  // void tearDown() { }

  bool isOne(const mpf_class& one) {
    if (!(one == 1.0)) return false;
    if (!(one != 0.0)) return false;
    if (!(one != -1.0)) return false;
    if (!(one != 1.00001)) return false;
    return true;
  }
  void testMakeOne() {
    mpf_class one("1.0");
    CPPUNIT_ASSERT(isOne(one));
  }

  void testMakeOneBits1() {
    mpf_class one("1.0", 1);
    CPPUNIT_ASSERT(isOne(one));
  }

  void testMakeOneBits129() {
    mpf_class one("1.0", 129);
    CPPUNIT_ASSERT(isOne(one));
  }

  void testMakeOnePlusEpsilon() {
    const size_t bits = 13357;
    mpf_class one("1.0", bits);
    mpf_class zer("0.0", bits);
    mpf_class eps("1.0e-51", bits);
    mpf_class tol("1.0e-4000");
    CPPUNIT_ASSERT(isOne(one));
    CPPUNIT_ASSERT(zer == 0.0);
    CPPUNIT_ASSERT(zer < 1.0);
    CPPUNIT_ASSERT(zer*one == zer);
    CPPUNIT_ASSERT(eps > 0.0);
    CPPUNIT_ASSERT(eps < 1.0);
    CPPUNIT_ASSERT(eps > zer);
    CPPUNIT_ASSERT(eps < one);
    CPPUNIT_ASSERT(one+eps > one);
    //std::cerr << (one*eps) << " " << eps << std::endl;
    //std::cerr << (one*eps)-eps << std::endl;
    CPPUNIT_ASSERT(fabs((one*eps)-eps) < tol);
    CPPUNIT_ASSERT(one+eps > 1.0);
  }

  void testGetDefaultPrecision(void) {
    const size_t prec = mpf_get_default_prec();
    mpf_set_default_prec(512);
    CPPUNIT_ASSERT(mpf_get_default_prec() == 512);
    mpf_set_default_prec(prec);
    CPPUNIT_ASSERT(mpf_get_default_prec() == prec);
  }

}; //GmpBasicTest

#endif //GMP_TEST_H

