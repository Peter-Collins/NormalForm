//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef STL_INIT_TEST_H
#define STL_INIT_TEST_H

//system
#include <assert.h>
#include <vector>
#include <iostream>

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "StlInit.h"

class StlInitTest : public CppUnit::TestFixture {

  //the following create static CppUnit::TestSuite *suite()

  CPPUNIT_TEST_SUITE(StlInitTest);
  CPPUNIT_TEST(test_do_nothing);
  CPPUNIT_TEST(test_assign);
  CPPUNIT_TEST(test_assign_one);
  CPPUNIT_TEST(test_assign_eq);
  CPPUNIT_TEST(test_comma);
  CPPUNIT_TEST_SUITE_END();

private:

protected:

  // void setUp(void) { }

  // void tearDown(void) { }

  void test_do_nothing(void) {
    std::vector< float > v;
    init(v);
    CPPUNIT_ASSERT(v.size() == 0);
  }

  void test_assign(void) {
    std::vector< float > v;
    init(v) = 3.0;
    CPPUNIT_ASSERT(v.size() == 1);
    CPPUNIT_ASSERT(v[0] == float(3.0));
  }

  void test_assign_one(void) {
    std::vector< float > v;
    init(v) += 3.0;
    CPPUNIT_ASSERT(v.size() == 1);
    CPPUNIT_ASSERT(v[0] == float(3.0));
  }

  void test_comma(void) {
    std::vector< float > v;
    init(v) = 1.0, -3.4, 5.6, 77, 42.5;
    CPPUNIT_ASSERT(v[0] == float(1.0));
    CPPUNIT_ASSERT(v[1] == float(-3.4));
    CPPUNIT_ASSERT(v[2] == float(5.6));
    CPPUNIT_ASSERT(v[3] == float(77));
    CPPUNIT_ASSERT(v[4] == float(42.5));
    CPPUNIT_ASSERT(v.size() == 5);
  }

  void test_assign_eq(void) {
    std::vector< float > v;
    init(v) = 1.0, -3.4;
    init(v) += 5.6;
    init(v) += 77, 42.5;
    CPPUNIT_ASSERT(v[0] == float(1.0));
    CPPUNIT_ASSERT(v[1] == float(-3.4));
    CPPUNIT_ASSERT(v[2] == float(5.6));
    CPPUNIT_ASSERT(v[3] == float(77));
    CPPUNIT_ASSERT(v[4] == float(42.5));
    CPPUNIT_ASSERT(v.size() == 5);
  }

}; //StlInitTest

#endif //STL_INIT_TEST_H
