#ifndef BASIC_TEST_H
#define BASIC_TEST_H

/// @file BasicTest.h
///
/// @brief Some basic tests to check compilation environment, etc.
///
/// Check that we can compile code in the current environment,
/// and act as a reminder concerning, for example, the action of fabs.
///
/// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

//system
#include <iostream>
#include <utility> //for pair

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "StlInit.h"

///set a default numerical tolerance for testing purposes
static const Real BasicTestFabsTolerance(1.0e-15);

///text fixture class
class BasicTest : public CppUnit::TestFixture {

 public:
  ///create static CppUnit::TestSuite *suite():
  CPPUNIT_TEST_SUITE(BasicTest);
  //fabs:
  CPPUNIT_TEST(testFabs);
  CPPUNIT_TEST(testFabsOnCoeff);
  ///end the test suite:
  CPPUNIT_TEST_SUITE_END();

 protected:
  ///test the action of the fabs function on some non-integer values
  void testFabs(void) {
    CPPUNIT_ASSERT(fabs(1.01) == 1.01);
    CPPUNIT_ASSERT(fabs(0.01) == 0.01);
    CPPUNIT_ASSERT(fabs(-1.01) == 1.01);
    CPPUNIT_ASSERT(fabs(-0.01) == 0.01);
  }
  ///test that fabs has the required behaviour on our complex and real types
  void testFabsOnCoeff(void) {
    CPPUNIT_ASSERT(fabs(Complex(1.01)) - Real(1.01) < BasicTestFabsTolerance);
    CPPUNIT_ASSERT(fabs(Complex(0.01)) - Real(0.01) < BasicTestFabsTolerance);
    CPPUNIT_ASSERT(fabs(Complex(-1.01)) - Real(1.01) < BasicTestFabsTolerance);
    CPPUNIT_ASSERT(fabs(Complex(-0.01)) - Real(0.01) < BasicTestFabsTolerance);
  }

 private:
}; //BasicTest

#endif //BASIC_TEST_H
