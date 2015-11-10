#ifndef XXX_TEST_H
#define XXX_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: XxxTest
//
// PURPOSE:
//
// Tests for Xxx.
//
// NOTES:
//
//----------------------------------------------------------------------

//standard headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers

class XxxTest : public CppUnit::TestFixture {
 public:
  CPPUNIT_TEST_SUITE(XxxTest);
  CPPUNIT_TEST(testFixtureIncluded);
  CPPUNIT_TEST_SUITE_END();
 protected:
  void testFixtureIncluded(void) {
    CPPUNIT_ASSERT(false);
  }
 private:
}; //XxxTest

#endif //XXX_TEST_H
