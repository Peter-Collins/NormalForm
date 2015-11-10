//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialRingTestRunner implementation
//
//----------------------------------------------------------------------

//standard headers
#include <cstdlib>

//library headers
#include <cppunit/ui/text/TestRunner.h>

//project headers
#include "PolynomialRingTestRunner.h"
#include "PolynomialRingTest.h"

int main( int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( PolynomialRingTest::suite() );
  runner.run();

  return EXIT_SUCCESS;
}
