//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: TestRunner implementation
//
//----------------------------------------------------------------------

//standard headers
#include <cstdlib>

//library headers
#include <cppunit/ui/text/TestRunner.h>

//project headers
#include "TestRunner.h"
//
#include "TypesTest.h"
#include "UtilityTest.h"
#include "GmpTest.h"
#include "BasicTest.h"
#include "StlInitTest.h"
#include "DefaultNumericalPrecisionTest.h"
//
#include "MapPowersTest.h"
#include "PolynomialTest.h"
#include "PolynomialTestRing.h"
#include "PolynomialTestOthers.h"
#include "PolynomialTestSubstitute.h"
#include "PolynomialRingTest.h"
#include "PolynomialIOTest.h"
//
#include "ClassicalLieAlgebraTest.h"
#include "ClassicalLieAlgebraTestPoisson.h"
#include "SemiClassicalLieAlgebraTest.h"
#include "SemiClassicalLieAlgebraTestMoyal.h"
//
#include "DepritTriangleTest.h"
#include "CoordinateChangeTest.h"
#include "DepritTriangleKnownTest.h"

int main(int argc, char **argv) {

  //Create test runner to aggregate test suites
  CppUnit::TextUi::TestRunner runner;

  //Infrastructure
  runner.addTest(DefaultNumericalPrecisionTest::suite());
  runner.addTest(UtilityTest::suite());
  runner.addTest(TypesTest::suite());
  runner.addTest(GmpBasicTest::suite());
  runner.addTest(BasicTest::suite());
  runner.addTest(StlInitTest::suite());

  //Polynomials
  runner.addTest(MapPowersTest::suite());
  runner.addTest(PolynomialIOTest::suite());
  runner.addTest(PolynomialTest::suite());
  runner.addTest(PolynomialTestOthers::suite());
  runner.addTest(PolynomialTestRing::suite());
  runner.addTest(PolynomialTestSubstitute::suite());
  runner.addTest(PolynomialRingTest::suite());

  //Deprit's method
  runner.addTest(CoordinateChangeTest::suite());
  runner.addTest(DepritTriangleTest::suite());
  runner.addTest(DepritTriangleKnownTest::suite());

  //Classical and semiclassical Lie algebras
  runner.addTest(ClassicalLieAlgebraTest::suite());
  runner.addTest(ClassicalLieAlgebraTestPoisson::suite());
  runner.addTest(SemiClassicalLieAlgebraTest::suite());
  runner.addTest(SemiClassicalLieAlgebraTestMoyal::suite());

  //Run all tests
  runner.run();

  return EXIT_SUCCESS;
}
