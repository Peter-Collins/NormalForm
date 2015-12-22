#ifndef UTILITY_TEST_H
#define UTILITY_TEST_H

/// @file UtilityTest.h
///
/// @brief Unit tests for basic utility methods, e.g., binomial coeffs.
///
/// Test that simple functions such as binomial coefficients and
/// factorials are implemented without overflow for suitable ranges,
/// and that they work across the various integer types.
///
/// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

//system headers

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "Utility.h"

///test fixture class
class UtilityTest : public CppUnit::TestFixture {
 public:

  ///create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(UtilityTest);

  //factorial
  CPPUNIT_TEST_EXCEPTION(testFactorialNegative, UtilityNegativeArgumentError);
  CPPUNIT_TEST(testFactorialZero);
  CPPUNIT_TEST(testFactorialOne);
  CPPUNIT_TEST(testFactorialExamplesAccumulate);
  CPPUNIT_TEST(testFactorialExamplesRelation);

  //binomial
  CPPUNIT_TEST_EXCEPTION(testBinomialNegative0, UtilityNegativeArgumentError);
  CPPUNIT_TEST_EXCEPTION(testBinomialNegative1, UtilityNegativeArgumentError);
  CPPUNIT_TEST_EXCEPTION(testBinomialBounds, UtilityOutOfRangeError);
  CPPUNIT_TEST(testBinomialZero);
  CPPUNIT_TEST(testBinomialAnyZero);
  CPPUNIT_TEST(testBinomialAnyOne);
  CPPUNIT_TEST(testBinomialAnyNMinusOne);
  CPPUNIT_TEST(testBinomialSymmetry);
  CPPUNIT_TEST(testBinomialRelations);

  //intPow
  CPPUNIT_TEST(testIntPowZero);
  CPPUNIT_TEST(testIntPowOne);
  CPPUNIT_TEST(testIntPowRelation);
  CPPUNIT_TEST(testIntPowNegative);

  //coeffPow
  //CPPUNIT_TEST(testCoeffPowZero);
  //CPPUNIT_TEST(testCoeffPowOne);
  //CPPUNIT_TEST(testCoeffPowExamples);

  ///end the test suite
  CPPUNIT_TEST_SUITE_END();

 protected:
  ///factorial(x) should throw an exception for all negative x
  void testFactorialNegative(void) {
    factorial(-1);
  }
  ///factorial(zero) == 1 for all suitable types of zero and one
  void testFactorialZero(void) {
    CPPUNIT_ASSERT(Index(1) == factorial(Index(0)));
    CPPUNIT_ASSERT(Power(1) == factorial(Power(0)));
    CPPUNIT_ASSERT(size_t(1) == factorial(size_t(0)));
  }
  ///factorial(one) == one for all suitable types of one
  void testFactorialOne(void) {
    CPPUNIT_ASSERT(Index(1) == factorial(Index(1)));
    CPPUNIT_ASSERT(Power(1) == factorial(Power(1)));
    CPPUNIT_ASSERT(size_t(1) == factorial(size_t(1)));
  }
  ///test factorial values vs accumulation in a loop
  void testFactorialExamplesAccumulate(void) {
    Index indexR = 1;
    for (Index i = 1; i < 10; ++i ) {
      indexR *= i;
      CPPUNIT_ASSERT(indexR == factorial(i));
    }
    Power powerR = 1;
    for (Power i = 1; i < 10; ++i) {
      powerR *= i;
      CPPUNIT_ASSERT(powerR == factorial(i));
    }
    size_t sizeR = 1;
    for (size_t i = 1; i < 10; ++i) {
      sizeR *= i;
      CPPUNIT_ASSERT(sizeR == factorial(i));
    }
  }
  ///test the recursion relation n! == (n-1)!n
  void testFactorialExamplesRelation(void) {
    Index indexR = 1;
    for (Index i = 1; i < 10; ++i ) {
      indexR *= i;
      CPPUNIT_ASSERT(factorial(i) == i*factorial(i-1));
    }
    Power powerR = 1;
    for (Power i = 1; i < 10; ++i) {
      powerR *= i;
      CPPUNIT_ASSERT(factorial(i) == i*factorial(i-1));
    }
    size_t sizeR = 1;
    for (size_t i = 1; i < 10; ++i) {
      sizeR *= i;
      CPPUNIT_ASSERT(factorial(i) == i*factorial(i-1));
    }
  }
  ///binomial(x, y) should throw exception for all negative x
  void testBinomialNegative0(void) {
    binomial(-1, 0);
  }
  ///binomial(x, y) should throw exception for all negative y
  void testBinomialNegative1(void) {
    binomial(1, -1);
  }
  ///binomial(x, y) should throw exception if x < y
  void testBinomialBounds(void) {
    binomial(3, 4);
  }
  ///binomial(0, 0) == 1
  void testBinomialZero(void) {
    CPPUNIT_ASSERT(binomial(0, 0) == 1);
  }
  ///binomial(x, 0) == 1 for all nonnegative x
  void testBinomialAnyZero(void) {
    for (size_t i = 0; i < 100; ++i)
      CPPUNIT_ASSERT(binomial(i, 0) == 1);
  }
  ///binomial(x, 1) == x for all nonnegative x
  void testBinomialAnyOne(void) {
    for (size_t i = 1; i < 100; ++i)
      CPPUNIT_ASSERT(binomial(i, 1) == i);
  }
  ///binomial(x, x-1) == x for all nonnegative x
  void testBinomialAnyNMinusOne(void) {
    for (size_t i = 1; i < 100; ++i)
      CPPUNIT_ASSERT(binomial(i, i-1) == i);
  }
  ///test for symmetry, binomial(i,j) == binomial(i,i-j)
  void testBinomialSymmetry(void) {
    for (size_t i = 0; i < 10; ++i)
      for (size_t j = 0; j <= i; ++j)
        CPPUNIT_ASSERT(binomial(i, j) == binomial(i, i-j));
  }
  ///test recursion relation
  void testBinomialRelations(void) {
    for (size_t n = 1; n < 10; ++n)
      for (size_t k = 0; k <= n; ++k)
        if ((k == 0) || (k==n))
          CPPUNIT_ASSERT(binomial(n, k) == 1);
        else
          CPPUNIT_ASSERT(binomial(n, k) == binomial(n-1, k-1)+binomial(n-1, k));
  }
  ///integer raised to power zero should give one
  void testIntPowZero(void) {
    for (size_t n = 0; n < 10; ++n)
      CPPUNIT_ASSERT(intPow(n, 0) == 1);
  }
  ///integer raised to power one should return the integer
  void testIntPowOne(void) {
    for (size_t n = 0; n < 10; ++n)
      CPPUNIT_ASSERT(intPow(n, 1) == n);
  }
  ///test recursion relation for integer powers for a range of values
  ///9**10 is too big for 32 bit int
  void testIntPowRelation(void) {
    for (size_t n = 0; n < 10; ++n)
      for (size_t p = 0; p < 9; ++p) {
	///if (intPow(n, p+1) != n*intPow(n, p))
	///  std::cerr << n << p << intPow(n, p+1) << std::endl;
        CPPUNIT_ASSERT(intPow(n, p+1) == n*intPow(n, p));
      }
  }
  ///test basic properties of negative integers raised to powers
  void testIntPowNegative(void) {
    for (int n = 0; n > -10; --n)
      for (size_t p = 0; p < 10; ++p)
        {
          CPPUNIT_ASSERT(intPow(n, p+1) == n*intPow(n, p));
          if (p % 2)
            if (n == 0)
              CPPUNIT_ASSERT(intPow(n, p) == 0);
            else
              CPPUNIT_ASSERT(intPow(n, p) < 0);
          else
            CPPUNIT_ASSERT(intPow(n, p) >= 0);
        }
  }

 private:

}; //UtilityTest

#endif //UTILITY_TEST_H
