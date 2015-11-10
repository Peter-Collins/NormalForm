#ifndef POLYNOMIAL_RING_TEST_H
#define POLYNOMIAL_RING_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
//----------------------------------------------------------------------

//system headers
#include <utility>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "PolynomialRing.h"

static const Index PolyRingTestMaxNumVars = 100;
static const Power PolyRingTestMaxPower = 8;

std::vector< Coeff > makeVector(const Index numVars,
                                const Coeff val = CoeffZero);

std::vector< Coeff > makeVector(const Index numVars,
                                const Coeff val) {
  std::vector< Coeff > vec;
  for (Index i = 0; i < numVars; ++i)
    vec.push_back(val);
  return vec;
}

class PolynomialRingTest : public CppUnit::TestFixture {
 public:
  //the following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(PolynomialRingTest);
  //constructor from number of variables
  CPPUNIT_TEST_EXCEPTION(testConstructorGivenZeroVars, PolynomialRingIndexError);
  CPPUNIT_TEST(testConstructorGivenNumVars);
  //big four
  CPPUNIT_TEST(testDefaultConstructor);
  CPPUNIT_TEST(testCopyConstructor);
  //CPPUNIT_TEST(testAssign); //forbidden
  //basic object creation
  CPPUNIT_TEST(testZeroPolynomial);
  CPPUNIT_TEST(testOnePolynomial);
  CPPUNIT_TEST(testCoordinateMonomial);
  CPPUNIT_TEST(testCoordinateMonomialPower);
  //CPPUNIT_TEST(testMonomialPolynomial);
  CPPUNIT_TEST_SUITE_END();
protected:
  void testConstructorGivenZeroVars(void) {
    //this should throw an index error
    PolynomialRing ring(0);
  }
  void testConstructorGivenNumVars(void) {
    //construct rings with given number of variables
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      CPPUNIT_ASSERT(ring.getNumVars() == numVars);
      Polynomial p(numVars);
      CPPUNIT_ASSERT(ring.hasElt(p));
    }
  }
  void testCopyConstructor(void) {
    //construct a copy of a ring
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      PolynomialRing ringCopy(ring);
      CPPUNIT_ASSERT(ring.getNumVars() == numVars);
      CPPUNIT_ASSERT(ringCopy.getNumVars() == numVars);
    }
  }
  void testDefaultConstructor(void) {
    //test the default constructor gives a 1-variable ring
    PolynomialRing ring;
    CPPUNIT_ASSERT(ring.getNumVars() == 1);
    Polynomial p;
    CPPUNIT_ASSERT(ring.hasElt(p));
  }
  void testZeroPolynomial(void) {
    //properties of zero polynomial and compatibility with poly zero
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      Polynomial ringZero = ring.zero();
      Polynomial polyZero = Polynomial::Zero(numVars);
      CPPUNIT_ASSERT(ring.hasElt(ringZero));
      CPPUNIT_ASSERT(ring.hasElt(polyZero));
      CPPUNIT_ASSERT(ringZero.getNumVars() == ring.getNumVars());
      CPPUNIT_ASSERT(ringZero.isZero());
      CPPUNIT_ASSERT(ringZero.isConstant());
      CPPUNIT_ASSERT(ringZero.isHomogeneous());
      CPPUNIT_ASSERT(ringZero.isHomogeneous(0));
      CPPUNIT_ASSERT(ringZero == polyZero);
      CPPUNIT_ASSERT(ringZero(makeVector(numVars, CoeffZero)) == CoeffZero);
      CPPUNIT_ASSERT(ringZero(makeVector(numVars, CoeffOne)) == CoeffZero);
      CPPUNIT_ASSERT(ringZero(makeVector(numVars, CoeffJ)) == CoeffZero);
    }
  }
  void testOnePolynomial(void) {
    //properties of one polynomial and compatibility with poly one
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      Polynomial ringOne = ring.one();
      Polynomial polyOne = Polynomial::One(numVars);
      CPPUNIT_ASSERT(ring.hasElt(ringOne));
      CPPUNIT_ASSERT(ring.hasElt(polyOne));
      CPPUNIT_ASSERT(ringOne.getNumVars() == ring.getNumVars());
      CPPUNIT_ASSERT(!ringOne.isZero());
      CPPUNIT_ASSERT(ringOne.isConstant());
      CPPUNIT_ASSERT(ringOne.isHomogeneous());
      CPPUNIT_ASSERT(ringOne.isHomogeneous(0));
      CPPUNIT_ASSERT(ringOne == polyOne);
      CPPUNIT_ASSERT(ringOne(makeVector(numVars, CoeffZero)) == CoeffOne);
      CPPUNIT_ASSERT(ringOne(makeVector(numVars, CoeffOne)) == CoeffOne);
      CPPUNIT_ASSERT(ringOne(makeVector(numVars, CoeffJ)) == CoeffOne);
    }
  }
  void testCoordinateMonomial(void) {
    //test the coordinate monomials x_{i}
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      std::vector< Coeff > vec;
      for (Index j = 0; j < numVars; ++j)
        vec.push_back(Coeff(j));
      for (Index i = 0; i < numVars; ++i) {
        Polynomial xI = ring.coordinateMonomial(i);
        CPPUNIT_ASSERT(xI.getNumVars() == ring.getNumVars());
        CPPUNIT_ASSERT(!xI.isZero());
        CPPUNIT_ASSERT(!xI.isConstant());
        CPPUNIT_ASSERT(xI.isHomogeneous());
        CPPUNIT_ASSERT(xI.isHomogeneous(1));
        CPPUNIT_ASSERT(xI(vec) == Coeff(i));
      }
    }
  }
  void testCoordinateMonomialPower(void) {
    //test powers of coordinate monomials x_{i}^{p}
    for (Index numVars = 1; numVars < PolyRingTestMaxNumVars; ++numVars) {
      PolynomialRing ring(numVars);
      std::vector< Coeff > vec;
      for (Index j = 0; j < numVars; ++j)
        vec.push_back(Coeff(j));
      for (Index i = 0; i < numVars; ++i) {
        for (Power p = 0; p < PolyRingTestMaxPower; ++p) {
          Polynomial xI = ring.coordinateMonomial(i, p);
          CPPUNIT_ASSERT(xI.getNumVars() == ring.getNumVars());
          if (p == 0) {
            //x^0 gives one
            CPPUNIT_ASSERT(!xI.isZero());
            CPPUNIT_ASSERT(xI.isConstant());
          }
          else {
            CPPUNIT_ASSERT(!xI.isZero());
            CPPUNIT_ASSERT(!xI.isConstant());
          }
          CPPUNIT_ASSERT(xI.isHomogeneous());
          CPPUNIT_ASSERT(xI.isHomogeneous(p));
          CPPUNIT_ASSERT(xI(vec) == coeffPow(Coeff(i), p));
        }
      }
    }
  }
 private:
}; //PolynomialRingTest

#endif //POLYNOMIAL_RING_TEST_H
