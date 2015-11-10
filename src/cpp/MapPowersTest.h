//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef MAP_POWERS_TEST_H
#define MAP_POWERS_TEST_H

//system headers
#include <assert.h>
#include <cstdlib>
#include <vector>
#include <utility>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "Random.h"
#include "MapPowers.h"

//constants
static const Index MPTLengths = 50;
static const size_t MPTNumCases = 20;
static const Power MPTMaxDegree = 16;

class MapPowersTest : public CppUnit::TestFixture {

  //the following create static CppUnit::TestSuite *suite()

  CPPUNIT_TEST_SUITE(MapPowersTest);

  //Basic construction:
  CPPUNIT_TEST(testDefaultNumVars);
  CPPUNIT_TEST(testConstructByNumVarsNumVars);
  CPPUNIT_TEST(testConstructByNumVarsEltsNumVars);
  CPPUNIT_TEST(testDefaultDegree);

  //Copy construction:
  CPPUNIT_TEST(testDefaultCopyDegree);

  //Assignment:
  CPPUNIT_TEST(testDefaultAssignDegree);
  CPPUNIT_TEST(testDefaultAssignSelfDegree);

  //Construction from length and mapping:
  CPPUNIT_TEST(testConstructFromLenMap1);
  CPPUNIT_TEST(testConstructFromLenMap2);
  CPPUNIT_TEST(testConstructFromLenMap3);
  CPPUNIT_TEST(testConstructFromLenMap4);
  CPPUNIT_TEST(testConstructFromLenMap5);

  //Multiplication:
  CPPUNIT_TEST(testDefaultMulAssignSelfDegree);
  CPPUNIT_TEST(testDefaultMulAssignDegree);
  CPPUNIT_TEST(testMulSelf);
  CPPUNIT_TEST(testSquare);
  CPPUNIT_TEST(testMul);

  //Ordering:
  CPPUNIT_TEST(testDefaultNotLessSelf);
  CPPUNIT_TEST(testTotalOrderingExampleLess);
  CPPUNIT_TEST(testTotalOrderingExampleMore);
  CPPUNIT_TEST(testTotalOrderingExampleLessEq);
  CPPUNIT_TEST(testTotalOrderingExampleMoreEq);
  CPPUNIT_TEST(testDefaultGreaterEqualSelf);
  CPPUNIT_TEST(testExampleGreaterEqualSelf);
  CPPUNIT_TEST(testNotGreaterSelf);
  
  CPPUNIT_TEST_EXCEPTION(testIncompatibleEq, MapPowersSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testIncompatibleNe, MapPowersSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testIncompatibleLe, MapPowersSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testIncompatibleLt, MapPowersSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testIncompatibleGe, MapPowersSizeMismatchError);
  CPPUNIT_TEST_EXCEPTION(testIncompatibleGt, MapPowersSizeMismatchError);

  //Element set:
  CPPUNIT_TEST(testAssignSomePowers);
  CPPUNIT_TEST(testSetZero);
  CPPUNIT_TEST(testSetPositive);
#ifdef SIGNED_POWER_TYPE
  CPPUNIT_TEST_EXCEPTION(testSetNegative, MapPowersValueError);
#endif //SIGNED_POWER_TYPE

  //Differentiation:
  CPPUNIT_TEST(testDefaultDiffCoeff);
  CPPUNIT_TEST(testExampleDiffCoeff);
  CPPUNIT_TEST(testDiff);
#ifdef SIGNED_POWER_TYPE
  CPPUNIT_TEST_EXCEPTION(testDiffNegativeIndex1, MapPowersIndexError);
  CPPUNIT_TEST_EXCEPTION(testDiffNegativeIndex2, MapPowersIndexError);
#endif //SIGNED_POWER_TYPE
  CPPUNIT_TEST_EXCEPTION(testDiffBadIndex1, MapPowersIndexError);
  CPPUNIT_TEST_EXCEPTION(testDiffBadIndex2, MapPowersIndexError);

  //Pow:
  CPPUNIT_TEST(testOnePows);
  CPPUNIT_TEST(testOnePowsAccumulated);
  CPPUNIT_TEST(testCoordPows);
  CPPUNIT_TEST(testCoordPowsAccumulated);
  CPPUNIT_TEST(testExamplePowZeroIsOne);
  CPPUNIT_TEST(testExamplePow);
#ifdef SIGNED_POWER_TYPE
  CPPUNIT_TEST_EXCEPTION(testNegativeIndexPow, MapPowersIndexError);
#endif //SIGNED_POWER_TYPE

  //Call:
  CPPUNIT_TEST_EXCEPTION(testDefaultCallWrongLen, MapPowersSizeMismatchError);
  CPPUNIT_TEST(testDefaultCallAllZero);
  CPPUNIT_TEST_EXCEPTION(testExampleCallWrongLen, MapPowersSizeMismatchError);
  CPPUNIT_TEST(testExampleCallAllZero);
  CPPUNIT_TEST(testExampleCallValues);
  CPPUNIT_TEST(testExampleCall_anyZero);
  CPPUNIT_TEST(testExampleCallNegativeCounts);

  CPPUNIT_TEST_SUITE_END();

private:

protected:

  // void setUp(void) { }

  // void tearDown(void) { }

  //    """Test the properties of the default constructor."""

  void testDefaultDegree(void) {
    MapPowers m;
    CPPUNIT_ASSERT(m.degree() == 0);
    CPPUNIT_ASSERT(m == m);
  }

  void testDefaultCopyDegree(void) {
    MapPowers m;
    MapPowers m2(m);
    CPPUNIT_ASSERT(m2.degree() == 0);
    CPPUNIT_ASSERT(m2 == m);
  }

  void testDefaultAssignDegree(void) {
    MapPowers m;
    MapPowers m2;
    m2 = m;
    CPPUNIT_ASSERT(m2.degree() == 0);
    CPPUNIT_ASSERT(m2 == m);
  }

  void testDefaultAssignSelfDegree(void) {
    MapPowers m;
    m = m;
    CPPUNIT_ASSERT(m.degree() == 0);
  }

  void testDefaultMulAssignSelfDegree(void) {
    MapPowers m;
    m *= m;
    CPPUNIT_ASSERT(m.degree() == 0);
  }

  void testDefaultMulAssignDegree(void) {
    MapPowers m;
    MapPowers m2;
    m *= m2;
    CPPUNIT_ASSERT(m.degree() == 0);
  }

  void testConstructFromLenMap1(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    CPPUNIT_ASSERT(p == p);
  }

  void testConstructFromLenMap2(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    MapPowers p2(p);
    CPPUNIT_ASSERT(p2.degree() == 6);
    CPPUNIT_ASSERT(p2 == p);
  }

  void testConstructFromLenMap3(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    MapPowers p3 = p;
    CPPUNIT_ASSERT(p3.degree() == 6);
    CPPUNIT_ASSERT(p3 == p);
  }

  void testConstructFromLenMap4(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    MapPowers p4(6);
    CPPUNIT_ASSERT(!(p4 == p));
    CPPUNIT_ASSERT(p4 != p);
    p4 = p;
    CPPUNIT_ASSERT(p4.degree() == 6);
    CPPUNIT_ASSERT(p4 == p);
  }

  void testConstructFromLenMap5(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    MapPowers p5(len, m);
    CPPUNIT_ASSERT(p5.degree() == 6);
    CPPUNIT_ASSERT(p5 == p);
  }

  void testAssignSomePowers(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    CPPUNIT_ASSERT(p.degree() == 6);
    CPPUNIT_ASSERT(p[2] == 0);
    CPPUNIT_ASSERT(p[3] == 0);
    CPPUNIT_ASSERT(p[5] == 0);
    MapPowers p2(len);
    p2.setPower(0, 1);
    p2.setPower(1, 0);
    p2.setPower(4, 5);
    CPPUNIT_ASSERT(p2[0] == p[0]);
    CPPUNIT_ASSERT(p2[1] == p[1]);
    CPPUNIT_ASSERT(p2[4] == p[4]);
    CPPUNIT_ASSERT(p2 == p);
    p.setPower(4, p[4]-1);
    CPPUNIT_ASSERT(p2 != p);
    p.setPower(4, p[4]+1);
    CPPUNIT_ASSERT(p2 == p);
  }

  void testMulSelf(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    p *= p;
    CPPUNIT_ASSERT(p[0] == 2);
    CPPUNIT_ASSERT(p[1] == 0);
    CPPUNIT_ASSERT(p[2] == 0);
    CPPUNIT_ASSERT(p[3] == 0);
    CPPUNIT_ASSERT(p[4] == 10);
    CPPUNIT_ASSERT(p[5] == 0);
    CPPUNIT_ASSERT(p.degree() == 12);
  }

  void testSquare(void) {
    const int len = 6;
    IndexToPowerMap m;
    m[0] = 1;
    m[1] = 0;
    m[4] = 5;
    MapPowers p(len, m);
    MapPowers q = p.pow(2);
    CPPUNIT_ASSERT(q[0] == 2);
    CPPUNIT_ASSERT(q[1] == 0);
    CPPUNIT_ASSERT(q[2] == 0);
    CPPUNIT_ASSERT(q[3] == 0);
    CPPUNIT_ASSERT(q[4] == 10);
    CPPUNIT_ASSERT(q[5] == 0);
    CPPUNIT_ASSERT(q.degree() == 12);
  }

  void testDefaultDiffCoeff(void) {
    MapPowers p;
    std::pair<Coeff, MapPowers> c_d = p.diff(0);
    CPPUNIT_ASSERT(fabs(c_d.first) == Real(0.0));
    CPPUNIT_ASSERT(c_d.second == p);
  }

  void testDefaultNumVars(void) {
    MapPowers p;
    CPPUNIT_ASSERT(p.getNumVars() == 1);
  }

  void testConstructByNumVarsNumVars(void) {
    MapPowers p(900);
    CPPUNIT_ASSERT(p.getNumVars() == 900);
  }

  void testConstructByNumVarsEltsNumVars(void) {
    IndexToPowerMap m;
    m[45] = 9;
    MapPowers p(900, m);
    CPPUNIT_ASSERT(p.getNumVars() == 900);
  }

  void testExampleDiffCoeff(void) {
    MapPowers p(6);
    p.setPower(0, 1);
    p.setPower(1, 0);
    p.setPower(4, 5);
    MapPowers p_orig = p;
    MapPowers d0(6);
    d0.setPower(4, 5);
    MapPowers d3(6);
    MapPowers d4(6);
    d4.setPower(0, 1);
    d4.setPower(4, 4);
    std::pair<Coeff, MapPowers> c_d = p.diff(0);
    CPPUNIT_ASSERT(c_d.first == CoeffOne);
    CPPUNIT_ASSERT(c_d.second == d0);
    c_d = p.diff(3);
    CPPUNIT_ASSERT(c_d.first == CoeffZero);
    CPPUNIT_ASSERT(c_d.second == d3);
    c_d = p.diff(4);
    CPPUNIT_ASSERT(c_d.first == Coeff(5.0));
    CPPUNIT_ASSERT(c_d.second == d4);
    CPPUNIT_ASSERT(p == p_orig);
  }

  void testOnePows(void) {
    MapPowers p(1);
    for (int i=0; i<10; ++i) {
      CPPUNIT_ASSERT(p.pow(i) == p);
    }
  }

  void testOnePowsAccumulated(void) {
    MapPowers p(1);
    MapPowers q(p);
    for (int i=0; i<10; ++i) {
      CPPUNIT_ASSERT(p.pow(i) == q);
      q *= p;
    }
  }

  void testCoordPows(void) {
    MapPowers p(1);
    p.setPower(0, 1);
    MapPowers q(1);
    for (int i=0; i<10; ++i) {
      q.setPower(0, i);
      CPPUNIT_ASSERT(p.pow(i) == q);
    }
  }

  void testCoordPowsAccumulated(void) {
    MapPowers p(1);
    p.setPower(0, 1);
    MapPowers q(p);
    for (int i=1; i<10; ++i) {
      CPPUNIT_ASSERT(p.pow(i) == q);
      q *= p;
    }
  }

  void testDefaultNotLessSelf(void) {
    MapPowers p;
    CPPUNIT_ASSERT(!(p < p));
    CPPUNIT_ASSERT(!(p > p));
  }

  void testTotalOrderingExampleLess(void) {
    IndexToPowerMap m;

    MapPowers p0(6, m);
    m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0; m[4] = 0; m[5] = 1;
    MapPowers p1(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 0; m[4] = 0; m[5] = 2;
    MapPowers p2(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 0; m[5] = 0;
    MapPowers p3(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 1; m[5] = 0;
    MapPowers p4(6, m);

    CPPUNIT_ASSERT(!(p0 < p0));
    CPPUNIT_ASSERT(p1 < p0);
    CPPUNIT_ASSERT(p2 < p0);
    CPPUNIT_ASSERT(p3 < p0);
    CPPUNIT_ASSERT(p4 < p0);

    CPPUNIT_ASSERT(!(p0 < p1));
    CPPUNIT_ASSERT(!(p1 < p1));
    CPPUNIT_ASSERT( (p2 < p1));
    CPPUNIT_ASSERT( (p3 < p1));
    CPPUNIT_ASSERT( (p4 < p1));

    CPPUNIT_ASSERT(!(p0 < p2));
    CPPUNIT_ASSERT(!(p1 < p2));
    CPPUNIT_ASSERT(!(p2 < p2));
    CPPUNIT_ASSERT( (p3 < p2));
    CPPUNIT_ASSERT( (p4 < p2));

    CPPUNIT_ASSERT(!(p0 < p3));
    CPPUNIT_ASSERT(!(p1 < p3));
    CPPUNIT_ASSERT(!(p2 < p3));
    CPPUNIT_ASSERT(!(p3 < p3));
    CPPUNIT_ASSERT( (p4 < p3));

    CPPUNIT_ASSERT(!(p0 < p4));
    CPPUNIT_ASSERT(!(p1 < p4));
    CPPUNIT_ASSERT(!(p2 < p4));
    CPPUNIT_ASSERT(!(p3 < p4));
    CPPUNIT_ASSERT(!(p4 < p4));
  }

  void testTotalOrderingExampleMore(void) {
    IndexToPowerMap m;

    MapPowers p0(6, m);
    m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0; m[4] = 0; m[5] = 1;
    MapPowers p1(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 0; m[4] = 0; m[5] = 2;
    MapPowers p2(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 0; m[5] = 0;
    MapPowers p3(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 1; m[5] = 0;
    MapPowers p4(6, m);

    CPPUNIT_ASSERT(!(p0 > p0));
    CPPUNIT_ASSERT( (p0 > p1));
    CPPUNIT_ASSERT( (p0 > p2));
    CPPUNIT_ASSERT( (p0 > p3));
    CPPUNIT_ASSERT( (p0 > p4));

    CPPUNIT_ASSERT(!(p1 > p0));
    CPPUNIT_ASSERT(!(p1 > p1));
    CPPUNIT_ASSERT( (p1 > p2));
    CPPUNIT_ASSERT( (p1 > p3));
    CPPUNIT_ASSERT( (p1 > p4));

    CPPUNIT_ASSERT(!(p2 > p0));
    CPPUNIT_ASSERT(!(p2 > p1));
    CPPUNIT_ASSERT(!(p2 > p2));
    CPPUNIT_ASSERT( (p2 > p3));
    CPPUNIT_ASSERT( (p2 > p4));

    CPPUNIT_ASSERT(!(p3 > p0));
    CPPUNIT_ASSERT(!(p3 > p1));
    CPPUNIT_ASSERT(!(p3 > p2));
    CPPUNIT_ASSERT(!(p3 > p3));
    CPPUNIT_ASSERT( (p3 > p4));

    CPPUNIT_ASSERT(!(p4 > p0));
    CPPUNIT_ASSERT(!(p4 > p1));
    CPPUNIT_ASSERT(!(p4 > p2));
    CPPUNIT_ASSERT(!(p4 > p3));
    CPPUNIT_ASSERT(!(p4 > p4));
  }

  void testTotalOrderingExampleLessEq(void) {
    IndexToPowerMap m;

    MapPowers p0(6, m);
    m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0; m[4] = 0; m[5] = 1;
    MapPowers p1(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 0; m[4] = 0; m[5] = 2;
    MapPowers p2(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 0; m[5] = 0;
    MapPowers p3(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 1; m[5] = 0;
    MapPowers p4(6, m);

    CPPUNIT_ASSERT( (p0 <= p0));
    CPPUNIT_ASSERT( (p1 <= p0));
    CPPUNIT_ASSERT( (p2 <= p0));
    CPPUNIT_ASSERT( (p3 <= p0));
    CPPUNIT_ASSERT( (p4 <= p0));

    CPPUNIT_ASSERT(!(p0 <= p1));
    CPPUNIT_ASSERT( (p1 <= p1));
    CPPUNIT_ASSERT( (p2 <= p1));
    CPPUNIT_ASSERT( (p3 <= p1));
    CPPUNIT_ASSERT( (p4 <= p1));

    CPPUNIT_ASSERT(!(p0 <= p2));
    CPPUNIT_ASSERT(!(p1 <= p2));
    CPPUNIT_ASSERT( (p2 <= p2));
    CPPUNIT_ASSERT( (p3 <= p2));
    CPPUNIT_ASSERT( (p4 <= p2));

    CPPUNIT_ASSERT(!(p0 <= p3));
    CPPUNIT_ASSERT(!(p1 <= p3));
    CPPUNIT_ASSERT(!(p2 <= p3));
    CPPUNIT_ASSERT( (p3 <= p3));
    CPPUNIT_ASSERT( (p4 <= p3));

    CPPUNIT_ASSERT(!(p0 <= p4));
    CPPUNIT_ASSERT(!(p1 <= p4));
    CPPUNIT_ASSERT(!(p2 <= p4));
    CPPUNIT_ASSERT(!(p3 <= p4));
    CPPUNIT_ASSERT( (p4 <= p4));
  }

  void testTotalOrderingExampleMoreEq(void) {
    IndexToPowerMap m;

    MapPowers p0(6, m);
    m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0; m[4] = 0; m[5] = 1;
    MapPowers p1(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 0; m[4] = 0; m[5] = 2;
    MapPowers p2(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 0; m[5] = 0;
    MapPowers p3(6, m);
    m[0] = 1; m[1] = 0; m[2] = 5; m[3] = 1; m[4] = 1; m[5] = 0;
    MapPowers p4(6, m);

    CPPUNIT_ASSERT(p0 >= p0);
    CPPUNIT_ASSERT(p0 >= p1);
    CPPUNIT_ASSERT(p0 >= p2);
    CPPUNIT_ASSERT(p0 >= p3);
    CPPUNIT_ASSERT(p0 >= p4);

    CPPUNIT_ASSERT(!(p1 >= p0));
    CPPUNIT_ASSERT( (p1 >= p1));
    CPPUNIT_ASSERT( (p1 >= p2));
    CPPUNIT_ASSERT( (p1 >= p3));
    CPPUNIT_ASSERT( (p1 >= p4));

    CPPUNIT_ASSERT(!(p2 >= p0));
    CPPUNIT_ASSERT(!(p2 >= p1));
    CPPUNIT_ASSERT( (p2 >= p2));
    CPPUNIT_ASSERT( (p2 >= p3));
    CPPUNIT_ASSERT( (p2 >= p4));

    CPPUNIT_ASSERT(!(p3 >= p0));
    CPPUNIT_ASSERT(!(p3 >= p1));
    CPPUNIT_ASSERT(!(p3 >= p2));
    CPPUNIT_ASSERT( (p3 >= p3));
    CPPUNIT_ASSERT( (p3 >= p4));

    CPPUNIT_ASSERT(!(p4 >= p0));
    CPPUNIT_ASSERT(!(p4 >= p1));
    CPPUNIT_ASSERT(!(p4 >= p2));
    CPPUNIT_ASSERT(!(p4 >= p3));
    CPPUNIT_ASSERT( (p4 >= p4));
  }

  void testDefaultGreaterEqualSelf(void) {
    MapPowers m;
    CPPUNIT_ASSERT(m >= m);
    CPPUNIT_ASSERT(m <= m);
  }

  void testExampleGreaterEqualSelf(void) {
    MapPowers m(999);
    m.setPower(421, 30);
    m.setPower(998, 2);
    MapPowers n(999);
    CPPUNIT_ASSERT(n != m);
    CPPUNIT_ASSERT(m.getNumVars() == 999);
    CPPUNIT_ASSERT(m.degree() == 32);
    CPPUNIT_ASSERT(m >= m);
    CPPUNIT_ASSERT(m <= m);
  }

  void testNotGreaterSelf(void) {
    MapPowers m;
    CPPUNIT_ASSERT(!(m < m));
    CPPUNIT_ASSERT(!(m > m));
    MapPowers n(999);
    CPPUNIT_ASSERT(!(n < n));
    CPPUNIT_ASSERT(!(n > n));
    n.setPower(421, 30);
    n.setPower(998, 2);
    CPPUNIT_ASSERT(!(n < n));
    CPPUNIT_ASSERT(!(n > n));
  }

  void testSetZero(void) {
    MapPowers m;
    m.setPower(0, 0);
    CPPUNIT_ASSERT(m.getNumPowersStored() == 0);
  }

  void testSetPositive(void) {
    MapPowers m;
    m.setPower(0, 1);
    CPPUNIT_ASSERT(m.getNumPowersStored() == 1);
    CPPUNIT_ASSERT(m[0] == 1);
    CPPUNIT_ASSERT(m.degree() == 1);
  }

#ifdef SIGNED_POWER_TYPE
  void testSetNegative(void) {
    MapPowers m;
    m.setPower(0, -1);
  }
#endif //SIGNED_POWER_TYPE

  void testIncompatibleEq(void) {
    MapPowers m;
    MapPowers n(3);
    m == n; //should throw index error
  }

  void testIncompatibleNe(void) {
    MapPowers m;
    MapPowers n(3);
    m != n; //should throw index error
  }

  void testIncompatibleLe(void) {
    MapPowers m;
    MapPowers n(3);
    m <= n; //should throw index error
  }

  void testIncompatibleGe(void) {
    MapPowers m;
    MapPowers n(3);
    m >= n; //should throw index error
  }

  void testIncompatibleLt(void) {
    MapPowers m;
    MapPowers n(3);
    m < n; //should throw index error
  }

  void testIncompatibleGt(void) {
    MapPowers m;
    MapPowers n(3);
    m > n; //should throw index error
  }

  void testDefaultCallWrongLen(void) {
    MapPowers m;
    std::vector< Coeff > v;
    v.push_back(Coeff(3.13));
    v.push_back(Coeff(2.2));
    m(v);
  }

  void testDefaultCallAllZero(void) {
    MapPowers m;
    std::vector< Coeff > v;
    v.push_back(CoeffZero);
    CPPUNIT_ASSERT(m(v) == CoeffOne);
  }

  void testExampleCallWrongLen(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 0);
    std::vector< Coeff > v;
    v.push_back(Coeff(3.13));
    v.push_back(Coeff(2.2));
    v.push_back(CoeffZero);
    v.push_back(Coeff(-5.7));
    m(v);
  }

  void testExampleCallAllZero(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 0);
    std::vector< Coeff > v(3);
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    CPPUNIT_ASSERT(v.size() == 3);
    CPPUNIT_ASSERT(v[0] == CoeffZero);
    CPPUNIT_ASSERT(v[1] == CoeffZero);
    CPPUNIT_ASSERT(v[2] == CoeffZero);
    CPPUNIT_ASSERT(m(v) == CoeffZero);
  }

  void testExampleCallValues(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 0);
    std::vector< Coeff > v(3);
    v[0] = -5.7;
    v[1] = +3.13;
    v[2] = +2.2;
    CPPUNIT_ASSERT(v.size() == 3);
    CPPUNIT_ASSERT(v[0] == Coeff(-5.7));
    CPPUNIT_ASSERT(v[1] == Coeff(+3.13));
    CPPUNIT_ASSERT(v[2] == Coeff(+2.2));
    const Real tol(1.0e-14); //will work for both double and mpf
    CPPUNIT_ASSERT(fabs(m(v) - Coeff(-5.7*3.13*3.13)) < tol);
  }

  void testExampleCallAllOne(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 0);
    std::vector< Coeff > v(3);
    v[0] = 1.0;
    v[1] = 1.0;
    v[2] = 1.0;
    CPPUNIT_ASSERT(m(v) == CoeffOne);
  }

  void testExampleCall_anyZero(void) {
    MapPowers m(3);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 7);
    std::vector< Coeff > v(3);
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
	v[j] = 1.0;
      }
      CPPUNIT_ASSERT(m(v) == CoeffOne);
      v[i] = 0.0;
      CPPUNIT_ASSERT(m(v) == CoeffZero);
    }
  }

  void testExampleCallNegativeCounts(void) {
    MapPowers m(6);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 7);
    m.setPower(3, 19);
    m.setPower(4, 0);
    m.setPower(5, 3);
    CPPUNIT_ASSERT(m.degree() == 1+2+7+19+3);
    std::vector< Coeff > ones(6);
    for (int i=0; i<6; ++i) {
      ones[i] = 1.0;
    }
    CPPUNIT_ASSERT(m(ones) == CoeffOne);
    std::vector< Coeff > alts(6);
    for (int i=0; i<6; ++i) {
      alts[i] = (i%2) ? +1.0 : -1.0;
    }
    CPPUNIT_ASSERT(m(alts) == CoeffOne);
    float z = 1;
    for (int i=0; i<6; ++i) {
      alts[i] = (i%2) ? -1.0 : +1.0;
      if ((m[i]%2) && (i%2)) {
	z *= -1;
      }
    }
    CPPUNIT_ASSERT(m(alts) == Coeff(z));
  }

  void testExamplePowZeroIsOne(void) {
    MapPowers m(6);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 7);
    m.setPower(3, 19);
    m.setPower(4, 0);
    m.setPower(5, 3);
    CPPUNIT_ASSERT(m.pow(0) == MapPowers(6));
  }

  void testExamplePowOneIsSame(void) {
    MapPowers m(6);
    const MapPowers n(m);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 7);
    m.setPower(3, 19);
    m.setPower(4, 0);
    m.setPower(5, 3);
    CPPUNIT_ASSERT(m.pow(1) == m);
    CPPUNIT_ASSERT(m.pow(1) == n);
  }

  void testExamplePow(void) {
    MapPowers m(6);
    m.setPower(0, 1);
    m.setPower(1, 2);
    m.setPower(2, 7);
    m.setPower(3, 19);
    m.setPower(4, 0);
    m.setPower(5, 3);
    MapPowers n(6);
    for (int p; p<10; ++p) {
      n = m.pow(p);
      for (int i; i<6; ++i) {
	CPPUNIT_ASSERT(n[i] == m[i]*p);
      }
    }
  }

#ifdef SIGNED_POWER_TYPE
  void testNegativeIndexPow(void) {
    MapPowers m(99);
    m.pow(-1);
  }
#endif //SIGNED_POWER_TYPE

  void testMul(void) {
    for (Index length=1; length<MPTLengths; ++length) {
      MapPowers z(length);
      for (size_t ncase=0; ncase<MPTNumCases; ++ncase) {
	std::vector< Power > pt = randomVectorPowers(length, MPTMaxDegree);
	MapPowers p(pt);
	std::vector< Power > qt = randomVectorPowers(length, MPTMaxDegree);
	MapPowers q(qt);
	std::vector< Power > rt(length);
	for (Index i=0; i<length; ++i) {
	  rt[i] = pt[i] + qt[i];
	}
	MapPowers r(rt);
	CPPUNIT_ASSERT(p*q == r);
	CPPUNIT_ASSERT(q*p == r);
	CPPUNIT_ASSERT(q*z == q);
	CPPUNIT_ASSERT(p*z == p);
	CPPUNIT_ASSERT(z*p == p);
	CPPUNIT_ASSERT(z*q == q);
	IndexToPowerMap d;
	for (Index i=0; i<length; ++i) {
	  if (rt[i] > 0) {
	    d[i] = rt[i];
	  }
	}
	CPPUNIT_ASSERT(MapPowers(length, d) == p*q);
      }
    }
  }

  void testDiff(void) {
    for (Index length=1; length<MPTLengths; ++length) {
      MapPowers z(length);
      for (size_t ncase=0; ncase<MPTNumCases; ++ncase) {
	std::vector< Power > pt = randomVectorPowers(length, MPTMaxDegree);
	const MapPowers p(pt);
	for (Index c=0; c<length; ++c) {
	  std::pair< Coeff, MapPowers > c_d = p.diff(c);
	  if (p[c] > 0) {
	    CPPUNIT_ASSERT(c_d.first == Coeff(p[c]));
	    for (Index d=0; d<length; ++d)
              if (d == c)
                CPPUNIT_ASSERT((c_d.second)[d] == (p[c]-1));
              else
                CPPUNIT_ASSERT((c_d.second)[d] == p[d]);
	  }
          else
	    CPPUNIT_ASSERT(c_d.first == CoeffZero);
	}
      }
    }
  }

#ifdef SIGNED_POWER_TYPE

  void testDiffNegativeIndex1(void) {
    MapPowers p;
    p.diff(-1); //throws index error
  }

  void testDiffNegativeIndex2(void) {
    MapPowers p;
    p.diff(-10); //throws index error
  }

#endif //SIGNED_POWER_TYPE

  void testDiffBadIndex1(void) {
    MapPowers p;
    p.diff(1); //throws index error
  }

  void testDiffBadIndex2(void) {
    MapPowers p(1);
    p.diff(2); //throws index error
  }

}; //MapPowersTest

#endif //MAP_POWERS_TEST_H
