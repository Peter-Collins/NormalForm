//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef POLYNOMIAL_TEST_OTHERS_H
#define POLYNOMIAL_TEST_OTHERS_H

//system
#include <utility> //for pair
#include <set>

//library
#include <cppunit/extensions/HelperMacros.h>

//project
#include "Types.h"
#include "StlInit.h"
#include "MapPowers.h"
#include "Polynomial.h"

//constants
static const size_t PolyOtherNumCases = 32;
static const Index PolyOtherMaxNumVars = 25;
static const Real PolyOtherTolerance = 1.0e-12;

inline bool isOddDegree(const MapPowers& p, const Coeff& c) {
  return (p.degree())%2;
}

class PolynomialTestOthers : public CppUnit::TestFixture {

 public:

  //The following create static CppUnit::TestSuite *suite():
  CPPUNIT_TEST_SUITE(PolynomialTestOthers);

  //Structural integrity:
  CPPUNIT_TEST(testInPlaceAddSelf);
  CPPUNIT_TEST(testAddMonomials);
  CPPUNIT_TEST(testSetMonomials);
  CPPUNIT_TEST(testAddHomogeneous);
  CPPUNIT_TEST(testCallSum);
  CPPUNIT_TEST(testHasTerms);
  CPPUNIT_TEST(testPowersAndCoefficients);

  //Homogeneous structure:
  CPPUNIT_TEST(testHomogeneousZeroIs);
  CPPUNIT_TEST(testHomogeneousOneIs);
  CPPUNIT_TEST(testHomogeneousExampleIs);
  CPPUNIT_TEST(testHomogeneousExampleIsNot);

  //Logical/contract integrity:
  CPPUNIT_TEST(testMonomialsNeverRepeated);
  CPPUNIT_TEST(testCoeffsNeverZero);
  CPPUNIT_TEST(testCoeffsSetZeroErasesTerm);
  //CPPUNIT_TEST(testMulByTinyCoeffErasesTerms);

  //Coordinate monomials:
  CPPUNIT_TEST(testCoordinateMonomialsCorrect);

  //Partitions by predicate:
  CPPUNIT_TEST(testPartitionOddsAndEvens);
  CPPUNIT_TEST(testPartitionOddsAndEvensTermsPresent);

  //Pow:
  CPPUNIT_TEST(testPowZeroGivesOne);
  CPPUNIT_TEST(testPowOneGivesSame);
  CPPUNIT_TEST(testPowTwoGivesMulSelf);
  CPPUNIT_TEST(testPowNGivesMulSelfN);
  CPPUNIT_TEST(testPowAssociative);

  //End the test suite:
  CPPUNIT_TEST_SUITE_END();

  void testAddMonomials(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = a;
      PowersToCoeffMap::const_iterator pc;
      Polynomial rhs(vars);
      for (pc = a.getPowersAndCoeffs().begin();
	   pc != a.getPowersAndCoeffs().end();
	   ++pc) {
	rhs += Polynomial(pc->first, pc->second);
      }
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testInPlaceAddSelf(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = a;
      lhs += a;
      Polynomial rhs(a + a);
      PowersToCoeffMap::const_iterator pc;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testSetMonomials(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      Polynomial a = randomPoly(vars);
      Polynomial lhs = a;
      PowersToCoeffMap::const_iterator pc;
      Polynomial rhs(vars);
      for (pc = a.getPowersAndCoeffs().begin();
	   pc != a.getPowersAndCoeffs().end();
	   ++pc) {
	rhs.setMonomial(pc->first, pc->second);
      }
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testAddHomogeneous(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      const Polynomial lhs = a;
      Polynomial rhs(vars);
      for (Power deg=0; deg<a.degree()+1; ++deg) {
	rhs += a.homogeneousPart(deg);
      }
      CPPUNIT_ASSERT(rhs.degree() == lhs.degree());
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testCallSum(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      std::vector< Coeff > ones;
      for (Index i=0; i<vars; ++i) {
	init(ones) += CoeffOne;
      }
      CPPUNIT_ASSERT(ones.size() == vars);
      Coeff lhs = a(ones);
      Coeff rhs(0.0);
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	rhs += pc->second; //coefficient
      }
      CPPUNIT_ASSERT(fabs(rhs-lhs) < PolyOtherTolerance);
    }
  }

  void testHasTerms(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(a.hasPowers(pc->first));
      }
    }
  }

  void testPowersAndCoefficients(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial lhs = randomPoly(vars);
      Polynomial rhs(vars);
      PowersToCoeffMap::const_iterator pc;
      for (pc=lhs.getPowersAndCoeffs().begin();
	   pc!=lhs.getPowersAndCoeffs().end();
	   ++pc) {
	rhs.setMonomial(pc->first, pc->second);
      }
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testPartitionOddsAndEvens(void) { 
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial lhs = randomPoly(vars);
      std::pair< Polynomial, Polynomial > oddEven = lhs.partition(&isOddDegree);
      PowersToCoeffMap::const_iterator pc;
      for (pc=oddEven.first.getPowersAndCoeffs().begin();
	   pc!=oddEven.first.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(((pc->first).degree())%2);
      }
      for (pc=oddEven.second.getPowersAndCoeffs().begin();
	   pc!=oddEven.second.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(!(((pc->first).degree())%2));
      }
      Polynomial rhs(vars);
      rhs += oddEven.first;
      rhs += oddEven.second;
      CPPUNIT_ASSERT((rhs-lhs).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testPartitionOddsAndEvensTermsPresent(void) { 
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial lhs = randomPoly(vars);
      std::pair< Polynomial, Polynomial > oddEven = lhs.partition(&isOddDegree);
      PowersToCoeffMap::const_iterator pc;
      for (pc=lhs.getPowersAndCoeffs().begin();
	   pc!=lhs.getPowersAndCoeffs().end();
	   ++pc) {
	if (oddEven.first.hasPowers(pc->first)) {
	  CPPUNIT_ASSERT(!(oddEven.second.hasPowers(pc->first)));
	} else {
	  CPPUNIT_ASSERT(oddEven.second.hasPowers(pc->first));
	}
      }
    }
  }

  void testMonomialsNeverRepeated(void) { 
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      std::set< MapPowers > seen;
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(seen.find(pc->first) == seen.end());
	seen.insert(pc->first);
      }
    }
  }

  void testCoeffsNeverZero(void) {
    //this is not a great test, since we are unlikely to set 0.0 via rand.
    //therefore, we deliberately set a zero in a copy.
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	CPPUNIT_ASSERT(pc->second != CoeffZero);
      }
    }
  }

  void testCoeffsSetZeroErasesTerm(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      PowersToCoeffMap::const_iterator pc;
      for (pc=a.getPowersAndCoeffs().begin();
	   pc!=a.getPowersAndCoeffs().end();
	   ++pc) {
	Polynomial b(a);
	b.setMonomial(pc->first, CoeffZero);
	CPPUNIT_ASSERT(!(b.hasPowers(pc->first)));
      }
    }
  }

/*   void testMulByTinyCoeffErasesTerms(void) { */
/*     for (size_t n = 0; n < PolyOtherNumCases; ++n) { */
/*       const Index vars(1 + randomIndex(PolyOtherMaxNumVars)); */
/*       Polynomial a(randomPoly(vars)); */
/*       a *= Coeff(1.0e-12); */
/*       a *= Coeff(1.0e-12); */
/*       a *= Coeff(1.0e-12); */
/*       a *= Coeff(1.0e-12); */
/*       CPPUNIT_ASSERT(a.isZero()); */
/*     } */
/*   } */

  void testPowZeroGivesOne(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      const Polynomial one = Polynomial::One(vars);
      CPPUNIT_ASSERT(a.pow(0) == one);
    }
  }

  void testPowOneGivesSame(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      const Polynomial b(a);
      CPPUNIT_ASSERT(a.pow(1) == b);
    }
  }

  void testPowTwoGivesMulSelf(void) {
    for (size_t n=0; n<PolyOtherNumCases; ++n) {
      Index vars = 1+randomIndex(PolyOtherMaxNumVars);
      const Polynomial a = randomPoly(vars);
      const Polynomial b(a);
      CPPUNIT_ASSERT((a.pow(2) - b*b).lInfinityNorm() < PolyOtherTolerance);
    }
  }

  void testPowNGivesMulSelfN(void) {
    //this is an expensive test to run
    for (size_t n=0; n<PolyOtherNumCases/2; ++n) {
      Index vars = 1+randomIndex(4);
      const Polynomial a = randomPoly(vars);
      Polynomial b(a);
      for (Power p=1; p<4; ++p) {
	CPPUNIT_ASSERT((a.pow(p) - b).lInfinityNorm() < PolyOtherTolerance);
	b *= a;
      }
    }
  }

  void testPowAssociative(void) {
    //this is an expensive test to run
    for (size_t n=0; n<PolyOtherNumCases/2; ++n) {
      Index vars = 1+randomIndex(4);
      const Polynomial a = randomPoly(vars);
      for (Power p=1; p<4; ++p) {
	Polynomial b = a.pow(p-1);
	Polynomial c = a.pow(p);
	CPPUNIT_ASSERT((a*b - c).lInfinityNorm() < PolyOtherTolerance);
      }
    }
  }

  void testHomogeneousZeroIs(void) {
    Polynomial p(5);
    CPPUNIT_ASSERT(p.isHomogeneous() );
    CPPUNIT_ASSERT(p.isHomogeneous(0));
    CPPUNIT_ASSERT(p.isHomogeneous(1));
    CPPUNIT_ASSERT(p.isHomogeneous(2));
    CPPUNIT_ASSERT(p.isHomogeneous(3));
    CPPUNIT_ASSERT(p.isHomogeneous(4));
    CPPUNIT_ASSERT(p.isHomogeneous(5));
    CPPUNIT_ASSERT(p.isHomogeneous(6));
  }

  void testHomogeneousOneIs(void) {
    Polynomial p = Polynomial::One(5);
    CPPUNIT_ASSERT(p.isHomogeneous());
    CPPUNIT_ASSERT(p.isHomogeneous(0));
    CPPUNIT_ASSERT(!(p.isHomogeneous(1)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(2)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(3)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(4)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(5)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(6)));
  }

  void testHomogeneousExampleIs(void) {
    Polynomial p(2);
    std::vector< Power > powers;
    init(powers) = 1, 2;
    p.setMonomial(powers, Coeff(0.5));
    init(powers) = 1, 1;
    p.setMonomial(powers, CoeffZero);
    CPPUNIT_ASSERT(p.isHomogeneous());
    CPPUNIT_ASSERT(!(p.isHomogeneous(0)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(1)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(2)));
    CPPUNIT_ASSERT(p.isHomogeneous(3));
    CPPUNIT_ASSERT(!(p.isHomogeneous(4)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(5)));
    CPPUNIT_ASSERT(!(p.isHomogeneous(6)));
  }

  void testHomogeneousExampleIsNot(void) {
    Polynomial p(2);
    std::vector< Power > powers;
    init(powers) = 1, 2;
    p.setMonomial(powers, Coeff(0.5));
    init(powers) = 1, 1;
    p.setMonomial(powers, Coeff(0.1));
    CPPUNIT_ASSERT(!p.isHomogeneous());
    CPPUNIT_ASSERT(!p.isHomogeneous(0));
    CPPUNIT_ASSERT(!p.isHomogeneous(1));
    CPPUNIT_ASSERT(!p.isHomogeneous(2));
    CPPUNIT_ASSERT(!p.isHomogeneous(3));
    CPPUNIT_ASSERT(!p.isHomogeneous(4));
    CPPUNIT_ASSERT(!p.isHomogeneous(5));
    CPPUNIT_ASSERT(!p.isHomogeneous(6));
  }

  void testCoordinateMonomialsCorrect(void) {
    for (Index dim=1; dim<10; ++dim) {
      for (Index i=0; i<dim; ++i) {
	Polynomial m = Polynomial::CoordinateMonomial(dim, i);
	CPPUNIT_ASSERT(m.getNumTerms() == 1);
	PowersToCoeffMap::const_iterator pc;
	pc = m.getPowersAndCoeffs().begin();
	CPPUNIT_ASSERT(pc->first.degree() == 1);
	CPPUNIT_ASSERT(pc->second == CoeffOne);
	for (Index j=0; j<dim; ++j) {
	  if (i == j) {
	    CPPUNIT_ASSERT(pc->first[j] == 1);
	  } else {
	    CPPUNIT_ASSERT(pc->first[j] == 0);
	  }
	}
      }
    }
  }

}; //PolynomialTestOthers

#endif //POLYNOMIAL_TEST_OTHERS_H
