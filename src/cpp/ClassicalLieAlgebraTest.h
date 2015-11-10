#ifndef CLASSICAL_LIE_ALGEBRA_TEST_H
#define CLASSICAL_LIE_ALGEBRA_TEST_H

/// @file ClassicalLieAlgebraTest.h
///
/// @brief Test fixture for Classical version of Lie algebra
///
/// Performs tests of construction, copy construction, helper methods
/// for the creation of polynomials within the algebra, the embedding
/// of the q (position) and p (conjugate momentum) components into the
/// polynomial ring, and the graded structure, including methods for
/// extracting isograde components, etc.
///
/// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

//system headers
#include <set>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "LieAlgebraBase.h"
#include "ClassicalLieAlgebra.h"
#include "LieAlgebraTestBase.h"

///
/// NOTE: we must use =, not () to initialize static data members!  DO
/// NOT PUT THESE DECLARATIONS IN A CLASS DEFINITION!
///
static const Index CLAMaxDof = 50;
static const Power CLAMaxPower = 8;

/// test fixture class
class ClassicalLieAlgebraTest : public CppUnit::TestFixture {

 public:

  /// create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(ClassicalLieAlgebraTest);

  //constructor from number of variables
  CPPUNIT_TEST_EXCEPTION(testConstructorGivenZeroVars, LieAlgebraIndexError);
  CPPUNIT_TEST(testConstructorGivenNumVars);
  CPPUNIT_TEST_EXCEPTION(testConstructorGivenZeroVarsX, LieAlgebraIndexError);
  CPPUNIT_TEST(testConstructorGivenNumVarsX);
  
  //big four
  CPPUNIT_TEST(testDefaultConstructor);
  CPPUNIT_TEST(testCopyConstructor);
  CPPUNIT_TEST(testDefaultConstructorX);
  CPPUNIT_TEST(testCopyConstructorX);

  //basic object creation
  CPPUNIT_TEST(testZeroPolynomial);
  CPPUNIT_TEST(testOnePolynomial);
  CPPUNIT_TEST(testZeroPolynomialX);
  CPPUNIT_TEST(testOnePolynomialX);

  CPPUNIT_TEST(testCoordinateMonomialQ);
  CPPUNIT_TEST(testCoordinateMonomialP);
  CPPUNIT_TEST(testCoordinateMonomialQX);
  CPPUNIT_TEST(testCoordinateMonomialPX);

  CPPUNIT_TEST(testCoordinateMonomialQPower);
  CPPUNIT_TEST(testCoordinateMonomialPPower);
  CPPUNIT_TEST(testCoordinateMonomialQPowerX);
  CPPUNIT_TEST(testCoordinateMonomialPPowerX);

  //test the q, p encoding
  CPPUNIT_TEST(testNoIndiciesAreInConflict);
  CPPUNIT_TEST(testAllIndiciesAreUsed);
  CPPUNIT_TEST(testNoIndiciesAreInConflictX);
  CPPUNIT_TEST(testAllIndiciesAreUsedX);

  //graded structure
  CPPUNIT_TEST(testGradeZero);
  CPPUNIT_TEST(testGradeOne);
  CPPUNIT_TEST(testGradeCoordinates);
  CPPUNIT_TEST(testGradePurePowers);
  CPPUNIT_TEST(testSumOfIsoGradesIsPolynomial);
  CPPUNIT_TEST(testSumOfIsoGradesIsPolynomialExample);

  /// end the test suite
  CPPUNIT_TEST_SUITE_END();

 protected:

  /// attempt to construct algebra over zero dof should throw an index error
  void testConstructorGivenZeroVarsX(void) {
    ClassicalLieAlgebra alg(0);
  }
  /// @brief construct algebras with given number of degrees of freedom.
  ///
  /// test their relationship to the corresponding polynomial ring
  void testConstructorGivenNumVarsX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof);
      Polynomial p(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(p));
      CPPUNIT_ASSERT(alg.grade(p) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(p));
      CPPUNIT_ASSERT(alg.isoGradePart(p, 0) == p);
    }
  }
  /// test copy construction of an algegra
  void testCopyConstructorX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      ClassicalLieAlgebra algCopy(alg);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof);
      CPPUNIT_ASSERT(algCopy.getDof() == dof);
      CPPUNIT_ASSERT(algCopy.getNumVars() == 2*dof);
    }
  }
  /// test that the default constructor gives a one degree of freedom algebra
  void testDefaultConstructorX(void) {
    ClassicalLieAlgebra alg;
    CPPUNIT_ASSERT(alg.getDof() == 1);
    CPPUNIT_ASSERT(alg.getNumVars() == 2);
    Polynomial p(2);
    CPPUNIT_ASSERT(alg.hasElt(p));
    CPPUNIT_ASSERT(alg.grade(p) == 0);
    CPPUNIT_ASSERT(alg.isIsoGrade(p));
    CPPUNIT_ASSERT(alg.isoGradePart(p, 0) == p);
  }
  /// @brief thoroughly test properties of zero polynomial.
  ///
  /// Test the zero polynomial from the algebra and the one given by
  /// explicit creation of a zero polynomial in the polynomial ring
  /// against one another: test for compatibility of the corresponding
  /// algebras with the zero polynomials. note that taking the grade
  /// of the zero polynomial will return zero, but that asking whether
  /// it is isograde of grade g will be true for all g; this behaviour
  /// is in accordance with everyday mathematical convention.
  void testZeroPolynomialX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      Polynomial algZero = alg.zero();
      Polynomial polyZero = Polynomial::Zero(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(algZero));
      CPPUNIT_ASSERT(alg.hasElt(polyZero));
      CPPUNIT_ASSERT(algZero.getNumVars() == 2*alg.getDof());
      CPPUNIT_ASSERT(algZero.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(algZero.isZero());
      CPPUNIT_ASSERT(algZero.isConstant());
      CPPUNIT_ASSERT(algZero.isHomogeneous());
      CPPUNIT_ASSERT(algZero.isHomogeneous(0));
      CPPUNIT_ASSERT(alg.grade(algZero) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(algZero));
      CPPUNIT_ASSERT(alg.isoGradePart(algZero, 0) == algZero);
      CPPUNIT_ASSERT(algZero == polyZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffZero)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffOne)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffJ)) == CoeffZero);
    }
  }
  /// @brief test properties of the constant one polynomial.
  ///
  /// and compatibility of the algebras with it.
  void testOnePolynomialX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      Polynomial algOne = alg.one();
      Polynomial polyOne = Polynomial::One(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(algOne));
      CPPUNIT_ASSERT(alg.hasElt(polyOne));
      CPPUNIT_ASSERT(algOne.getNumVars() == 2*alg.getDof());
      CPPUNIT_ASSERT(algOne.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(!algOne.isZero());
      CPPUNIT_ASSERT(algOne.isConstant());
      CPPUNIT_ASSERT(algOne.isHomogeneous());
      CPPUNIT_ASSERT(algOne.isHomogeneous(0));
      CPPUNIT_ASSERT(alg.grade(algOne) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(algOne));
      CPPUNIT_ASSERT(alg.isoGradePart(algOne, 0) == algOne);
      CPPUNIT_ASSERT(algOne == polyOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffZero)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffOne)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffJ)) == CoeffOne);
    }
  }
  /// test the creation of a coordinate monomial and its grade properties
  void testCoordinateMonomialQX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial qI = alg.q(i);
        Polynomial xQ = alg.coordinateMonomial(alg.iQ(i));
        CPPUNIT_ASSERT(qI == xQ);
        CPPUNIT_ASSERT(alg.hasElt(qI));
        CPPUNIT_ASSERT(qI.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(!qI.isZero());
        CPPUNIT_ASSERT(!qI.isConstant());
        CPPUNIT_ASSERT(qI.isHomogeneous());
        CPPUNIT_ASSERT(qI.isHomogeneous(1));
        CPPUNIT_ASSERT(qI(vec) == Coeff(i));
        CPPUNIT_ASSERT(alg.grade(qI) == 1);
        CPPUNIT_ASSERT(alg.isIsoGrade(qI));
        CPPUNIT_ASSERT(alg.isoGradePart(qI, 1) == qI);
      }
    }
  }
  /// test the creation of a momentum monomial and its grade properties
  void testCoordinateMonomialPX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial pI = alg.p(i);
        Polynomial xP = alg.coordinateMonomial(alg.iP(i));
        CPPUNIT_ASSERT(pI == xP);
        CPPUNIT_ASSERT(alg.hasElt(pI));
        CPPUNIT_ASSERT(pI.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(!pI.isZero());
        CPPUNIT_ASSERT(!pI.isConstant());
        CPPUNIT_ASSERT(pI.isHomogeneous());
        CPPUNIT_ASSERT(pI.isHomogeneous(1));
        CPPUNIT_ASSERT(pI(vec) == Coeff(i));
        CPPUNIT_ASSERT(alg.grade(pI) == 1);
        CPPUNIT_ASSERT(alg.isIsoGrade(pI));
        CPPUNIT_ASSERT(alg.isoGradePart(pI, 1) == pI);
      }
    }
  }
  /// test higher powers in single coordinate monomials
  void testCoordinateMonomialQPowerX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial qI = alg.q(i);
        Polynomial xQ = alg.coordinateMonomial(alg.iQ(i));
        for (Power p = 0; p < CLAMaxPower; ++p) {
          Polynomial qIP = alg.q(i, p);
          Polynomial xQP = alg.coordinateMonomial(alg.iQ(i), p);
          CPPUNIT_ASSERT(alg.hasElt(qIP));
          CPPUNIT_ASSERT(alg.hasElt(xQP));
          CPPUNIT_ASSERT(qIP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(xQP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(!qIP.isZero());
          CPPUNIT_ASSERT(!xQP.isZero());
          if (p == 0) {
            CPPUNIT_ASSERT(qIP.isConstant());
            CPPUNIT_ASSERT(xQP.isConstant());
          }
          else {
            CPPUNIT_ASSERT(!qIP.isConstant());
            CPPUNIT_ASSERT(!xQP.isConstant());
          }
          CPPUNIT_ASSERT(qIP.isHomogeneous());
          CPPUNIT_ASSERT(xQP.isHomogeneous());
          CPPUNIT_ASSERT(qIP.isHomogeneous(p));
          CPPUNIT_ASSERT(xQP.isHomogeneous(p));
          CPPUNIT_ASSERT(qIP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(xQP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(qIP == qI.pow(p));
          CPPUNIT_ASSERT(xQP == xQ.pow(p));
          CPPUNIT_ASSERT(alg.grade(qIP) == p);
          CPPUNIT_ASSERT(alg.grade(xQP) == p);
          CPPUNIT_ASSERT(alg.isIsoGrade(qIP));
          CPPUNIT_ASSERT(alg.isIsoGrade(xQP));
          CPPUNIT_ASSERT(alg.isoGradePart(qIP, p) == qIP);
          CPPUNIT_ASSERT(alg.isoGradePart(xQP, p) == xQP);
        }
      }
    }
  }
  /// test higher powers in single momentum monomials
  void testCoordinateMonomialPPowerX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial pI = alg.p(i);
        Polynomial xP = alg.coordinateMonomial(alg.iP(i));
        for (Power p = 0; p < CLAMaxPower; ++p) {
          Polynomial pIP = alg.p(i, p);
          Polynomial xPP = alg.coordinateMonomial(alg.iP(i), p);
          CPPUNIT_ASSERT(alg.hasElt(pIP));
          CPPUNIT_ASSERT(alg.hasElt(xPP));
          CPPUNIT_ASSERT(pIP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(xPP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(!pIP.isZero());
          CPPUNIT_ASSERT(!xPP.isZero());
          if (p == 0) {
            CPPUNIT_ASSERT(pIP.isConstant());
            CPPUNIT_ASSERT(xPP.isConstant());
          }
          else {
            CPPUNIT_ASSERT(!pIP.isConstant());
            CPPUNIT_ASSERT(!xPP.isConstant());
          }
          CPPUNIT_ASSERT(pIP.isHomogeneous());
          CPPUNIT_ASSERT(xPP.isHomogeneous());
          CPPUNIT_ASSERT(pIP.isHomogeneous(p));
          CPPUNIT_ASSERT(xPP.isHomogeneous(p));
          CPPUNIT_ASSERT(pIP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(xPP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(pIP == pI.pow(p));
          CPPUNIT_ASSERT(xPP == xP.pow(p));
          CPPUNIT_ASSERT(alg.grade(pIP) == p);
          CPPUNIT_ASSERT(alg.grade(xPP) == p);
          CPPUNIT_ASSERT(alg.isIsoGrade(pIP));
          CPPUNIT_ASSERT(alg.isIsoGrade(xPP));
          CPPUNIT_ASSERT(alg.isoGradePart(pIP, p) == pIP);
          CPPUNIT_ASSERT(alg.isoGradePart(xPP, p) == xPP);
        }
      }
    }
  }
  /// @brief test injectivity of coordinate embedding
  ///
  /// ensure that the embedding of the coordinates and the momenta
  /// into a single polynomial ring is really an embedding; first
  /// step: ensure that no indices are in conflict (the index mapping
  /// is injective).
  void testNoIndiciesAreInConflictX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      for (Index d = 0; d < dof; ++d)
        for (Index e = 0; e < dof; ++e) {
          Index iQ = alg.iQ(d);
          Index iP = alg.iP(e);
          CPPUNIT_ASSERT(0 <= iQ && iQ < alg.getNumVars());
          CPPUNIT_ASSERT(0 <= iP && iP < alg.getNumVars());
          CPPUNIT_ASSERT(iQ != iP);
          Polynomial qD = alg.q(d);
          Polynomial pE = alg.p(e);
          CPPUNIT_ASSERT(qD != pE);
        }
    }
  }
  /// @brief test that coordinate embedding is onto
  ///
  /// ensure that the embedding of the coordinates and the momenta
  /// into a single polynomial ring is really an embedding; first
  /// step: ensure that all indices are used (the index mapping is
  /// onto).
  void testAllIndiciesAreUsedX(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::set< Index > available;
      for (Index i = 0; i < alg.getNumVars(); ++i)
        available.insert(i);
      CPPUNIT_ASSERT(available.size() == alg.getNumVars());
      for (Index d = 0; d < dof; ++d) {
        Index iQ = alg.iQ(d);
        CPPUNIT_ASSERT(0 <= iQ && iQ < alg.getNumVars());
        CPPUNIT_ASSERT(available.erase(iQ) == 1);
      }
      for (Index e = 0; e < dof; ++e) {
        Index iP = alg.iP(e);
        CPPUNIT_ASSERT(0 <= iP && iP < alg.getNumVars());
        CPPUNIT_ASSERT(available.erase(iP) == 1);
      }
      CPPUNIT_ASSERT(available.size() == 0);
    }
  }
  /// @brief zero dof algebras are not allowed
  ///
  /// construction of an algebra with zero degrees of freedom should
  /// throw an index error
  void testConstructorGivenZeroVars(void) {
    ClassicalLieAlgebra alg(0);
  }
  /// construct algs with given number of variables
  void testConstructorGivenNumVars(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof);
      Polynomial p(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(p));
    }
  }
  /// construct a copy of an algebra
  void testCopyConstructor(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      ClassicalLieAlgebra algCopy(alg);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof);
      CPPUNIT_ASSERT(algCopy.getDof() == dof);
      CPPUNIT_ASSERT(algCopy.getNumVars() == 2*dof);
    }
  }
  ///test that the default constructor gives a one degree of freedom algebra
  void testDefaultConstructor(void) {
    ClassicalLieAlgebra alg;
    CPPUNIT_ASSERT(alg.getDof() == 1);
    CPPUNIT_ASSERT(alg.getNumVars() == 2);
    Polynomial p(2);
    CPPUNIT_ASSERT(alg.hasElt(p));
  }
  /// properties of zero polynomial and compatibility with poly zero (see above)
  void testZeroPolynomial(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      Polynomial algZero = alg.zero();
      Polynomial polyZero = Polynomial::Zero(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(algZero));
      CPPUNIT_ASSERT(alg.hasElt(polyZero));
      CPPUNIT_ASSERT(algZero.getNumVars() == 2*alg.getDof());
      CPPUNIT_ASSERT(algZero.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(algZero.isZero());
      CPPUNIT_ASSERT(algZero.isConstant());
      CPPUNIT_ASSERT(algZero.isHomogeneous());
      CPPUNIT_ASSERT(algZero.isHomogeneous(0));
      CPPUNIT_ASSERT(algZero == polyZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffZero)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffOne)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(2*dof, CoeffJ)) == CoeffZero);
    }
  }
  /// properties of the unit (one) polynomial and compatibility with algebra
  void testOnePolynomial(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      Polynomial algOne = alg.one();
      Polynomial polyOne = Polynomial::One(2*dof);
      CPPUNIT_ASSERT(alg.hasElt(algOne));
      CPPUNIT_ASSERT(alg.hasElt(polyOne));
      CPPUNIT_ASSERT(algOne.getNumVars() == 2*alg.getDof());
      CPPUNIT_ASSERT(algOne.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(!algOne.isZero());
      CPPUNIT_ASSERT(algOne.isConstant());
      CPPUNIT_ASSERT(algOne.isHomogeneous());
      CPPUNIT_ASSERT(algOne.isHomogeneous(0));
      CPPUNIT_ASSERT(algOne == polyOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffZero)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffOne)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(2*dof, CoeffJ)) == CoeffOne);
    }
  }
  /// test basic coordinate monomials (see above)
  void testCoordinateMonomialQ(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial qI = alg.q(i);
        Polynomial xQ = alg.coordinateMonomial(alg.iQ(i));
        CPPUNIT_ASSERT(qI.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(xQ.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(!qI.isZero());
        CPPUNIT_ASSERT(!xQ.isZero());
        CPPUNIT_ASSERT(!qI.isConstant());
        CPPUNIT_ASSERT(!xQ.isConstant());
        CPPUNIT_ASSERT(qI.isHomogeneous());
        CPPUNIT_ASSERT(xQ.isHomogeneous());
        CPPUNIT_ASSERT(qI.isHomogeneous(1));
        CPPUNIT_ASSERT(xQ.isHomogeneous(1));
        CPPUNIT_ASSERT(qI(vec) == Coeff(i));
        CPPUNIT_ASSERT(xQ(vec) == Coeff(i));
      }
    }
  }
  /// basic momentum monomials (see above)
  void testCoordinateMonomialP(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial pI = alg.p(i);
        Polynomial xP = alg.coordinateMonomial(alg.iP(i));
        CPPUNIT_ASSERT(pI.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(xP.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(!pI.isZero());
        CPPUNIT_ASSERT(!xP.isZero());
        CPPUNIT_ASSERT(!pI.isConstant());
        CPPUNIT_ASSERT(!xP.isConstant());
        CPPUNIT_ASSERT(pI.isHomogeneous());
        CPPUNIT_ASSERT(xP.isHomogeneous());
        CPPUNIT_ASSERT(pI.isHomogeneous(1));
        CPPUNIT_ASSERT(xP.isHomogeneous(1));
        CPPUNIT_ASSERT(pI(vec) == Coeff(i));
        CPPUNIT_ASSERT(xP(vec) == Coeff(i));
      }
    }
  }
  /// powers of coordinate monomials (see above)
  void testCoordinateMonomialQPower(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial qI = alg.q(i);
        Polynomial xQ = alg.coordinateMonomial(alg.iQ(i));
        for (Power p = 0; p < CLAMaxPower; ++p) {
          Polynomial qIP = alg.q(i, p);
          Polynomial xQP = alg.coordinateMonomial(alg.iQ(i), p);
          CPPUNIT_ASSERT(qIP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(xQP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(!qIP.isZero());
          CPPUNIT_ASSERT(!xQP.isZero());
          if (p == 0) {
            CPPUNIT_ASSERT(qIP.isConstant());
            CPPUNIT_ASSERT(xQP.isConstant());
          }
          else {
            CPPUNIT_ASSERT(!qIP.isConstant());
            CPPUNIT_ASSERT(!xQP.isConstant());
          }
          CPPUNIT_ASSERT(qIP.isHomogeneous());
          CPPUNIT_ASSERT(xQP.isHomogeneous());
          CPPUNIT_ASSERT(qIP.isHomogeneous(p));
          CPPUNIT_ASSERT(xQP.isHomogeneous(p));
          CPPUNIT_ASSERT(qIP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(xQP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(qIP == qI.pow(p));
          CPPUNIT_ASSERT(xQP == xQ.pow(p));
        }
      }
    }
  }
  /// powers of momentum monomials (see above)
  void testCoordinateMonomialPPower(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial pI = alg.p(i);
        Polynomial xP = alg.coordinateMonomial(alg.iP(i));
        for (Power p = 0; p < CLAMaxPower; ++p) {
          Polynomial pIP = alg.p(i, p);
          Polynomial xPP = alg.coordinateMonomial(alg.iP(i), p);
          CPPUNIT_ASSERT(pIP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(xPP.getNumVars() == alg.getNumVars());
          CPPUNIT_ASSERT(!pIP.isZero());
          CPPUNIT_ASSERT(!xPP.isZero());
          if (p == 0) {
            CPPUNIT_ASSERT(pIP.isConstant());
            CPPUNIT_ASSERT(xPP.isConstant());
          }
          else {
            CPPUNIT_ASSERT(!pIP.isConstant());
            CPPUNIT_ASSERT(!xPP.isConstant());
          }
          CPPUNIT_ASSERT(pIP.isHomogeneous());
          CPPUNIT_ASSERT(xPP.isHomogeneous());
          CPPUNIT_ASSERT(pIP.isHomogeneous(p));
          CPPUNIT_ASSERT(xPP.isHomogeneous(p));
          CPPUNIT_ASSERT(pIP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(xPP(vec) == coeffPow(Coeff(i), p));
          CPPUNIT_ASSERT(pIP == pI.pow(p));
          CPPUNIT_ASSERT(xPP == xP.pow(p));
        }
      }
    }
  }
  /// conflicts between coordinate and momentum indices in polynomial embedding
  void testNoIndiciesAreInConflict(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      for (Index d = 0; d < dof; ++d)
        for (Index e = 0; e < dof; ++e) {
          Index iQ = alg.iQ(d);
          Index iP = alg.iP(e);
          CPPUNIT_ASSERT(0 <= iQ && iQ < alg.getNumVars());
          CPPUNIT_ASSERT(0 <= iP && iP < alg.getNumVars());
          CPPUNIT_ASSERT(iQ != iP);
          Polynomial qD = alg.q(d);
          Polynomial pE = alg.p(e);
          CPPUNIT_ASSERT(qD != pE);
        }
    }
  }
  /// ensure that all indices are used in polynomial embedding
  void testAllIndiciesAreUsed(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      ClassicalLieAlgebra alg(dof);
      std::set< Index > available;
      for (Index i = 0; i < alg.getNumVars(); ++i)
        available.insert(i);
      CPPUNIT_ASSERT(available.size() == alg.getNumVars());
      for (Index d = 0; d < dof; ++d) {
        Index iQ = alg.iQ(d);
        CPPUNIT_ASSERT(0 <= iQ && iQ < alg.getNumVars());
        CPPUNIT_ASSERT(available.erase(iQ) == 1);
      }
      for (Index e = 0; e < dof; ++e) {
        Index iP = alg.iP(e);
        CPPUNIT_ASSERT(0 <= iP && iP < alg.getNumVars());
        CPPUNIT_ASSERT(available.erase(iP) == 1);
      }
      CPPUNIT_ASSERT(available.size() == 0);
    }
  }
  /// @brief test the grade structure of the zero polynomial.
  ///
  /// mathematically, the zero polynomial belongs to all grades and,
  /// as such, it will answer true when queried whether it is isograde
  /// of grade g for all g. however, asking for its grade returns
  /// zero.
  void testGradeZero(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      const ClassicalLieAlgebra alg(dof);
      const Polynomial zero(alg.zero());
      CPPUNIT_ASSERT(alg.isIsoGrade(zero));
      CPPUNIT_ASSERT(alg.grade(zero) == 0);
      for (Power p = 0; p < 10; ++p) {
        CPPUNIT_ASSERT(alg.isIsoGrade(zero, p));
        CPPUNIT_ASSERT(alg.isoGradePart(zero, p) == zero);
        for (Power q = 0; q < 10; ++q) {
          if (p < q) {
            CPPUNIT_ASSERT(alg.isoGradePart(zero, p, q) == zero);
          }
        }
      }
    }
  }
  /// test the grade structure of the one polynomial.
  void testGradeOne(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      const ClassicalLieAlgebra alg(dof);
      const Polynomial one(alg.one());
      CPPUNIT_ASSERT(alg.isIsoGrade(one));
      CPPUNIT_ASSERT(alg.grade(one) == 0);
      for (Power p = 0; p < 10; ++p) {
        if (p == 0) {
          CPPUNIT_ASSERT(alg.isIsoGrade(one, p));
          CPPUNIT_ASSERT(alg.isoGradePart(one, p) == one);
        }
        else {
          CPPUNIT_ASSERT(!alg.isIsoGrade(one, p));
          CPPUNIT_ASSERT(alg.isoGradePart(one, p).isZero());
        }
        for (Power q = 0; q < 10; ++q) {
          if (p < q) {
            if ((p <= 0) && (q > 0)) {
              CPPUNIT_ASSERT(alg.isoGradePart(one, p, q) == one);
            }
            else {
              CPPUNIT_ASSERT(alg.isoGradePart(one, p, q).isZero());
            }
          }
        }
      }
    }
  }
  /// test the graded structure of single coordinate monomials
  void testGradeCoordinates(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      const ClassicalLieAlgebra alg(dof);
      for (Index i = 0; i < alg.getNumVars(); ++i) {
        const Polynomial xI(alg.coordinateMonomial(i));
        CPPUNIT_ASSERT(alg.isIsoGrade(xI));
        CPPUNIT_ASSERT(alg.grade(xI) == 1);
        for (Power p = 0; p < 10; ++p) {
          if (p == 1) {
            CPPUNIT_ASSERT(alg.isIsoGrade(xI, p));
            CPPUNIT_ASSERT(alg.isoGradePart(xI, p) == xI);
          }
          else {
            CPPUNIT_ASSERT(!alg.isIsoGrade(xI, p));
            CPPUNIT_ASSERT(alg.isoGradePart(xI, p).isZero());
          }
          for (Power q = 0; q < 10; ++q) {
            if (p < q) {
              if ((p <= 1) && (q > 1)) {
                CPPUNIT_ASSERT(alg.isoGradePart(xI, p, q) == xI);
              }
              else {
                CPPUNIT_ASSERT(alg.isoGradePart(xI, p, q).isZero());
              }
            }
          }
        }
      }
    }
  }
  /// test the graded structure of pure power monomials
  void testGradePurePowers(void) {
    for (Index dof = 1; dof < CLAMaxDof; ++dof) {
      const ClassicalLieAlgebra alg(dof);
      for (Index i = 0; i < alg.getNumVars(); ++i) {
        for (Index j = 0; j < 10; ++j) {
          const Polynomial xI(alg.coordinateMonomial(i, j));
          CPPUNIT_ASSERT(alg.isIsoGrade(xI));
          CPPUNIT_ASSERT(alg.grade(xI) == j);
          for (Power p = 0; p < 10; ++p) {
            if (p == j) {
              CPPUNIT_ASSERT(alg.isIsoGrade(xI, p));
              CPPUNIT_ASSERT(alg.isoGradePart(xI, p) == xI);
            }
            else {
              CPPUNIT_ASSERT(!alg.isIsoGrade(xI, p));
              CPPUNIT_ASSERT(alg.isoGradePart(xI, p).isZero());
            }
            for (Power q = 0; q < 10; ++q) {
              if (p < q) {
                if ((p <= j) && (q > j)) {
                  CPPUNIT_ASSERT(alg.isoGradePart(xI, p, q) == xI);
                }
                else {
                  CPPUNIT_ASSERT(alg.isoGradePart(xI, p, q).isZero());
                }
              }
            }
          }
        }
      }
    }
  }
  /// @brief test isograde part extraction on random examples
  ///
  /// ensure that a polynomial is the sum of its isograde parts
  /// (random examples)
  void testSumOfIsoGradesIsPolynomial(void) {
    for (size_t c = 0; c < 100; ++c) {
      const Index dof(1 + randomIndex(25));
      const ClassicalLieAlgebra alg(dof);
      const Polynomial p(randomPoly(alg.getNumVars()));
      Polynomial sumIsoGrades(alg.zero());
      for (Power g = 0; g <= alg.grade(p); ++g) {
        sumIsoGrades += alg.isoGradePart(p, g);
      }
      CPPUNIT_ASSERT(sumIsoGrades == p);
    }
  }
  /// @brief test isograde extraction on a particular example.
  ///
  /// ensure that a specific example polynomial is the sum of its
  /// isograde parts (avoid potential degeneracy in the random
  /// examples)
  void testSumOfIsoGradesIsPolynomialExample(void) {
    const Index dof(2);
    const ClassicalLieAlgebra alg(dof);
    const Polynomial q0(alg.q(0));
    const Polynomial p0(alg.p(0));
    const Polynomial q1(alg.q(1));
    const Polynomial p1(alg.p(1));
    const Polynomial one(alg.one());
    Polynomial p(alg.zero());
    p += Coeff(-5) * one + Coeff(2) * q0;
    p += Coeff(-3) * q0 * p0 * p1;
    p *= one + Coeff(4) * q1;
    CPPUNIT_ASSERT(alg.grade(p) == 4);
    CPPUNIT_ASSERT(alg.isoGradePart(p, 0) == Coeff(-5) * one);
    CPPUNIT_ASSERT(alg.isoGradePart(p, 1) == Coeff(2)* q0 + Coeff(-20) * q1);
    CPPUNIT_ASSERT(alg.isoGradePart(p, 2) == Coeff(8)* q0 * q1);
    CPPUNIT_ASSERT(alg.isoGradePart(p, 3) == Coeff(-3) * q0 * p0 * p1);
    CPPUNIT_ASSERT(alg.isoGradePart(p, 4) == Coeff(-12) * q0 * p0 * q1 * p1);
    for (Power po = 5; po < 10; ++po) {
      CPPUNIT_ASSERT(alg.isoGradePart(p, po) == alg.zero());
    }
  }

 private:
  
}; //ClassicalLieAlgebraTest

#endif //CLASSICAL_LIE_ALGEBRA_TEST_H
