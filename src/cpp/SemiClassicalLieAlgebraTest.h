#ifndef SEMI_CLASSICAL_LIE_ALGEBRA_TEST_H
#define SEMI_CLASSICAL_LIE_ALGEBRA_TEST_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
//----------------------------------------------------------------------

//system headers
#include <set>

//library headers
#include <cppunit/extensions/HelperMacros.h>

//project headers
#include "Types.h"
#include "LieAlgebraBase.h"
#include "SemiClassicalLieAlgebra.h"
#include "LieAlgebraTestBase.h"

//we must use =, not () to initialize static data members!  never do
//this in a class declaration!
static const Index SCLAMaxDof = 50;
static const Power SCLAMaxPower = 8;

class SemiClassicalLieAlgebraTest : public CppUnit::TestFixture {

 public:

  //the following create static CppUnit::TestSuite *suite()
  CPPUNIT_TEST_SUITE(SemiClassicalLieAlgebraTest);

  //constructor from number of variables
  CPPUNIT_TEST_EXCEPTION(testConstructorGivenZeroVars, LieAlgebraIndexError);
  CPPUNIT_TEST(testConstructorGivenNumVars);
  
  //big four
  CPPUNIT_TEST(testDefaultConstructor);
  CPPUNIT_TEST(testCopyConstructor);
  //CPPUNIT_TEST(testAssign); //forbidden

  //basic object creation
  CPPUNIT_TEST(testZeroPolynomial);
  CPPUNIT_TEST(testOnePolynomial);
  CPPUNIT_TEST(testCoordinateMonomialQ);
  CPPUNIT_TEST(testCoordinateMonomialP);
  CPPUNIT_TEST(testCoordinateMonomialQPower);
  CPPUNIT_TEST(testCoordinateMonomialPPower);
  CPPUNIT_TEST(testCoordinateMonomialHBar);
  CPPUNIT_TEST(testCoordinateMonomialHBarPower);

  //test the q, p, h-bar encoding
  CPPUNIT_TEST(testNoIndiciesAreInConflict);
  CPPUNIT_TEST(testAllIndiciesAreUsed);

  CPPUNIT_TEST_SUITE_END();

protected:

  void testConstructorGivenZeroVars(void) {
    //this should throw an index error
    SemiClassicalLieAlgebra alg(0);
  }

  void testConstructorGivenNumVars(void) {
    //construct algs with given number of variables
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof+1);
      Polynomial p(2*dof+1);
      CPPUNIT_ASSERT(alg.hasElt(p));
      CPPUNIT_ASSERT(alg.grade(p) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(p));
      CPPUNIT_ASSERT(alg.isoGradePart(p, 0) == p);
    }
  }

  void testCopyConstructor(void) {
    //construct a copy of a alg
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      SemiClassicalLieAlgebra algCopy(alg);
      CPPUNIT_ASSERT(alg.getDof() == dof);
      CPPUNIT_ASSERT(alg.getNumVars() == 2*dof+1);
      CPPUNIT_ASSERT(algCopy.getDof() == dof);
      CPPUNIT_ASSERT(algCopy.getNumVars() == 2*dof+1);
    }
  }

  void testDefaultConstructor(void) {
    //test the default constructor gives a 1-dof alg
    SemiClassicalLieAlgebra alg;
    CPPUNIT_ASSERT(alg.getDof() == 1);
    CPPUNIT_ASSERT(alg.getNumVars() == 3);
    Polynomial p(3);
    CPPUNIT_ASSERT(alg.hasElt(p));
    CPPUNIT_ASSERT(alg.grade(p) == 0);
    CPPUNIT_ASSERT(alg.isIsoGrade(p));
    CPPUNIT_ASSERT(alg.isoGradePart(p, 0) == p);
  }

  void testZeroPolynomial(void) {
    //properties of zero polynomial and compatibility with poly zero
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      Polynomial algZero = alg.zero();
      Polynomial polyZero = Polynomial::Zero(2*dof+1);
      CPPUNIT_ASSERT(alg.hasElt(algZero));
      CPPUNIT_ASSERT(alg.hasElt(polyZero));
      CPPUNIT_ASSERT(algZero.getNumVars() == 2*alg.getDof()+1);
      CPPUNIT_ASSERT(algZero.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(algZero.isZero());
      CPPUNIT_ASSERT(algZero.isConstant());
      CPPUNIT_ASSERT(algZero.isHomogeneous());
      CPPUNIT_ASSERT(algZero.isHomogeneous(0));
      CPPUNIT_ASSERT(alg.grade(algZero) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(algZero));
      CPPUNIT_ASSERT(alg.isoGradePart(algZero, 0) == algZero);
      CPPUNIT_ASSERT(algZero == polyZero);
      CPPUNIT_ASSERT(algZero(makeVector(1+2*dof, CoeffZero)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(1+2*dof, CoeffOne)) == CoeffZero);
      CPPUNIT_ASSERT(algZero(makeVector(1+2*dof, CoeffJ)) == CoeffZero);
    }
  }

  void testOnePolynomial(void) {
    //properties of one polynomial and compatibility with poly one
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      Polynomial algOne = alg.one();
      Polynomial polyOne = Polynomial::One(2*dof+1);
      CPPUNIT_ASSERT(alg.hasElt(algOne));
      CPPUNIT_ASSERT(alg.hasElt(polyOne));
      CPPUNIT_ASSERT(algOne.getNumVars() == 2*alg.getDof()+1);
      CPPUNIT_ASSERT(algOne.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(!algOne.isZero());
      CPPUNIT_ASSERT(algOne.isConstant());
      CPPUNIT_ASSERT(algOne.isHomogeneous());
      CPPUNIT_ASSERT(algOne.isHomogeneous(0));
      CPPUNIT_ASSERT(alg.grade(algOne) == 0);
      CPPUNIT_ASSERT(alg.isIsoGrade(algOne));
      CPPUNIT_ASSERT(alg.isoGradePart(algOne, 0) == algOne);
      CPPUNIT_ASSERT(algOne == polyOne);
      CPPUNIT_ASSERT(algOne(makeVector(1+2*dof, CoeffZero)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(1+2*dof, CoeffOne)) == CoeffOne);
      CPPUNIT_ASSERT(algOne(makeVector(1+2*dof, CoeffJ)) == CoeffOne);
    }
  }

  void testCoordinateMonomialQ(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
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

  void testCoordinateMonomialP(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
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

  void testCoordinateMonomialHBar(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      vec[alg.iHBar()] = Coeff(-1);
      Polynomial hBar = alg.hBar();
      CPPUNIT_ASSERT(alg.hasElt(hBar));
      CPPUNIT_ASSERT(hBar.getNumVars() == alg.getNumVars());
      CPPUNIT_ASSERT(!hBar.isZero());
      CPPUNIT_ASSERT(!hBar.isConstant());
      CPPUNIT_ASSERT(hBar.isHomogeneous());
      CPPUNIT_ASSERT(hBar.isHomogeneous(1));
      CPPUNIT_ASSERT(hBar(vec) == Coeff(-1));
      CPPUNIT_ASSERT(alg.grade(hBar) == 2);
      CPPUNIT_ASSERT(alg.isIsoGrade(hBar));
      CPPUNIT_ASSERT(alg.isoGradePart(hBar, 2) == hBar);
    }
  }

  void testCoordinateMonomialQPower(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial qI = alg.q(i);
        Polynomial xQ = alg.coordinateMonomial(alg.iQ(i));
        for (Power p = 0; p < SCLAMaxPower; ++p) {
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
  void testCoordinateMonomialPPower(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      for (Index i = 0; i < dof; ++i) {
        Polynomial pI = alg.p(i);
        Polynomial xP = alg.coordinateMonomial(alg.iP(i));
        for (Power p = 0; p < SCLAMaxPower; ++p) {
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

  void testCoordinateMonomialHBarPower(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      std::vector< Coeff > vec;
      for (Index j = 0; j < 2*dof+1; ++j)
        vec.push_back(CoeffZero);
      for (Index d = 0; d < dof; ++d) {
        vec[alg.iQ(d)] = Coeff(d);
        vec[alg.iP(d)] = Coeff(d);
      }
      vec[alg.iHBar()] = Coeff(-1);
      Polynomial hBar = alg.hBar();
      for (Power p = 0; p < SCLAMaxPower; ++p) {
        Polynomial hBarP = alg.hBar(p);
        CPPUNIT_ASSERT(alg.hasElt(hBarP));
        CPPUNIT_ASSERT(hBarP.getNumVars() == alg.getNumVars());
        CPPUNIT_ASSERT(hBarP == hBar.pow(p));
        CPPUNIT_ASSERT(!hBarP.isZero());
        if (p == 0)
          CPPUNIT_ASSERT(hBarP.isConstant());
        else
          CPPUNIT_ASSERT(!hBarP.isConstant());
        CPPUNIT_ASSERT(hBarP.isHomogeneous());
        CPPUNIT_ASSERT(hBarP.isHomogeneous(p));
        CPPUNIT_ASSERT(hBarP(vec) == coeffPow(Coeff(-1), p));
      }
    }
  }

  void testNoIndiciesAreInConflict(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
      const Index iHBar = alg.iHBar();
      const Polynomial hBar = alg.hBar();
      CPPUNIT_ASSERT(hBar == alg.coordinateMonomial(iHBar));
      CPPUNIT_ASSERT(0 <= iHBar && iHBar < alg.getNumVars());
      for (Index d = 0; d < dof; ++d)
        for (Index e = 0; e < dof; ++e) {
          Index iQ = alg.iQ(d);
          Index iP = alg.iP(e);
          CPPUNIT_ASSERT(0 <= iQ && iQ < alg.getNumVars());
          CPPUNIT_ASSERT(0 <= iP && iP < alg.getNumVars());
          CPPUNIT_ASSERT(iQ != iP);
          CPPUNIT_ASSERT(iQ != iHBar);
          CPPUNIT_ASSERT(iP != iHBar);
          Polynomial qD = alg.q(d);
          Polynomial pE = alg.p(e);
          CPPUNIT_ASSERT(qD != pE);
          CPPUNIT_ASSERT(qD != hBar);
          CPPUNIT_ASSERT(pE != hBar);
        }
    }
  }

  void testAllIndiciesAreUsed(void) {
    for (Index dof = 1; dof < SCLAMaxDof; ++dof) {
      SemiClassicalLieAlgebra alg(dof);
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
      Index iHBar = alg.iHBar();
      CPPUNIT_ASSERT(0 <= iHBar && iHBar < alg.getNumVars());
      CPPUNIT_ASSERT(available.erase(iHBar) == 1);
      CPPUNIT_ASSERT(available.size() == 0);
    }
  }

 private:
  
}; //SemiClassicalLieAlgebraTest

#endif //SEMI_CLASSICAL_LIE_ALGEBRA_TEST_H
