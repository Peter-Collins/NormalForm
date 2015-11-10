//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: CoordinateChange implementation
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers
#include "CoordinateChange.h"
#include "Utility.h"

CoordinateChange::~CoordinateChange(void) {
}

CoordinateChange::CoordinateChange(const LieAlgebraBase& algebra,
                                   const std::vector< Polynomial >& wI,
                                   const bool isReverse)
  : algebra_(algebra), wI_(wI), isReverse_(isReverse), currentGrade_(0) {
  checkInnerGeneratorTerms_();
}

void CoordinateChange::checkInnerGeneratorTerms_(void) const {
  std::vector< Polynomial >::const_iterator wTerm;
  Power count;
  for (wTerm = wI_.begin(), count = 0;
       wTerm != wI_.end();
       ++wTerm, ++count) {
    //the generator should begin with zero quadratic term
    if (!algebra_.isIsoGrade(*wTerm, count + 2)) {
      throw CoordinateChangeGradeError();
    }
    if (count == 0) {
      if (!(algebra_.isoGradePart(*wTerm, count + 2).isZero())) {
        throw CoordinateChangeGradeError();
      }
    }
  }
}

void CoordinateChange::computeNextTerm(const Polynomial& nextInputTerm,
                                       Polynomial& result) {
  if (isReverse_) {
    reverseCoordinateChange_(nextInputTerm, result);
  }
  else {
    forwardCoordinateChange_(nextInputTerm, result);
  }
}

void CoordinateChange::checkGrades_(const Power i,
                                    const Power j,
                                    const Polynomial& fTerm,
                                    const Polynomial& wTerm) const {
  assert(i + j == (currentGrade_ - 1));
  const Power gradeF(algebra_.grade(fTerm));
  const Power gradeW(algebra_.grade(wTerm));
  assert((gradeF + gradeW - 2) <= (i + j + 1));
  if ((gradeF + gradeW - 2) != (i + j + 1))
    assert((gradeF == 0) || (gradeW == 0));
}

//forward: f:diag->scalar becomes f':norm->scalar (direct)
void
CoordinateChange::forwardCoordinateChange_(const Polynomial& fTerm,
                                           Polynomial& resultFTerm) {
  //prepare the triangle row
  currentGrade_ += 1;
  const Power i = currentGrade_ - 1;
  assert(currentGrade_ >= 1);
  assert(algebra_.isIsoGrade(fTerm, currentGrade_));
  assert(i >= 0);
  assert(wI_.size() >= i);
  assert(fIJ_.size() == ((i * (i + 1)) / 2));
  //prepare the original function term, removing factorial divisor
  //row i is grade i + 1 associated with factorial (i - 1)
  //const Coeff preFactor((i == 0) ? 1 : factorial(i - 1));
  const Coeff preFactor(factorial(i));
  fIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(i, 0),
                                                       fTerm * preFactor));
  //compute the row
  const Polynomial innerResultFTerm(forwardTriangle_());
  //return the polynomial terms, correctly weighted by the preFactor
  const Coeff divFactor = CoeffOne / preFactor;
  resultFTerm = innerResultFTerm * divFactor;
  assert(algebra_.isIsoGrade(resultFTerm, currentGrade_));
}

const Polynomial&
CoordinateChange::forwardTriangle_() {
  std::map< TriangleAddress, Polynomial >::const_iterator iterFIJ;
  const size_t wSize(wI_.size());
  for (Power j = 1; j <= currentGrade_ - 1; ++j) {
    for (Power i = 0; i <= currentGrade_ - 1 - j; ++i) {
      iterFIJ = fIJ_.find(TriangleAddress(i, j));
      if (iterFIJ == fIJ_.end()) {
        iterFIJ = fIJ_.find(TriangleAddress(i + 1, j - 1));
        assert(iterFIJ != fIJ_.end());
        Polynomial temp(iterFIJ->second);
        for (Power k = 0; k <= i; ++k) {
          assert((i + 1 - k) < wSize);
          const Polynomial& wTerm(wI_.at(i + 1 - k));
          iterFIJ = fIJ_.find(TriangleAddress(k, j - 1));
          assert(iterFIJ != fIJ_.end());
          const Polynomial& inner(iterFIJ->second);
          checkGrades_(i, j, inner, wTerm);
          temp += Coeff(binomial(i, k)) * algebra_.lieBracket(inner, wTerm);
        }
        fIJ_.insert(std::map< TriangleAddress, Polynomial >::value_type(TriangleAddress(i, j), temp));
      }
    }
  }
  iterFIJ = fIJ_.find(TriangleAddress(0, currentGrade_ - 1));
  assert(iterFIJ != fIJ_.end());
  return (iterFIJ->second);
}

//reverse: f:diag->scalar becomes f':norm->scalar (inverse)
//
// note that the fundamental theory of lie series does not apply and
// one therefore _cannot_ use the negative of the generator as the
// generator of the inverse!!!
//
void
CoordinateChange::reverseCoordinateChange_(const Polynomial& fTerm,
                                           Polynomial& resultFTerm) {
  //prepare the triangle row
  currentGrade_ += 1;
  const Power i = currentGrade_ - 1;
  assert(currentGrade_ >= 1);
  assert(algebra_.isIsoGrade(fTerm, currentGrade_));
  assert(i >= 0);
  assert(wI_.size() >= i);
  assert(fIJ_.size() == ((i * (i + 1)) / 2));
  //prepare the original function term, removing factorial divisor
  //the divisor always matches the power on epsilon
  const Coeff preFactor(factorial(i));
  fIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(0, i),
                                                       fTerm * preFactor));
  //compute the row
  const Polynomial innerResultFTerm(reverseTriangle_());
  //return the polynomial terms, correctly weighted by the preFactor
  const Coeff divFactor = CoeffOne / preFactor;
  resultFTerm = innerResultFTerm * divFactor;
  assert(algebra_.isIsoGrade(resultFTerm, currentGrade_));
}

const Polynomial&
CoordinateChange::reverseTriangle_() {
  std::map< TriangleAddress, Polynomial >::const_iterator iterFIJ;
  const size_t wSize(wI_.size());
  for (Power i = 1; i <= currentGrade_ - 1; ++i) {
    for (Power j = 0; j <= currentGrade_ - 1 - i; ++j) {
      iterFIJ = fIJ_.find(TriangleAddress(i, j));
      if (iterFIJ == fIJ_.end()) {
        iterFIJ = fIJ_.find(TriangleAddress(i - 1, j + 1));
        assert(iterFIJ != fIJ_.end());
        Polynomial temp(iterFIJ->second);
        for (Power k = 0; k <= i - 1; ++k) { //j - 1?
          assert((k + 1) < wSize);
          const Polynomial& wTerm(wI_.at(k + 1));
          iterFIJ = fIJ_.find(TriangleAddress(i - k - 1, j));
          assert(iterFIJ != fIJ_.end());
          const Polynomial& inner(iterFIJ->second);
          checkGrades_(i, j, inner, wTerm);
          temp -= Coeff(binomial(i - 1, k)) * algebra_.lieBracket(inner, wTerm);
        }
        fIJ_.insert(std::map< TriangleAddress, Polynomial >::value_type(TriangleAddress(i, j), temp));
      }
    }
  }
  iterFIJ = fIJ_.find(TriangleAddress(currentGrade_ - 1, 0));
  assert(iterFIJ != fIJ_.end());
  return (iterFIJ->second);
}
