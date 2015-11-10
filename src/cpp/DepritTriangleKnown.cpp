//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangleKnown implementation
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers
#include "DepritTriangleKnown.h"
#include "Utility.h"
#include "Polynomial.h"
#include "LieAlgebraBase.h"

DepritTriangleKnown::~DepritTriangleKnown(void) {
}

DepritTriangleKnown::DepritTriangleKnown(const LieAlgebraBase& algebra,
                                         const std::vector< Polynomial >& innerGeneratorTerms)
  : algebra_(algebra), wI_(innerGeneratorTerms), currentGrade_(1) {
  checkInnerGeneratorTerms_();
}

void
DepritTriangleKnown::checkInnerGeneratorTerms_(void) const {
  for (size_t i = 0; i < wI_.size() ; ++i) {
    if (!algebra_.isIsoGrade(wI_.at(i), i + 2)) {
      throw DepritTriangleKnownGeneratorError();
    }
  }
  if (wI_.size() >= 1) {
    if (wI_.at(0) != algebra_.zero()) {
      throw DepritTriangleKnownGeneratorError();
    }
  }
}

void
DepritTriangleKnown::checkGrades_(const Power i,
                                  const Power j,
                                  const Polynomial& fTerm,
                                  const Polynomial& wTerm) const {
  assert(i + j == (currentGrade_ - 2));
  const Power gradeF(algebra_.grade(fTerm));
  const Power gradeW(algebra_.grade(wTerm));
  assert((gradeF + gradeW - 2) <= (i + j + 2));
  if ((gradeF + gradeW - 2) != (i + j + 2))
    assert((gradeF == 0) || (gradeW == 0));
}

void
DepritTriangleKnown::transform_(const Polynomial& fTerm,
                                Polynomial& resultFTerm,
                                const bool isBackward) {
  //prepare Deprit triangle row
  currentGrade_ += 1;
  Power i = currentGrade_ - 2;
  assert(currentGrade_ >= 2);
  assert(i >= 0);
  assert(wI_.size() > i);

  //prepare input term, removing factorial divisor to give hI
  assert(algebra_.isIsoGrade(fTerm, currentGrade_));
  const Coeff preFactor(factorial(i));
  assert(fIJ_.size() == ((i * (i + 1)) / 2));
  fIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(i, 0),
                                                       preFactor * fTerm));
  const Polynomial innerResultFTerm(triangle_(0, i, isBackward));

  //return the polynomial terms, correctly weighted by the preFactor
  const Coeff divFactor = CoeffOne / preFactor;
  resultFTerm = innerResultFTerm * divFactor;
  assert(algebra_.isIsoGrade(resultFTerm, currentGrade_));
}

Polynomial
DepritTriangleKnown::triangle_(const Power i,
                               const Power j,
                               const bool isBackward) {
  std::map< TriangleAddress, Polynomial >::iterator iterFIJ;
  iterFIJ = fIJ_.find(TriangleAddress(i, j));
  if (iterFIJ != fIJ_.end()) {
    return (iterFIJ->second);
  }
  //term not found
  assert(j != 0);
  Polynomial temp(triangle_(i + 1, j - 1, isBackward));
  for (Power k = 0; k <= i; ++k) {
    const Polynomial inner(triangle_(i - k, j - 1, isBackward));
    checkGrades_(i, j, inner, wI_.at(k+1));
    if (isBackward) {
      temp -= Coeff(binomial(i, k)) * algebra_.lieBracket(inner, wI_.at(k+1));
    }
    else {
      temp += Coeff(binomial(i, k)) * algebra_.lieBracket(inner, wI_.at(k+1));
    }
  }
  fIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(i, j),
                                                       temp));
  iterFIJ = fIJ_.find(TriangleAddress(i, j));
  assert(iterFIJ != fIJ_.end());
  return (iterFIJ->second);
}
