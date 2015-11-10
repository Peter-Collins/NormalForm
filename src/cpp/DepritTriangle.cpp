//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangle implementation
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers
#include "DepritTriangle.h"
#include "Utility.h"
#include "Polynomial.h"
#include "LieAlgebraBase.h"

DepritTriangle::~DepritTriangle(void) {
}

DepritTriangle::DepritTriangle(const LieAlgebraBase& algebra,
                               const Polynomial& quadratic)
  : algebra_(algebra), quadratic_(quadratic), currentGrade_(1) {
  std::cerr << "Deprit triangle construction" << std::endl;
  checkQuadraticPart_();
  extractFundamentalFrequencies_();
}

void
DepritTriangle::checkQuadraticPart_(void) const {
  std::cerr << "Check quadratic part is in correct algebra" << std::endl;
  assert(algebra_.hasElt(quadratic_));
  std::cerr << "Check quadratic part has zero constant and linear terms" << std::endl;
  assert(algebra_.isoGradePart(quadratic_, 0, 2).isZero());
  std::cerr << "Check quadratic part is diagonal" << std::endl;
  assert(algebra_.isDiagonal(quadratic_));
}

void
DepritTriangle::extractFundamentalFrequencies_(void) {
  std::cerr << "Extract fundamental frequencies" << std::endl;
  Coeff freq;
  for (Index d = 0; d < algebra_.getDof(); ++d) {
    Powers powers(algebra_.getNumVars());
    powers.setPower(algebra_.iQ(d), 1);
    powers.setPower(algebra_.iP(d), 1);
    freq = quadratic_[powers];
    std::cerr << freq << std::endl;
    frequencies_.push_back(freq);
  }
  assert(frequencies_.size() == algebra_.getDof());
}

void
DepritTriangle::partitionKnownTerms_(const Polynomial& knownTerms,
                                     Polynomial& resultNormal,
                                     Polynomial& resultRemainder) const {
  std::cerr << "Partitioning the known terms..." << std::endl;
  resultNormal = algebra_.zero();
  resultRemainder = algebra_.zero();
  const PowersToCoeffMap& terms = knownTerms.getPowersAndCoeffs();
  for (PowersToCoeffMap::const_iterator pc = terms.begin();
       pc != terms.end();
       ++pc)
    if (isInNormalForm_(pc->first))
      resultNormal.setMonomial(pc->first, pc->second);
    else
      resultRemainder.setMonomial(pc->first, pc->second);
  assert(resultNormal.getNumTerms() + resultRemainder.getNumTerms() == knownTerms.getNumTerms());
}

void
DepritTriangle::computeNormalFormAndGenerator(const Polynomial& hTerm,
                                              Polynomial& resultKTerm,
                                              Polynomial& resultWTerm) {
  //prepare Deprit triangle row
  currentGrade_ += 1;
  Power i = currentGrade_ - 2;
  assert(currentGrade_ >= 2);
  assert(i >= 0);
  assert(wI_.size() == i);

  //report starting the row
  for (int m = 0; m < 72; ++m)
      std::cerr << "=";
  std::cerr << std::endl;
  std::cerr << "Deprit triangle row: " << i << " ";
  std::cerr << "Grade: " << currentGrade_ << std::endl;

  //prepare hamiltonian term, removing factorial divisor to give hI
  assert(algebra_.isIsoGrade(hTerm, currentGrade_));
  const Coeff preFactor(factorial(i));
  assert(hIJ_.size() == ((i * (i + 1)) / 2));
  //std container indexing default constructs, which prevents assign poly!
  //this next bit is a bit stupid; why create the pair (which uses copy?)
  hIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(i, 0),
                                                       preFactor * hTerm));

  //partition the known terms into normalised part and remainder
  std::cerr << "Computing deprit triangle row without unknown..." << std::endl;
  //could the following actually be a reference?
  const Polynomial knownTerms(triangleMinusUnknownTerm_(0, i));
  Polynomial normal(algebra_.zero());
  Polynomial remainder(algebra_.zero());
  partitionKnownTerms_(knownTerms, normal, remainder);
  correctRowUsingPartition_(normal, remainder);

  //solve the homological equation
  wI_.push_back(solveHomologicalEquation_(remainder));
  checkHomologicalEquation_(remainder);

  //report on some sizes
  std::cerr << "InputTerms: " << hTerm.getNumTerms() << std::endl;
  std::cerr << "KnownTerms: " << knownTerms.getNumTerms() << std::endl;
  std::cerr << "NormalForm: " << normal.getNumTerms() << std::endl;
  std::cerr << "Generator : " << remainder.getNumTerms() << std::endl;

  //return the polynomial terms, correctly weighted by the preFactor
  const Coeff divFactor = CoeffOne / preFactor;
  resultKTerm = normal*divFactor;
  resultWTerm = wI_.at(i)*divFactor;

  std::cerr << "Confirming isograde structure..." << std::endl;
  assert(algebra_.isIsoGrade(resultKTerm, currentGrade_));
  assert(algebra_.isIsoGrade(resultWTerm, currentGrade_));
}

Polynomial
DepritTriangle::solveHomologicalEquation_(const Polynomial& remainder) const {
  std::cerr << "Solving the homological equation..." << std::endl;
  Polynomial result(algebra_.zero());
  const PowersToCoeffMap& terms = remainder.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms.begin(); pc != terms.end(); ++pc) {
    //the following could be done more efficiently for sparse powers; later
    Coeff innerProd(0.0);
    const Powers& powers = (pc->first);
    for (Index d = 0; d < algebra_.getDof(); ++d)
      innerProd += Coeff(powers[algebra_.iP(d)] - powers[algebra_.iQ(d)]) * frequencies_[d];
    if (fabs(innerProd) < 1.0e-6)
      std::cerr << "WARNING! small inner product " << innerProd << std::endl;
    result.setMonomial(powers, -(pc->second)/innerProd);
  }
  return result;
}

void
DepritTriangle::checkHomologicalEquation_(const Polynomial& remainder) const {
  std::cerr << "Computing the homological equation error..." << std::endl;
  const Power i(currentGrade_ - 2);
  const Polynomial bracket(algebra_.lieBracket(quadratic_, wI_.at(i)));
  const Real err((bracket + remainder).lInfinityNorm());
  std::cerr << "Error in the homological equation: " << err << std::endl;
  assert(err < 5.0e-6);
}

void
DepritTriangle::checkGrades_(const Power i,
                             const Power j,
                             const Polynomial& hTerm,
                             const Polynomial& wTerm) const {
  assert(2 + (i + j) == currentGrade_);
  const Power gradeH(algebra_.grade(hTerm));
  const Power gradeW(algebra_.grade(wTerm));
  assert((gradeH + gradeW - 2) <= (2 + (i + j)));
  if ((gradeH + gradeW - 2) != (2 + (i + j)))
    assert((gradeH == 0) || (gradeW == 0));
}

Polynomial
DepritTriangle::triangleMinusUnknownTerm_(const Power i, const Power j) {
  std::cerr << "Triangle " << i << " " << j << std::endl;
  std::map< TriangleAddress, Polynomial >::const_iterator iterHIJ;
  iterHIJ = hIJ_.find(TriangleAddress(i, j));
  if (iterHIJ != hIJ_.end()) {
    return (iterHIJ->second);
  }
  //term not found
  assert(j != 0);
  Polynomial temp(triangleMinusUnknownTerm_(i + 1, j - 1));
  for (Power k = 0; k <= i; ++k) {
    if ((k + 1) == (i + j)) {
      //unknown term
      assert((i - k) == 0);
      assert((j - 1) == 0);
    }
    else {
      //should the following not be a reference?
      const Polynomial& inner(triangleMinusUnknownTerm_(i - k, j - 1));
      checkGrades_(i, j, inner, wI_.at(k+1));
      temp += Coeff(binomial(i, k)) * algebra_.lieBracket(inner, wI_.at(k+1));
    }
  }
  //the above must be corrected after solution of the homological equation
  hIJ_.insert(std::pair< TriangleAddress, Polynomial >(TriangleAddress(i, j),
                                                       temp));
  iterHIJ = hIJ_.find(TriangleAddress(i, j));
  assert(iterHIJ != hIJ_.end());
  return (iterHIJ->second);
}

void
DepritTriangle::correctRowUsingPartition_(const Polynomial& normal,
                                          const Polynomial& remainder) {
  std::cerr << "Correcting the Deprit triangle row using remainder." << std::endl;
  const Power i(currentGrade_ - 2);
  for (Power j = 1; j <= i; ++j)
    hIJ_[TriangleAddress(i - j, j)] -= remainder;
  const Real err((hIJ_[TriangleAddress(0, i)] - normal).lInfinityNorm());
  std::cerr << "Error in predicted normal form term: " << err << std::endl;
  assert(err < 1.0e-6);
}
