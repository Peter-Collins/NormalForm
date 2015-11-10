//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: LieAlgebraBase implementation
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers
#include "LieAlgebraBase.h"

LieAlgebraBase::LieAlgebraBase(void)
  : dof_(1) {
}

LieAlgebraBase::LieAlgebraBase(const LieAlgebraBase& that)
  : dof_(that.dof_) {
}

LieAlgebraBase::~LieAlgebraBase(void) {
  //nothing
}

LieAlgebraBase::LieAlgebraBase(const Index dof)
  : dof_(dof) {
  if (!(dof > 0))
    throw LieAlgebraIndexError();
}

bool LieAlgebraBase::isDiagonal(const Powers& powers) const {
  if (!(powers.getNumVars() == getNumVars()))
    throw LieAlgebraSizeMismatchError();
  for (Index d = 0; d < dof_; ++d) {
    if (powers[iQ(d)] != powers[iP(d)]) {
      return false;
    }
  }
  return true;
}

bool LieAlgebraBase::isDiagonal(const Polynomial& pol) const {
  if (!hasElt(pol))
    throw LieAlgebraSizeMismatchError();
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms.begin(); pc != terms.end(); ++pc) {
    if (!isDiagonal(pc->first)) {
      return false;
    }
  }
  return true;
}

void
LieAlgebraBase::diagonalAndNonDiagonal(const Polynomial& pol,
                                       std::pair< Polynomial, Polynomial >& dAndN) const {
  if (!hasElt(pol))
    throw LieAlgebraSizeMismatchError();
  dAndN.first = zero();
  dAndN.second = zero();
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms.begin(); pc != terms.end(); ++pc) {
    if (isDiagonal(pc->first)) {
      (dAndN.first).addMonomial(pc->first, pc->second);
    }
    else {
      (dAndN.second).addMonomial(pc->first, pc->second);
    }
  }
}

//
// implement the graded algebra structure in terms of the single pure
// virtual function, gradeOfPowers_, rather than trying to do it in
// terms of the three pure virtuals in the old implementation.
//

Power LieAlgebraBase::grade(const Polynomial& pol) const {
  if (!hasElt(pol))
    throw LieAlgebraSizeMismatchError();
  Power maxGra(0);
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  Power gra;
  for (pc = terms.begin(); pc != terms.end(); ++pc) {
    gra = gradeOfPowers_(pc->first);
    if (gra > maxGra) {
      maxGra = gra;
    }
  }
  return maxGra;
}

bool LieAlgebraBase::isIsoGrade(const Polynomial& pol) const {
  if (!hasElt(pol))
    throw LieAlgebraSizeMismatchError();
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc = terms.begin();
  if (pc == terms.end())
    return true;
  Power gra(gradeOfPowers_(pc->first));
  for (; pc != terms.end(); ++pc) {
    if (gra != gradeOfPowers_(pc->first))
      return false;
  }
  return true;
}

bool LieAlgebraBase::isIsoGrade(const Polynomial& pol, const Power gra) const {
  if (!hasElt(pol))
    throw LieAlgebraSizeMismatchError();
  if (!isNonnegativePower(gra))
    throw LieAlgebraPowerError();
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  for (PowersToCoeffMap::const_iterator pc = terms.begin();
       pc != terms.end();
       ++pc)
    if (gra != gradeOfPowers_(pc->first))
      return false;
  return true;
}

Polynomial LieAlgebraBase::isoGradePart(const Polynomial& pol,
                                        const Power grade,
                                        const Power lessThan) const {
  if (grade >= lessThan) {
    throw LieAlgebraPowerError();
  }
  Polynomial res = zero();
  Power gra;
  const PowersToCoeffMap& terms = pol.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator pc;
  for (pc = terms.begin(); pc != terms.end(); ++pc) {
    gra = gradeOfPowers_(pc->first);
    if ((grade <= gra) && (gra < lessThan))
      res.setMonomial(pc->first, pc->second);
  }
  return res;
}
