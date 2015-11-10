//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: SemiClassicalLieAlgebra implementation
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>

//library headers

//project headers
#include "Utility.h"
#include "Polynomial.h"
#include "SemiClassicalLieAlgebra.h"

//big four

SemiClassicalLieAlgebra::SemiClassicalLieAlgebra(void)
  : LieAlgebraBase() {
}

SemiClassicalLieAlgebra::SemiClassicalLieAlgebra(const SemiClassicalLieAlgebra& that)
  : LieAlgebraBase(that) {
}

SemiClassicalLieAlgebra::~SemiClassicalLieAlgebra(void) {
}

//constructors

SemiClassicalLieAlgebra::SemiClassicalLieAlgebra(const Index dof)
  : LieAlgebraBase(dof) {
}

Polynomial SemiClassicalLieAlgebra::lieBracket_(const Polynomial& polA,
                                                const Polynomial& polB) const {
  return moyalBracketHolger(polA, polB);
}

Polynomial
SemiClassicalLieAlgebra::moyalBracketHolger(const Polynomial& polA,
                                            const Polynomial& polB) const {
  //
  // This implementation of the Moyal bracket in terms of the product
  // involves a division by h-bar.  We actually achieve this by
  // testing that the power of h-bar is non-zero, and then
  // differentiating with respect to h-bar (which has the effect of
  // reducing the power by one), ignoring the resulting multiplier.
  //
  const Real tolerance(1.0e-15);
  Polynomial ab = moyalProductHolger(polA, polB);
  Polynomial ba = moyalProductHolger(polB, polA);
  Polynomial hBarResult = CoeffJ * (ab - ba);
  //now divide by hbar
  Polynomial result = zero();
  const PowersToCoeffMap& pc = hBarResult.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator iter;
  const Index iH = iHBar();
  for (iter = pc.begin(); iter != pc.end(); ++iter) {
    assert(isValidCoeff(iter->second));
    if (fabs(iter->second) > tolerance) {
      Power powerHBar = (iter->first)[iH];
      if (!(powerHBar > 0)) {
        std::cerr << "moyalBracketHolger: ";
        std::cerr << "zero power on h-bar, coeff: ";
        std::cerr << (iter->second);
      }
      //differentiate with respect to h-bar and discard any multiplier
      Powers powersDivHBar = ((iter->first).diff(iH)).second;
      result.setMonomial(powersDivHBar, iter->second);
    }
  }
  return result;
}

Polynomial
SemiClassicalLieAlgebra::moyalProductHolger(const Polynomial& polA,
                                            const Polynomial& polB) const {
  //
  // This implementation of the Moyal product is adapted from one by
  // Dr. Holger Waalkens, 2005.
  //
  if (!(hasElt(polA) && hasElt(polB)))
    throw LieAlgebraSizeMismatchError();
  const Power gradeA = grade(polA);
  const Power gradeB = grade(polB);
  const Power nMax = (gradeA <= gradeB) ? gradeA : gradeB;
  Polynomial result = zero();
  const PowersToCoeffMap& pcA = polA.getPowersAndCoeffs();
  const PowersToCoeffMap& pcB = polB.getPowersAndCoeffs();
  PowersToCoeffMap::const_iterator aIter;
  PowersToCoeffMap::const_iterator bIter;
  Coeff c;
  const Index iH = iHBar();
  for (Power n = 0; n <= nMax; ++n) {
    for (aIter = pcA.begin(); aIter != pcA.end(); ++aIter) {
      for (bIter = pcB.begin(); bIter != pcB.end(); ++bIter) {
        Powers m(getNumVars_());
        const Powers& poA = aIter->first;
        const Coeff& coA = aIter->second;
        const Powers& poB = bIter->first;
        const Coeff& coB = bIter->second;
        if (n == 0) {
          m = poA * poB;
          c = coA * coB;
          result.setMonomial(m, result[m] + c);
        }
        else {
          for (size_t N = 0; N < intPow((n+1), 2*dof_); ++N) {
            size_t Nnew = N;
            for (Index d = 0; d < dof_; ++d) {
              Index k = 2*dof_ - d;
              Index ii = (dof_ - d)-1;
              assert((0 <= ii) && (ii < dof_));
              m.setPower(iP(ii), Nnew / intPow(n+1, k-1));
              Nnew -= m[iP(ii)] * intPow(n+1, k-1);
            }
            for (Index d = 0; d < dof_; ++d) {
              Index k = dof_ - d;
              Index ii = k - 1;
              assert((0 <= ii) && (ii < dof_));
              m.setPower(iQ(ii), Nnew / intPow(n+1, k-1));
              Nnew -= m[iQ(ii)] * intPow(n+1, k-1);
            }
            Power diffOrder = 0;
            for (Index d = 0; d < dof_; ++d) {
              diffOrder += m[iQ(d)];
              diffOrder += m[iP(d)];
            }
            if (diffOrder == n) {
              c = coeffPow(Coeff(0.0, 0.5), n);
              Powers nMono(getNumVars_());
              nMono.setPower(iH, n + poA[iH] + poB[iH]);
              c *= coA * coB;
              for (Index d = 0; d < dof_; ++d) {
                c /= Coeff(factorial(m[iQ(d)]));
                c /= Coeff(factorial(m[iP(d)]));
              }
              Power signExp = 0;
              for (Index d = 0; d < dof_; ++d)
                signExp += m[iP(d)];
              c *= Coeff(intPow((-1), signExp));
              for (Index d = 0; d < dof_; ++d) {
                Power poAP = poA[iP(d)];
                Power poBQ = poB[iQ(d)];
                if ((m[iQ(d)] <= poAP) && (m[iQ(d)] <= poBQ)) {
                  nMono.setPower(iP(d), nMono[iP(d)] + poAP - m[iQ(d)]);
                  nMono.setPower(iQ(d), nMono[iQ(d)] + poBQ - m[iQ(d)]);
                  c *= Coeff(factorial(poAP))/Coeff(factorial(poAP-m[iQ(d)]));
                  c *= Coeff(factorial(poBQ))/Coeff(factorial(poBQ-m[iQ(d)]));
                }
                else
                  c = CoeffZero;
                Power poAQ = poA[iQ(d)];
                Power poBP = poB[iP(d)];
                if ((m[iP(d)] <= poAQ) && (m[iP(d)] <= poBP)) {
                  nMono.setPower(iQ(d), nMono[iQ(d)] + poAQ - m[iP(d)]);
                  nMono.setPower(iP(d), nMono[iP(d)] + poBP - m[iP(d)]);
                  c *= Coeff(factorial(poAQ))/Coeff(factorial(poAQ-m[iP(d)]));
                  c *= Coeff(factorial(poBP))/Coeff(factorial(poBP-m[iP(d)]));
                }
                else
                  c = CoeffZero;
              } //d
              if (isValidCoeff(c))
                result.addMonomial(nMono, c);
            } //diff_order==n
          } //for N
        } //case n!=0
      } //bIter
    } //aIter
  } //for n
  return result;
}

//def grade(self, elt):
//    """Overrides GradedInterface."""
//    self.check_elt(elt)
//    gra_pow = self.grade_of_powers
//    res = 0
//    for m, c in elt.powers_and_coefficients():
//        assert c != 0.0 #security
//        gra = gra_pow(m)
//        if gra > res:
//            res = gra
//    return res
//def is_isograde(self, elt, grade=-1):
//    """Overrides GradedInterface."""
//    self.check_elt(elt)
//    gra_pow = self.grade_of_powers
//    gra = grade
//    for m, c in elt.powers_and_coefficients():
//        assert c != 0.0 #security
//        if gra == -1:
//            gra = gra_pow(m)
//        else:
//            if gra != gra_pow(m):
//                return False
//    return True
//def isograde(self, elt, grade, up_to=None):
//    """Overrides GradedInterface."""
//    self.check_elt(elt)
//    gra_pow = self.grade_of_powers
//    if up_to==None:
//        up_to = grade+1
//    assert grade < up_to
//    res = self.zero()
//    for m, c in elt.powers_and_coefficients():
//        assert c != 0.0 #security
//        if grade <= gra_pow(m) < up_to:
//            res[m] = c
//    return res

//def h_bar(self, pow=1):
//    """Semiclassical h-bar variable."""
//    return self.coordinate_monomial(-1, pow)
