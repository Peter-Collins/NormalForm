#ifndef SEMI_CLASSICAL_LIE_ALGEBRA_H
#define SEMI_CLASSICAL_LIE_ALGEBRA_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: SemiClassicalLieAlgebra
//
// PURPOSE:
//
// Implement the semi-classical Lie algebra by embedding the q and p
// components into a polynomial ring with an additional h-bar variable
// and defining the graded structure in terms of the homogeneous
// polynomial structure of that ring plus counting h-bar twice.
//
// NOTES:
//
// We could add a constructor from a classical algebra to do the
// quantization; think about this, and what mechanism would be used
// for conversion between the two.  In principle, it should allow for
// different variable encodings?
//
// Note: we should define a gradeOfPowers public method in the
// LieAlgebraBase which calls a pure virtual gradeOfPowers_.  We then
// implement grade, isIsoGrade, isoGradePart, all in terms of this
// gradeOfPowers_ function.  This refactors the graded structure code
// up into the LieAlgebraBase.
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers
#include "LieAlgebraBase.h"

class SemiClassicalLieAlgebra : public LieAlgebraBase {
 public:
  //big four
  SemiClassicalLieAlgebra(void);
  SemiClassicalLieAlgebra(const SemiClassicalLieAlgebra& that);
  ~SemiClassicalLieAlgebra(void);

  //constructors
  explicit SemiClassicalLieAlgebra(const Index dof);

  //semi-classical structure via h-bar variable
  inline Index iHBar(void) const;
  inline Polynomial hBar(const Power p = 1) const;

  //some lie bracket implementations
  Polynomial poissonBracket(const Polynomial& polA,
                            const Polynomial& polB) const;
  Polynomial moyalBracketHolger(const Polynomial& polA,
                                const Polynomial& polB) const;
  Polynomial moyalProductHolger(const Polynomial& polA,
                                const Polynomial& polB) const;
 protected:
  //forbid assignment
  SemiClassicalLieAlgebra& operator=(const SemiClassicalLieAlgebra& that);

  //implement PolynomialRingInterface
  inline Index getNumVars_(void) const;

  //implement LieAlgebraBase; graded structure
  inline Power gradeOfPowers_(const Powers& powers) const;

  //implement LieAlgebraBase; lie bracket
  Polynomial lieBracket_(const Polynomial& polA,
                         const Polynomial& polB) const;
 private:
};

inline
Index SemiClassicalLieAlgebra::getNumVars_(void) const {
  return 2*dof_ + 1;
}

inline
Index SemiClassicalLieAlgebra::iHBar(void) const {
  //index of the h-bar variable; we use the final index
  return 2*dof_;
}

inline
Polynomial SemiClassicalLieAlgebra::hBar(const Power p) const {
  return coordinateMonomial(iHBar(), p);
}

inline
Power
SemiClassicalLieAlgebra::gradeOfPowers_(const Powers& powers) const {
  return powers.degree() + powers[iHBar()];
}

#endif //SEMI_CLASSICAL_LIE_ALGEBRA_H
