#ifndef CLASSICAL_LIE_ALGEBRA_H
#define CLASSICAL_LIE_ALGEBRA_H

/// @file ClassicalLieAlgebra.h
///
/// @brief Subclass Classical Lie algebras from the general base.
///
/// Implement the classical Lie algebra by embedding the q and p
/// components into a polynomial ring and defining the graded structure
/// in terms of the homogeneous polynomial structure of that ring.
///
/// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

//standard headers

//library headers

//project headers
#include "LieAlgebraBase.h"

/// Implement Classical Lie algebras by subclassing from a general Lie algebra base.
class ClassicalLieAlgebra : public LieAlgebraBase {
 public:
  //big four
  ClassicalLieAlgebra(void);
  ClassicalLieAlgebra(const ClassicalLieAlgebra& that);
  ~ClassicalLieAlgebra(void);
  //constructors
  explicit ClassicalLieAlgebra(const Index dof);
  //poisson bracket
  Polynomial poissonBracket(const Polynomial& polA, const Polynomial& polB) const;
 protected:
  //forbid assignment
  ClassicalLieAlgebra& operator=(const ClassicalLieAlgebra& that);
  //implement PolynomialRingInterface
  inline Index getNumVars_(void) const;
  //implement LieAlgebraBase
  inline Power gradeOfPowers_(const Powers& powers) const;
  inline Polynomial lieBracket_(const Polynomial& polA,
                                const Polynomial& polB) const;
 private:
};

/// there are twice as many variables (positions and momenta) as degrees of freedom
inline
Index ClassicalLieAlgebra::getNumVars_(void) const {
  return 2*dof_;
}

/// the graded structure here is simply that of the degree of homogeneous polynomials
inline
Power ClassicalLieAlgebra::gradeOfPowers_(const Powers& powers) const {
  return powers.degree();
}

/// the Lie bracket is the usual Poisson bracket
inline
Polynomial ClassicalLieAlgebra::lieBracket_(const Polynomial& polA,
                                            const Polynomial& polB) const {
  return poissonBracket(polA, polB);
}

#endif //CLASSICAL_LIE_ALGEBRA_H
