/// @file ClassicalLieAlgebra.cpp
///
/// @brief ClassicalLieAlgebra implementation
///
/// The Lie bracket here is the usual Poisson bracket.  The same
/// implementation will work for semiclassical Poisson bracket which
/// could be used to test the lowest order term of the semiclassical
/// Moyal bracket; this would make a useful additional test.
///
/// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

//standard headers

//library headers

//project headers
#include "ClassicalLieAlgebra.h"

/// the big four: default constructor
ClassicalLieAlgebra::ClassicalLieAlgebra(void)
  : LieAlgebraBase() {
}
/// the big four: copy constructor
ClassicalLieAlgebra::ClassicalLieAlgebra(const ClassicalLieAlgebra& that)
  : LieAlgebraBase(that) {
}
/// the big four: destructor
ClassicalLieAlgebra::~ClassicalLieAlgebra(void) {
}

/// constructor from the number of degrees of freedom
ClassicalLieAlgebra::ClassicalLieAlgebra(const Index dof)
  : LieAlgebraBase(dof) {
}

/// the Poisson bracket of two polynomials in this algebra
Polynomial ClassicalLieAlgebra::poissonBracket(const Polynomial& polA,
                                               const Polynomial& polB) const {
  if (!(hasElt(polA) && hasElt(polB)))
    throw LieAlgebraSizeMismatchError();
  assert(polA.hasSameNumVars(polB));
  Polynomial pb(zero());
  Polynomial pbA(zero());
  Polynomial pbB(zero());
  for (Index d = 0; d < dof_; ++d) {
    pbA = (polA.diff(iQ(d)) * polB.diff(iP(d)));
    pbB = (polA.diff(iP(d)) * polB.diff(iQ(d)));
    pb += (pbA - pbB);
  }
  return pb;
}
