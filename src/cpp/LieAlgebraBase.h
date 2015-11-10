#ifndef LIE_ALGEBRA_BASE
#define LIE_ALGEBRA_BASE

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: LieAlgebraBase
//
// Define the base class (_not_ interface; it has data members!) to
// all _graded_ lie algebras in our code.  This differs from the
// Python implementation: for simplicity, I do not introduce the
// concept of grade into PolynomialRing, only at this stage in the
// LieAlgebra.  This removes an extra level of the inheritance
// hierarchy.
//
//----------------------------------------------------------------------

//standard headers
#include <utility>

//library headers

//project headers
#include "Types.h"
#include "PolynomialRingInterface.h"

class LieAlgebraError {
};
class LieAlgebraIndexError : public LieAlgebraError {
};
class LieAlgebraPowerError : public LieAlgebraError {
};
class LieAlgebraSizeMismatchError : public LieAlgebraError {
};

class LieAlgebraBase : public PolynomialRingInterface {

 public:

  //inspectors
  inline Index getDof(void) const;

  //encoding of q, p into x; for now, we implement here
  inline Index iQ(const Index d) const;
  inline Index iP(const Index d) const;

  //coordinate monomials;
  inline Polynomial q(const Index d, const Power power = 1) const;
  inline Polynomial p(const Index d, const Power power = 1) const;

  //diagonality
  bool isDiagonal(const Powers& powers) const;
  bool isDiagonal(const Polynomial& pol) const;
  void diagonalAndNonDiagonal(const Polynomial& pol,
                              std::pair< Polynomial, Polynomial >& dAndN) const;

  //graded structure; public interface to protected pure virtual
  inline Power gradeOfPowers(const Powers& powers) const;

  //the following can actually be implemented in terms of the above!
  Power grade(const Polynomial& pol) const;
  bool isIsoGrade(const Polynomial& pol) const;
  bool isIsoGrade(const Polynomial& pol, const Power gra) const;
  Polynomial isoGradePart(const Polynomial& pol,
                          const Power p,
                          const Power lessThan) const;
  inline Polynomial isoGradePart(const Polynomial& pol,
                                 const Power p) const;

  //the lie bracket
  inline Polynomial lieBracket(const Polynomial& polA,
                               const Polynomial& polB) const;
  
 protected:

  //the big four with virtual destructor
  LieAlgebraBase(void);
  LieAlgebraBase(const LieAlgebraBase& that);
  LieAlgebraBase& operator=(const LieAlgebraBase& that);
  virtual ~LieAlgebraBase(void);

  //protected constructor from dof
  explicit LieAlgebraBase(const Index dof);

  //graded structure; pure virtual
  virtual Power gradeOfPowers_(const Powers& powers) const = 0;

  //lie bracket; pure virtual
  virtual Polynomial lieBracket_(const Polynomial& polA,
                                 const Polynomial& polB) const = 0;

  //data members
  const Index dof_;

 private:
};

inline
Index LieAlgebraBase::getDof(void) const {
  return dof_;
}

inline
Index LieAlgebraBase::iQ(const Index d) const {
  if (!(d >= 0 && d < dof_))
    throw LieAlgebraIndexError();
  return 2*d;
}

inline
Index LieAlgebraBase::iP(const Index d) const {
  if (!(d >= 0 && d < dof_))
    throw LieAlgebraIndexError();
  return 2*d + 1;
}

inline
Polynomial LieAlgebraBase::q(const Index d, const Power power) const {
  return coordinateMonomial(iQ(d), power);
}

inline
Polynomial LieAlgebraBase::p(const Index d, const Power power) const {
  return coordinateMonomial(iP(d), power);
}

inline
Power LieAlgebraBase::gradeOfPowers(const Powers& powers) const {
  return gradeOfPowers_(powers);
}

inline
Polynomial LieAlgebraBase::lieBracket(const Polynomial& polA,
                                      const Polynomial& polB) const {
  if (!(hasElt(polA) && hasElt(polB)))
    throw LieAlgebraSizeMismatchError();
  return lieBracket_(polA, polB);
}

inline
Polynomial LieAlgebraBase::isoGradePart(const Polynomial& pol,
                                        const Power p) const {
  return isoGradePart(pol, p, p + 1);
}

#endif //LIE_ALGEBRA_BASE
