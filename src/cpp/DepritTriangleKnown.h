#ifndef DEPRIT_TRIANGLE_KNOWN__H
#define DEPRIT_TRIANGLE_KNOWN__H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangleKnown ---> not complete; pending python iter version.
//
// PURPOSE:
//
// Implements the Deprit TriangleKnown for normalization of a complex diagonal
// Hamiltonian with unknown generating function.
//
// NOTES:
//
// The Deprit triangle is one of the methods for computing the
// normalisation and coordinate changes.  There are others.  These
// various methods are referred to by Murdock as Normal Form Formats.
// Each format is either (a) iterative, or (b) recursive.  Each format
// is either (i) generative (working via a generating function), or
// (ii) direct.  The Deprit triangle, presented here, is of type
// (b)(i); recursive and generative.
//
// For now, this is a straightforward re-implementation of the Python
// version.
//
// TODO/FUTURE:
// 
// (1) I suggest making the deprit triangle accept an iterator-like
// object that provides a next() method for the input Hamiltonian
// (giving the terms of the correct grade --- note that we must still
// handle the factorial pre-factor within the triangle), and would
// itself provide a next() method so that the computation can be
// driven on-demand.
// 
// (2) Only the previous and current rows of the triangle (i.e., the
// intermediate quantities $H_{i}^{j}$ for (previous <= (i + j) <=
// current) are needed at each stage.  We can therefore safely delete
// (or release --- resource-wise) the earlier rows.
// 
// (3) The generating function list (i.e., the removed factorial
// pre-factor list of $W_i$'s) must be kept in order to perform the
// coordinate changes, etc.  However, for the Deprit triangle itself,
// we only require the previous generating function term [check this],
// so we could release earlier ones to a disk-based resource.
//
// (4) Ideally, we would pass a NormalFormStyle object which will define
// the partition functions, generator, solver, etc.
//
//----------------------------------------------------------------------

//standard headers
#include <utility>
#include <vector>

//library headers

//project headers
#include "Types.h"
#include "LieAlgebraBase.h"

//forward class declarations
class Polynomial;

//typedefs
typedef std::pair< Power, Power > TriangleAddress;

class DepritTriangleKnownError {
};
class DepritTriangleKnownGeneratorError : public DepritTriangleKnownError {
};

class DepritTriangleKnown {
 public:
  //public big four
  ~DepritTriangleKnown(void);
  //constructor
  DepritTriangleKnown(const LieAlgebraBase& algebra,
                      const std::vector< Polynomial >& innerGeneratorTerms);
  //inspectors
  inline const LieAlgebraBase& getAlgebra(void) const;
  inline const Power getCurrentGrade(void) const;
  inline const std::vector< Polynomial >& getInnerGeneratorTerms(void) const;
  //forward and backward lie transform
  inline void forward(const Polynomial& fTerm, Polynomial& resultFTerm);
  inline void reverse(const Polynomial& fTerm, Polynomial& resultFTerm);
 protected:
  //protected big four
  DepritTriangleKnown(void);
  DepritTriangleKnown(const DepritTriangleKnown& that);
  DepritTriangleKnown& operator=(const DepritTriangleKnown& that);
  //helper methods
  void makeReadyToComputeNextRow_(void);
  void checkGrades_(const Power i,
                    const Power j,
                    const Polynomial& fTerm,
                    const Polynomial& wTerm) const;
  void checkInnerGeneratorTerms_(void) const;
  void transform_(const Polynomial& fTerm,
                  Polynomial& resultfTerm,
                  const bool isBackward = 0);
  Polynomial triangle_(const Power i, const Power j, const bool isBackward);
 private:
  const LieAlgebraBase& algebra_;
  Power currentGrade_;
  std::map< TriangleAddress, Polynomial > fIJ_;
  std::vector< Polynomial > wI_;
}; //DepritTriangleKnown

//inline function implementations
const LieAlgebraBase& DepritTriangleKnown::getAlgebra(void) const {
  return algebra_;
}
const Power DepritTriangleKnown::getCurrentGrade(void) const {
  return currentGrade_;
}
const
std::vector< Polynomial >& DepritTriangleKnown::getInnerGeneratorTerms(void) const {
  return wI_;
}
void
DepritTriangleKnown::forward(const Polynomial& fTerm, Polynomial& resultFTerm) {
  transform_(fTerm, resultFTerm, false);
}
void
DepritTriangleKnown::reverse(const Polynomial& fTerm, Polynomial& resultFTerm) {
  transform_(fTerm, resultFTerm, true);
}

#endif //DEPRIT_TRIANGLE_KNOWN__H
