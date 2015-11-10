#ifndef DEPRIT_TRIANGLE_H
#define DEPRIT_TRIANGLE_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DepritTriangle
//
// PURPOSE:
//
// Implements the Deprit Triangle for normalization of a complex
// diagonal Hamiltonian via solution of the corresponding
// inhomogeneous homological equation, producing as a side-effect the
// "generator" (see important notes, below).
//
// NOTES:
//
// [a(i)] The Deprit triangle is one of the methods for computing the
// normalisation and coordinate changes.  There are others.  These
// various methods are referred to by Murdock as Normal Form Formats.
// Each format is either (a) iterative, or (b) recursive.  Each format
// is either (i) generative (working via a generating function, _NOT_
// necessarily in the sense of the generator of the coordinate-change
// flow), or (ii) direct.
//
// [a(ii)] The Deprit triangle, presented here, is of type (b)(i);
// recursive and generative.
//
// [b] For now, this is a re-implementation of the Python version with
// the added benefit of iterator-like structure for efficiency (we
// pass in each term of successively higher grade and yield the
// corresponding normalised term and generator term).
//
// The following observations note _IMPORTANT_ differences between the
// method of Lie series exemplified by Deprit's method, and other
// methods (for example, Hori's method, the method of Dragt-Finn,
// etc.)  See reference [1], below, for more details.
//
// It is especially important to understand that expectations
// involving the generator, that one might have from knowledge of the
// other methods, often do _NOT_ hold for the "generator" in Deprit's
// method!
//
// [c(i)] In Deprit's method, the "generator" w appearing in the
// Deprit triangle is _NOT_ the generator of the flow between the two
// coordinate systems in the usual sense.
//
// [c(ii)] In particular, the "generator" w is _NOT_ in general
// invariant under the coordinate transformation defined by the Deprit
// triangle!  It is very important to understand this point.
//
// [c(iii)] Instead, the proper test for correctness of implementation
// is one that tests both the normalisation Deprit triangle _and_ the
// direct coordinate change triangle _together_.  In other words, let
// the normalisation triangle applied to h yield h' with "generator"
// w, and the direct coordinate change triangle applied to x yield x'
// under the same "generator"; then the test for correctness is h'(x')
// == h(x).  Suppose that the inverse triangle (using the same
// "direct" generator) maps y' to y, then we test this via y(x') == x.
//
// [d(i)] In Deprit's method, the negative (-w) of the "generator" (w)
// is _NOT_ the generator of the inverse (transform/flow).
//
// [d(ii)] Instead, one must apply the Deprit triangle recurrence "in
// reverse" using the original (direct) "generator" in order to get
// the inverse transformation.  This is the reason for the methods
// forwardTriangle_ and reverseTriangle_ in what follows.
//
// REFERENCES:
//
// [1] For references to the above issues with Deprit's method, please
// see, for example, Murdock "Normal forms and unfoldings for local
// dynamical systems", pp. 175, 176, 182, 464, 470.
//
// [2] Further comments on the above issues with Deprit's method are
// in the module CoordinateChange.
//
// TODO/FUTURE:
// 
// [a] I suggest making the deprit triangle accept an iterator-like
// object that provides a next() method for the input Hamiltonian
// (giving the terms of the correct grade --- note that we must still
// handle the factorial pre-factor within the triangle), and would
// itself provide a next() method so that the computation can be
// driven on-demand.
// 
// [b(i)] Ideally, we would pass a NormalFormStyle object which will
// define the partition functions, generator, solver, etc.
//
// [b(ii)] In fact, styles that would be useful in our current work,
// on reaction dynamics, are:-
//
// - straightforward semisimple normal form,
//
// - inner-product / Belitskii / Elphick et. al.
//
// At present, the sl(2) style, although computationally more powerful
// (especially for large systems) does not look practicable for our
// purposes.  We should review the truth of this statement later.
//
// [c] In view of the issues involved with Deprit's method (lack of
// "generator" invariance, inability to use the negative of the
// generator as the generator of the inverse, and the proliferation of
// factorial and binomial-coefficient prefactors, it would be
// advisable at some stage to change to Hori's method (for the most
// similar implementation) or to Dragt-Finn (which is different in
// structure and has the benefit of fewer Lie bracket evaluations).
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

class DepritTriangle {
 public:

  //public big four
  ~DepritTriangle(void);

  //constructor
  DepritTriangle(const LieAlgebraBase& algebra, const Polynomial& quadratic);

  //inspectors
  inline const LieAlgebraBase& getAlgebra(void) const;
  inline const Polynomial& getQuadratic(void) const;
  inline const std::vector< Coeff >& getFrequencies(void) const;
  inline const Power getCurrentGrade(void) const;
  inline const std::vector< Polynomial >& getInnerGeneratorTerms(void) const;

  //methods to compute the normalisation and generating function
  void computeNormalFormAndGenerator(const Polynomial& hTerm,
                                     Polynomial& resultKTerm,
                                     Polynomial& resultWTerm);

 protected:

  //protected big four
  DepritTriangle(void);
  DepritTriangle(const DepritTriangle& that);
  DepritTriangle& operator=(const DepritTriangle& that);

  //helper methods
  void checkQuadraticPart_(void) const;
  void extractFundamentalFrequencies_(void);
  void makeReadyToComputeNextRow_(void);
  inline bool isInNormalForm_(const Powers& powers) const;
  void partitionKnownTerms_(const Polynomial& hTerm,
                            Polynomial& resultNormal,
                            Polynomial& resultRemainder) const;
  void checkGrades_(const Power i,
                    const Power j,
                    const Polynomial& hTerm,
                    const Polynomial& wTerm) const;
  Polynomial solveHomologicalEquation_(const Polynomial& remainder) const;
  void checkHomologicalEquation_(const Polynomial& remainder) const;
  Polynomial triangleMinusUnknownTerm_(const Power i, const Power j);
  void correctRowUsingPartition_(const Polynomial& normal,
                                 const Polynomial& remainder);

 private:
  const LieAlgebraBase& algebra_;
  const Polynomial& quadratic_;
  std::vector< Coeff > frequencies_;
  Power currentGrade_;
  std::map< TriangleAddress, Polynomial > hIJ_;
  std::vector< Polynomial > wI_;

}; //DepritTriangle

//inline function implementations
const LieAlgebraBase& DepritTriangle::getAlgebra(void) const {
  return algebra_;
}
const Polynomial& DepritTriangle::getQuadratic(void) const {
  return quadratic_;
}
const std::vector< Coeff >& DepritTriangle::getFrequencies(void) const {
  return frequencies_;
}
const Power DepritTriangle::getCurrentGrade(void) const {
  return currentGrade_;
}
const
std::vector< Polynomial >& DepritTriangle::getInnerGeneratorTerms(void) const {
  return wI_;
}
bool DepritTriangle::isInNormalForm_(const Powers& powers) const {
  for (Index i = 0; i < algebra_.getDof(); ++i)
    if (powers[algebra_.iQ(i)] != powers[algebra_.iP(i)])
      return false;
  return true;
}

#endif //DEPRIT_TRIANGLE_H
