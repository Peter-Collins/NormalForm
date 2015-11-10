#ifndef COORDINATE_CHANGE_H
#define COORDINATE_CHANGE_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: CoordinateChange
//
// PURPOSE:
//
// Implements the Deprit Triangle for a known generator, for the
// transformation of coordinate selector functions.
//
// NOTES:
//
// The following observations note _IMPORTANT_ differences between the
// method of Lie series exemplified by Deprit's method, and other
// methods (for example, Hori's method, the method of Dragt-Finn,
// etc.)  See reference [1], below, for more details.
//
// [a(i)] In Deprit's method, the "generator" w appearing in the
// Deprit triangle is _NOT_ the generator of the flow between the two
// coordinate systems in the usual sense.
//
// [a(ii)] In particular, the "generator" w is _NOT_ in general
// invariant under the coordinate transformation defined by the Deprit
// triangle!  It is very important to understand this point.
//
// [a(iii)] Instead, the proper test for correctness of implementation
// is one that tests both the normalisation Deprit triangle _and_ the
// direct coordinate change triangle _together_.  In other words, let
// the normalisation triangle applied to h yield h' with "generator"
// w, and the direct coordinate change triangle applied to x yield x'
// under the same "generator"; then the test for correctness is h'(x')
// == h(x).  Suppose that the inverse triangle (using the same
// "direct" generator) maps y' to y, then we test this via y(x') == x.
//
// [b(i)] In Deprit's method, the negative (-w) of the "generator" (w)
// is _NOT_ the generator of the inverse (transform/flow).
//
// [b(ii)] Instead, one must apply the Deprit triangle recurrence "in
// reverse" using the original (direct) "generator" in order to get
// the inverse transformation.  This is the reason for the methods
// forwardTriangle_ and reverseTriangle_ in what follows.
//
// [c(i)] In light of the above, one cannot transform arbitrary scalar
// functions using this method.  It is designed for the coordinate
// selectors only in a near-identity transformation (see next
// comment), where the "generator" (w) has been obtained specifically
// from the Deprit triangle used for normalisation.
//
// [c(ii)] As an example of the differences between the Deprit
// triangle and a direct computation of the flow under a generator,
// consider the generator (in the usual sense) q0p0.  The image of the
// coordinate selectors (q0, p0) under the forward time-1 flow for
// this generator are (q0.e, p0/e), where e = exp(1).  These are not
// near-identity maps, since they should take the form (q0 + h.o.t.,
// p0 + h.o.t.).  Also, note that it is perfectly fine to perform the
// coordinate selector transforms via [exp(ad_W)], giving the results
// above, _BUT_ this is not a grade-by-grade computation in the
// following sense: computation of the linear terms q0.e and q1/e
// actually requires access to _all_ terms of the series for
// exp(ad_W).
//
// [c(iii)] One can construct particular example of cubic generators
// for which the Deprit triangle does indeed compute the time-1 map of
// the flow; however, these are an exception rather than the rule;
// they have special structure since there is only a single order of
// epsilon (the "small parameter") present.
//
// REFERENCES:
//
// [1] For references to the above issues with Deprit's method, please
// see, for example, Murdock "Normal forms and unfoldings for local
// dynamical systems", pp. 175, 176, 182, 464, 470.
//
// [2] Please see the comments in DepritTriangle.h, many of which
// pertain equally well to this module.
//
// TODO:
//
// Think about the grade-by-grade nature of things here, with regard
// to reorganising and parallelising the code.
//
// Note that we need to clear the entries of the triangle for each
// transformation.  This hints that the design could be changed
// slightly, with the main class spawning a TriangleInstance whose
// constructor always creates blank triangle elements?
//
// Input for the "generator" terms could be a PolynomialGrades object
// which provides a next() method to iterate through grades?
// Particular implementations of this would then be, for example, an
// ordinary std::list, a filestream, a socket- (or other communication
// channel-) driven load-on-demand?
//
// Write now using std structures; template later!
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

//exceptions
class CoordinateChangeError { };
class CoordinateChangeGradeError : public CoordinateChangeError { };

class CoordinateChange {
 public:
  //public big four
  ~CoordinateChange(void);
  //constructor
  CoordinateChange(const LieAlgebraBase& algebra,
                   const std::vector< Polynomial >& wI,
                   const bool isReverse = false);
  //inspectors
  inline const LieAlgebraBase& getAlgebra(void) const;
  inline const Power getCurrentGrade(void) const;
  //forward transform and reverse transform of a scalar function
  void computeNextTerm(const Polynomial& nextInputTerm,
                       Polynomial& result);
 protected:
  //protected big four
  CoordinateChange(void);
  CoordinateChange(const CoordinateChange& that);
  CoordinateChange& operator=(const CoordinateChange& that);
  //helper methods
  void makeReadyToComputeNextRow_(void);
  void checkInnerGeneratorTerms_(void) const;
  void checkGrades_(const Power i,
                    const Power j,
                    const Polynomial& fTerm,
                    const Polynomial& wTerm) const;
  void forwardCoordinateChange_(const Polynomial& nextFTerm,
                                Polynomial& nextResultFTerm);
  const Polynomial& forwardTriangle_(void);
  void reverseCoordinateChange_(const Polynomial& nextFTerm,
                                Polynomial& nextResultFTerm);
  const Polynomial& reverseTriangle_(void);
 private:
  const LieAlgebraBase& algebra_;
  const std::vector< Polynomial >& wI_;
  const bool isReverse_;
  Power currentGrade_;
  std::map< TriangleAddress, Polynomial > fIJ_;
}; //CoordinateChange

//inline function implementations
const LieAlgebraBase& CoordinateChange::getAlgebra(void) const {
  return algebra_;
}
const Power CoordinateChange::getCurrentGrade(void) const {
  return currentGrade_;
}

#endif //COORDINATE_CHANGE_H
