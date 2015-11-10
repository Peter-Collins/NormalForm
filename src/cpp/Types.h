#ifndef TYPES_H
#define TYPES_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Types
//
// PURPOSE:
// 
// Here we define the basic types and utility functions that will be
// used throughout the code.  Some of this could be done more cleanly
// via templates, but for now this is a viable solution.
//
// NOTES:
//
// (1) We will assume that counting quantities are always size_t.
//
// (2) We will define a type Index for the number of variables, dof.
//
// (3) A type Power for the powers to which a variable can be raised,
//
// (4) A Real type, which may be a built-in or a GNU Multi-Precision,
//
// (5) A Complex type, derived from Real via std::complex,
//
// (6*) A Coeff type that will be used in Polynomials.
//
// Taking the trouble to stick to these types ensures consistent
// operation across all the modules.
//
// (6*) At present, the Coeff type is synonymous with Complex, but
// ultimately, one should probably template Polynomial over such a
// type, retrieving its Real type via std::complex< T >::value_type,
// or some other mechanism.  We assume that each coeff type can
// provide an absolute value function, Real fabs(const CoeffT&) and
// constants zero and one in order to allow comparisons.  Traits could
// be used to organise this.
//
// For example, template< CoeffT > class CoeffTraits { }; where we
// specialise to various types?  class CoeffTraits< Complex > {
// typedef Complex::value_type RealT; static const RealT
// RealZero(0.0); ...}?
//
// [1] My default random number generator:
//
// In the random-number generation routines, one must be careful with
// operator%, as this will throw a float error if the right operand is
// zero.  We return zero by default in this situation; an alternative
// would be to throw an exception for the unsuitable upper-bound.
//
// [2] GMP's random number generator:
//
// The following comment is relevant if we wish to use GMP's built-in
// random number generator.  Note that we must also be careful in
// generating random numbers for the GMP library; one needs to ensure
// that a correct random state object is defined and that we request
// enough bits.  At present, I have fudged this and just used ordinary
// random doubles, cast into GMP numbers.  It would be good to fix
// this.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>
#include <map>
#include <complex>

//library headers
#include <gmpxx.h>

//project headers

//flag to use gnu multi-precision library - set in Jamfile/compilation
#define USE_GMP

//---basic typdefs

//(use size_t for sizes)

//index of a variable in a monomial
typedef int Index;

//power of a variable
//at some points in the code, we take differences between powers!
//hence use signed type
typedef int Power;

//representation for reals
#ifdef USE_GMP
typedef mpf_class Real;
#else //USE_GMP
typedef double Real;
#endif //USE_GMP

//derived representation for complex
typedef std::complex<Real> Complex;

//coefficients in polynomials
typedef Complex Coeff;

//index to power map
typedef std::map< Index, Power > IndexToPowerMap;

//---static constants

static const Real RealZero(0.0);
static const Real RealOne(1.0);

static const Complex ComplexZero(0.0, 0.0);
static const Complex ComplexOne(1.0, 0.0);
static const Complex ComplexJ(0.0, 1.0);

static const Coeff CoeffZero(ComplexZero);
static const Coeff CoeffOne(ComplexOne);
static const Coeff CoeffJ(ComplexJ);

//maximum size random integer
static const size_t MaxRandInt = 54321;

//---function declarations

//pretty print coefficients?
std::ostream& operator<<(std::ostream& ostr, const Complex& c);

//generate random values for basic types
size_t randomSize(const size_t maxVal);
Index randomIndex(const Index maxVal);
Power randomPower(const Power maxVal);
Real randomReal(const Real maxAbsVal);
Complex randomComplex(const Real maxAbsVal);
Coeff randomCoeff(const Real maxAbsVal);

//absolute value of coefficients
#ifdef USE_GMP
inline Real fabs(const Real& r) {
  if (r >= RealZero)
    return r;
  else
    return -r;
}
#endif //USE_GMP

class TypesError {
};
class TypesOverflowError : public TypesError {
};

//---inline function declarations

inline Real square(const Real& r);
Real fabs(const Complex& c);
inline bool isNonnegativePower(const Power p);
inline bool isPositivePower(const Power p);
inline bool isValidCoeff(const Coeff& coeff);

//---inline function bodies

inline
bool isNonnegativePower(const Power p) {
  return (p >= 0);
}

inline
bool isPositivePower(const Power p) {
  return (p > 0);
}

//
// Here, a valid coeff means one that we take the trouble to store.
// the following can be modified to introduce a numerical tolerance
//
inline
bool isValidCoeff(const Coeff& coeff) {
  return (coeff != CoeffZero);
  //return (fabs(coeff) > 1.0e-14);
}

#endif //TYPES_H

