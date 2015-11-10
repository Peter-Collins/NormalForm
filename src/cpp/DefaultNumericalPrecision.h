#ifndef DEFAULT_NUMERICAL_PRECISION_H
#define DEFAULT_NUMERICAL_PRECISION_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DefaultNumericalPrecision
//
// PURPOSE:
// 
// Here, we provide a very simple facility for automatically-scoped
// precision changes to the GMP library.
//
// NOTES:
//
// [1] Usage:
//
// In order to use this facility, one simply creates a Precision
// object within a scope.  The constructor records the previous GMP
// precision and sets the new one.  As the object goes out of scope,
// the destructor automatically reverses this process, restoring the
// previous numerical precision.
//
// In the case where one only wants to use a single particular precision
// throughout an executable, one could instantiate a static variable in
// a single linked module.  The C++ standard demands that the
// compiler should perform all such intitialisations before execution of
// the main function begins.
//
// Just in case of compiler evils, we include a trace output to confirm
// that things are really operating as expected.
//
// [2] Using ordinary double precision numbers:
//
// If double precision is actually needed, one should unset the symbol
// USE_GMP in Types.h and recompile the code (among other things, this
// will instantiate the correct templates over doubles, including
// std::complex< T >).
//
// [3] Bits vs digits of precision in the mantissa of a float:
//
// For standard double precision numbers, we expect 64 bits, 53 of
// which are used in the mantissa [see, for example, "High-Precision
// Floating-Point Arithmetic in Scientific Computation", David
// H. Bailey, 28 January 2005, and the references therein].
//
// Converting to digits, via log(10)/log(2), this gives us
// approximately 15 digits.  I note this, here, because we
// need to know about the number of digits in order to
// perform certain formatted I/O operations without accidental
// truncation or output of spurious extra digits.
//
// At present, I assume that if I ask GMP for an N-bit float, then all
// N bits are used for the mantissa of the number.  I do the
// conversion to digits based on this assumption.  It would be advisable
// to confirm this in the documentation, but this assumption is
// certainly borne out by the code that I have used so far.
//
// [4] Miscellaneous:
//
// We will forbid the default constructor.
//
// GOTCHAS:
//
// A major gotcha here is thread safety.  I have not thought deeply
// about this.  For now, we are in single-threaded code, so that the
// issue does not arise.
//
// ASSOCIATIONS:
//
// Please see the associated Types module for reference to conditional
// compilation with the GNU Multi-Precision library and, for example,
// the associated ostream& operator<<(const Complex&) operator.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>

//library headers

//project headers

//exceptions
class DefaultNumericalPrecisionError {
};

class DefaultNumericalPrecision
{
 public:
  //big four
  ~DefaultNumericalPrecision(void);
  //constructor
  explicit DefaultNumericalPrecision(size_t desiredNumBits);
  //inspectors
  inline size_t getNumBits(void) const;
  inline size_t getNumDigits(void) const;
 protected:
  //big four
  DefaultNumericalPrecision(void);
  DefaultNumericalPrecision(const DefaultNumericalPrecision& that);
  DefaultNumericalPrecision& operator=(const DefaultNumericalPrecision& that);
 private:
  size_t previousNumBits_;
  size_t currentNumBits_;
};

//report
std::ostream& operator<<(std::ostream& ostr,
                         const DefaultNumericalPrecision& prec);

//inline function declarations
inline size_t numDigitsFromNumBits(size_t numBits);

//inline member function implementations
size_t DefaultNumericalPrecision::getNumBits(void) const {
  return currentNumBits_;
}
size_t DefaultNumericalPrecision::getNumDigits(void) const {
  return numDigitsFromNumBits(currentNumBits_);
}

//inline function implementations
inline size_t numDigitsFromNumBits(size_t numBits) {
  return size_t(double(numBits / 3.321928094887));
}

#endif //DEFAULT_NUMERICAL_PRECISION_H

