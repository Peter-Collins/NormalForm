//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Types implementation
//
//----------------------------------------------------------------------

//standard headers
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

//project headers
#include "Types.h"
#include "DefaultNumericalPrecision.h"

//pretty print coefficients
std::ostream& operator<<(std::ostream& ostr, const Complex& c) {
  ostr << std::setiosflags(std::ios_base::showpos | std::ios_base::left | std::ios_base::scientific);
#ifdef USE_GMP
  const size_t actualPrecisionDigits(numDigitsFromNumBits(c.real().get_prec()));
#else
  const size_t actualPrecisionDigits(15);
#endif //USE_GMP
  ostr << std::setprecision(actualPrecisionDigits);
  ostr << "(";
  //6 extra characters: +/-, ., e, +/-, 2.
  ostr << std::setw(actualPrecisionDigits+6) << c.real();
  ostr << " ";
  ostr << std::setw(actualPrecisionDigits+6) << c.imag();
  ostr << "J";
  ostr << ")";
  ostr << std::resetiosflags(std::ios_base::showpos);
  return ostr;
}

//
// Safe implementation of absolute value of complex number.  This
// avoids possible loss of precision and overflow errors.
//

Real fabs(const Complex& c) {
  const Real absR(fabs(c.real()));
  const Real absI(fabs(c.imag()));
  if (absR == RealZero)
    return absI;
  if (absI == RealZero)
    return absR;
  if (absR > absI) {
    const Real temp(absI / absR);
    return absR * sqrt(RealOne + temp * temp);
  }
  else {
    const Real temp(absR / absI);
    return absI * sqrt(RealOne + temp * temp);
  }
}

//random values

size_t randomSize(const size_t maxVal) {
  if (maxVal == 0)
    return 0;
  assert(maxVal > 0);
  size_t p = (rand()%(maxVal));
  assert(p >= 0);
  assert(p < maxVal);
  return p;
}

Index randomIndex(const Index maxVal) {
  if (maxVal == 0)
    return 0;
  assert(maxVal > 0);
  Index p = (rand()%(maxVal));
  assert(p >= 0);
  assert(p < maxVal);
  return p;
}

Power randomPower(const Power maxVal) {
  if (maxVal == 0)
    return 0;
  assert(maxVal > 0);
  Power p = (rand()%maxVal);
  assert(p >= 0);
  assert(p < maxVal);
  return p;
}

Real randomReal(const Real maxAbsVal) {
  assert(maxAbsVal >= 0.0);
  //ifndef USE_GMP
  Real a = Real(rand()%(MaxRandInt))/Real(MaxRandInt); //0.0--1.0
  Real result = (a-Real(0.5))*Real(2.0)*Real(maxAbsVal);
  //else
  //Real result;
  //mpf_urandomb(result, gmpRandState, actualPrecisionBits);
  //endif //USE_GMP
  assert(-maxAbsVal <= result);
  assert(result <= maxAbsVal);
  return result;
}

Complex randomComplex(const Real maxAbsVal) {
  Real re(randomReal(maxAbsVal));
  Real im(randomReal(maxAbsVal));
  return Complex(re, im);
}

Coeff randomCoeff(const Real maxAbsVal) {
  return randomComplex(maxAbsVal);
}
