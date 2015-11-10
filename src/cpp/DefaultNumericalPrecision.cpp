//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: DefaultNumericalPrecision implementation
//
//----------------------------------------------------------------------

//standard headers
#include <assert.h>
#include <iostream>

//project headers
#include "Types.h"
#include "DefaultNumericalPrecision.h"

//namespaces
using std::cerr;
using std::endl;

//constructor given desired number of bits
DefaultNumericalPrecision::DefaultNumericalPrecision(size_t desiredNumBits) {
#ifdef USE_GMP
  previousNumBits_ = mpf_get_default_prec();
  mpf_set_default_prec(desiredNumBits);
  currentNumBits_ = mpf_get_default_prec();
  assert(currentNumBits_ >= desiredNumBits);
#else //USE_GMP
  previousNumBits_ = 53;
  currentNumBits_ = previousNumBits_;
  if (currentNumBits_ < desiredNumBits)
    cerr << "WARNING: USING DOUBLE PRECISION!" << endl;
#endif //USE_GMP
}

DefaultNumericalPrecision::~DefaultNumericalPrecision(void) {
#ifdef USE_GMP
  mpf_set_default_prec(previousNumBits_);
  currentNumBits_ = mpf_get_default_prec();
  assert(currentNumBits_ >= previousNumBits_);
#endif  
}

std::ostream& operator<<(std::ostream& ostr,
                         const DefaultNumericalPrecision& prec) {
  ostr << "DefaultNumericalPrecision: ";
  ostr << prec.getNumBits() << " bits = approx ";
  ostr << numDigitsFromNumBits(prec.getNumBits()) << " digits." << endl;
  return ostr;
}

