#ifndef UTILITY_H
#define UTILITY_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Utility
//
// PURPOSE:
// 
// Here we define the basic utility functions that will be used
// throughout the code.  Some of this could be done more cleanly via
// templates, but for now this is a viable solution.
//
// NOTES:
//
// We should probably replace the separate coeffPow and intPow
// functions by a single templated function.  For now, the current
// solution is fine.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>

//library headers
#include <gmpxx.h>

//project headers
#include "Types.h"

//exceptions
class UtilityError { };
class UtilityNegativeArgumentError : UtilityError { };
class UtilityOutOfRangeError : UtilityError { };

//---function declarations

//efficiently raise a coefficient to a power by repeated squaring
Coeff coeffPow(const Coeff c, const Power p);

//do the same for integers; we must allow negative values for i
int intPow(const int i, const size_t p);

size_t factorial(const int i);

size_t binomial(const int n, const int k);

#endif //UTILITY_H
