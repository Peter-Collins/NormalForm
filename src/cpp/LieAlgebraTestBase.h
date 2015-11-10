#ifndef LIE_ALGEBRA_TEST_BASE
#define LIE_ALGEBRA_TEST_BASE

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
//----------------------------------------------------------------------

//standard headers
#include <vector>

//library headers

//project headers
#include "Types.h"

std::vector< Coeff > makeVector(const Index numVars,
                                const Coeff val);

#endif //LIE_ALGEBRA_TEST_BASE
