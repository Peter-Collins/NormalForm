//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: LieAlgebraTestBase implementation
//
//----------------------------------------------------------------------

//standard headers
#include <vector>
#include <assert.h>

//library headers

//project headers
#include "LieAlgebraTestBase.h"

std::vector< Coeff > makeVector(const Index numVars,
                                const Coeff val) {
  std::vector< Coeff > vec;
  for (Index i = 0; i < numVars; ++i)
    vec.push_back(val);
  assert(vec.size() == numVars);
  return vec;
}
