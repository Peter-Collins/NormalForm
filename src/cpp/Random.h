//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef NF_RANDOM_H
#define NF_RANDOM_H

//system
#include <vector>

//library

//project
#include "Types.h"
#include "MapPowers.h"
#include "Polynomial.h"

std::vector< Power > randomVectorPowers(const Index length,
					const Power maxDegree);

MapPowers randomPowers(const Index length, const Power maxDegree);

Polynomial randomPoly(const Index numVars=0,
		      const Index maxNumVars=40,
		      const Power maxDegree=12,
		      const size_t maxNumTerms=32,
		      const Real maxAbsCoeff=2.0);

#endif //NF_RANDOM_H
