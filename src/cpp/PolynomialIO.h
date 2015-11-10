#ifndef POLYNOMIAL_IO_H
#define POLYNOMIAL_IO_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialIO
//
// PURPOSE:
//
// Implement input/output routines for polynomials in
// s-expression-like format.
//
// NOTES:
//
// It is very important to have the correct numerical precision set if
// using the GMP library, in order to ensure that all the necessary
// digits are read.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>

//library headers

//project headers
#include "Types.h"
#include "Polynomial.h"

void seWritePowers(std::ostream& os, const Powers& powers);
void seWriteCoeff(std::ostream& os, const Coeff& coeff);
void seWriteMonomial(std::ostream& os, const Powers& powers, const Coeff& coeff);
void seWritePolynomial(std::ostream& os, const Polynomial& poly);
void seWriteVectorOfPolynomials(std::ostream& os, const std::vector< Polynomial >& polys);

Powers seReadPowers(std::istream& is);
Coeff seReadCoeff(std::istream& is);
Polynomial& seReadMonomial(std::istream& is, Polynomial& result);
Polynomial seReadPolynomial(std::istream& is);
void seReadVectorOfPolynomials(std::istream& os, std::vector< Polynomial >& result);

#endif //POLYNOMIAL_IO_H
