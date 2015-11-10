//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: PolynomialIO implementation.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>
#include <iomanip>
#include <vector>

//library headers

//project headers
#include "Types.h"
#include "DefaultNumericalPrecision.h"
#include "SExprIO.h"
#include "PolynomialIO.h"

//template instantiations
#include "SExprIO.cpp"
template
std::vector< Power > seReadVector< Power >(std::istream& is, const size_t len);

class PolynomialIOError { };

class PolynomialIOSizeError : PolynomialIOError { };

class PolynomialIOValueError : PolynomialIOError { };

class PolynomialIORanOutOfInputError : PolynomialIOError { };

void seWritePowers(std::ostream& os, const Powers& powers) {
  os << "(";
  for (Index i = 0; i < powers.getNumVars(); ++i) {
    if ((i >= 1) and (i < powers.getNumVars()))
      os << " ";
    os << powers[i];
  }
  os << ")";
}

Powers seReadPowers(std::istream& is, const int numVars) {
  Powers powers(seReadVector< Power >(is, numVars));
  if (!numVars == powers.getNumVars())
    throw PolynomialIOSizeError();
  return powers;
}

void seWriteCoeff(std::ostream& os, const Coeff& coeff) {
  //this needs to be written in a stable form for double / gmp
  os << std::setiosflags(std::ios_base::showpos);
  os << std::setiosflags(std::ios_base::left);
  os << std::setiosflags(std::ios_base::scientific);
#ifdef USE_GMP
  const size_t actualPrecisionDigits(numDigitsFromNumBits(coeff.real().get_prec()));
#else
  const size_t actualPrecisionDigits(15);
#endif //USE_GMP
  os << std::setprecision(actualPrecisionDigits);
  os << "(";
  //6 extra characters: +/-, ., e, +/-, 2.
  os << std::setw(actualPrecisionDigits+6) << coeff.real();
  os << " ";
  os << std::setw(actualPrecisionDigits+6) << coeff.imag();
  os << ")";
  os << std::resetiosflags(std::ios_base::showpos);
}

Coeff seReadCoeff(std::istream& is) {
  seReadOBrace(is);
  Real re;
  if (!is)
    throw PolynomialIORanOutOfInputError();
  is >> re;
  is >> std::ws;
  Real im;
  if (!is)
    throw PolynomialIORanOutOfInputError();
  is >> im;
  seReadCBrace(is);
  return Coeff(re, im);
}

void seWriteMonomial(std::ostream& os, const Powers& powers, const Coeff& coeff) {
  os << "(";
  seWritePowers(os, powers);
  os << " ";
  seWriteCoeff(os, coeff);
  os << ")\n";
}

Polynomial& seReadMonomial(std::istream& is, Polynomial& result) {
  seReadOBrace(is);
  const Powers powers(seReadPowers(is, result.getNumVars()));
  const Coeff coeff(seReadCoeff(is));
  seReadCBrace(is);
  result.setMonomial(powers, coeff);
  return result;
}

void seWritePolynomial(std::ostream& os, const Polynomial& poly) {
  os << "(polynomial\n";
  os << "(num-variables " << poly.getNumVars() << ")\n";
  os << "(num-monomials " << poly.getNumTerms() << ")\n";
  os << "(powers-format \"dense\")\n";
  os << "(monomials\n";
  for (PowersToCoeffMap::const_iterator term = poly.getPowersAndCoeffs().begin();
       term != poly.getPowersAndCoeffs().end();
       ++term) {
    seWriteMonomial(os, term->first, term->second);
  }
  os << ")\n)\n";
}

Polynomial seReadPolynomial(std::istream& is) {
  seReadStart(is, "polynomial");
  //
  seReadStart(is, "num-variables");
  if (!is)
    throw PolynomialIORanOutOfInputError();
  Index numVars;
  is >> numVars;
  if (!(numVars > 0))
    throw PolynomialIOValueError();
  seReadEnd(is, "num-variables");
  //
  seReadStart(is, "num-monomials");
  if (!is)
    throw PolynomialIORanOutOfInputError();
  size_t numMonomials;
  is >> numMonomials;
  if (!(numMonomials >= 0))
    throw PolynomialIOValueError();
  seReadEnd(is, "num-monomials");
  //
  seReadStart(is, "powers-format");
  std::string powersFormat(seReadString(is));
  seReadEnd(is, "powers-format");
  //
  assert(powersFormat == "dense");
  //
  seReadStart(is, "monomials");
  Polynomial poly(numVars);
  for (size_t i = 0; i < numMonomials; ++i) {
    seReadMonomial(is, poly);
  }
  seReadEnd(is, "monomials");
  //
  seReadEnd(is, "polynomial");
  return poly;
}

void seWriteVectorOfPolynomials(std::ostream& os,
                                const std::vector< Polynomial >& polys) {
  os << "(vector-of-polynomials\n";
  os << "(num-polynomials " << polys.size() << ")\n";
  os << "(polynomials\n";
  for (std::vector< Polynomial >::const_iterator poly = polys.begin();
       poly != polys.end();
       ++poly) {
    seWritePolynomial(os, *poly);
  }
  os << ")\n)\n";
}

void seReadVectorOfPolynomials(std::istream& is,
                               std::vector< Polynomial >& result) {
  seReadStart(is, "vector-of-polynomials");
  //
  seReadStart(is, "num-polynomials");
  if (!is)
    throw PolynomialIORanOutOfInputError();
  Index numPolys;
  is >> numPolys;
  if (!(numPolys > 0))
    throw PolynomialIOValueError();
  seReadEnd(is, "num-polynomials");
  //
  seReadStart(is, "polynomials");
  result.clear();
  assert(result.empty());
  for (size_t i = 0; i < numPolys; ++i) {
    result.push_back(seReadPolynomial(is));
  }
  assert(result.size() == numPolys);
  seReadEnd(is, "polynomials");
  //
  seReadEnd(is, "vector-of-polynomials");
}
