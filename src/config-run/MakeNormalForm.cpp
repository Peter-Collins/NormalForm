//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
//         Peter Collins 2006
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: MakeNormalForm implementation.
//
//----------------------------------------------------------------------

//standard headers
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

//library headers

//project headers
#include "Types.h"
#include "DefaultNumericalPrecision.h"
#include "Polynomial.h"
#include "PolynomialIO.h"
#include "ClassicalLieAlgebra.h"
#include "SemiClassicalLieAlgebra.h"
#include "DepritTriangle.h"
#include "Utility.h"
#include "CoordinateChange.h"
#include "MakeNormalForm.h"

//function implementations

//--------------------------------------------------------------------
//
// Substitute a vector of polynomials into another vector of polynomials.
//
//--------------------------------------------------------------------
std::vector< Polynomial > polyvecSubstitute(const std::vector< Polynomial > pv1, const std::vector< Polynomial > pv2) {

    std::vector< Polynomial > polyvec;
    for (std::vector< Polynomial >::const_iterator iterpv1 = pv1.begin();
	 iterpv1 != pv1.end();
	 ++iterpv1) {
	polyvec.push_back((*iterpv1)(pv2));
    }
    return polyvec;
}

int main(const int argc, const char* argv[]) {

  if (argc!=4) {
    std::cerr << "argv[0] systemPrefix n_vars maxGrade\n";
    exit(1);
  }

  std::string systemPrefix;
  Index numDof;
  Index n_vars;
  Power maxGrade;
  //parameters for the computations
  systemPrefix = argv[1];
  n_vars = atoi(argv[2]);
  maxGrade = atoi(argv[3]);
#ifdef USE_GMP
  const size_t numBitsPrecision(96);
#else //USE_GMP
  const size_t numBitsPrecision(53);
#endif //USE_GMP
  
  //create context
  const DefaultNumericalPrecision prec(numBitsPrecision);
  numDof = n_vars / 2;
  //Decide if we are Classical or SemiClassical
  const LieAlgebraBase * thealgebra;
  // These should be inside the 'if' but g++ looses it
  const ClassicalLieAlgebra algebra1(numDof);
  const SemiClassicalLieAlgebra algebra2(numDof);
  if ( n_vars % 2 == 0 ) {
    thealgebra = &algebra1;
  } else {
    thealgebra = &algebra2;
  }

  //report on context
  std::cerr << "System          : " << systemPrefix << "\n";
  std::cerr << "numDof          : " << numDof << "\n";
  std::cerr << "numVars         : " << thealgebra->getNumVars() << "\n";
  std::cerr << "maxGrade        : " << maxGrade << "\n";
  std::cerr << "numBitsPrecision: " << numBitsPrecision << "\n";
  std::cerr << std::endl;

  //create file names and prefixes
  std::stringstream dofStr;
  dofStr << numDof;
  const std::string hFromEPrefix(systemPrefix + "--h-from-e");
  const std::string eFromDFileName(systemPrefix + "--e-from-d.vpol");
  const std::string eFromDCFileName(systemPrefix + "--e-from-dc.vpol");
  const std::string rFromCFileName(systemPrefix + "--r-from-c.vpol");
  const std::string cFromRFileName(systemPrefix + "--c-from-r.vpol");
  const std::string hFromDCPrefix(systemPrefix + "--h-from-dc");
  const std::string kFromNCPrefix(systemPrefix + "--k-from-nc");
  const std::string wPrefix(systemPrefix + "--w");

  //load diagonalisation polynomials
  std::ifstream eFromDRStr(eFromDFileName.c_str());
  assert(eFromDRStr.is_open());
  std::vector< Polynomial > eFromDR;
  seReadVectorOfPolynomials(eFromDRStr, eFromDR);
  eFromDRStr.close();

  //sanity check: correct number of polynomials
  assert(eFromDR.size() / 2 == numDof);
  assert(eFromDR.size() == thealgebra->getNumVars());

  //sanity check: for each polynomial
  for (std::vector< Polynomial >::const_iterator iterPoly = eFromDR.begin();
       iterPoly != eFromDR.end();
       ++iterPoly) {
    //sanity check: all polynomials belong to the algebra
    assert(thealgebra->hasElt(*iterPoly));
    //sanity check: all polynomials are non-zero
    assert(!(iterPoly->isZero()));
    //sanity check: all polynomials are grade 1 only
    // SemiClassical transforms get degree 2 from h_bar
    //assert(thealgebra->isIsoGrade(*iterPoly, 1));
    //sanity check: polys have correct precision, regardless of input
    if (iterPoly->getNumTerms() > 0) {
      PowersToCoeffMap::const_iterator iterTerm;
      iterTerm = iterPoly->getPowersAndCoeffs().begin();
      assert(iterTerm != iterPoly->getPowersAndCoeffs().end());
#ifdef USE_GMP
      // doubles don't have .get_prec(), maybe should have it in polys
      assert((iterTerm->second).real().get_prec() >= numBitsPrecision);
      assert((iterTerm->second).imag().get_prec() >= numBitsPrecision);
#endif //USE_GMP
    }
  }

  //--------------------------------------------------------------------
  //
  // temporary: zero constant term on taylor series (special)
  //
  //--------------------------------------------------------------------

  std::ifstream hFromE0Str((hFromEPrefix + "--grade-0.pol").c_str());
  assert(hFromE0Str.is_open());
  const Polynomial hFromE0(seReadPolynomial(hFromE0Str));
  hFromE0Str.close();
  assert(thealgebra->hasElt(hFromE0));
  assert(hFromE0.isZero());

  //--------------------------------------------------------------------
  //
  // temporary: zero linear term on taylor series (special)
  //
  //--------------------------------------------------------------------

  std::ifstream hFromE1Str((hFromEPrefix + "--grade-1.pol").c_str());
  assert(hFromE1Str.is_open());
  const Polynomial hFromE1(seReadPolynomial(hFromE1Str));
  hFromE1Str.close();
  assert(thealgebra->hasElt(hFromE1));
  assert(hFromE1.isZero());

  //--------------------------------------------------------------------
  //
  // temporary: deal with quadratic part (special)
  //
  //--------------------------------------------------------------------

  std::ifstream hFromE2Str((hFromEPrefix + "--grade-2.pol").c_str());
  assert(hFromE2Str.is_open());
  const Polynomial hFromE2(seReadPolynomial(hFromE2Str));
  hFromE2Str.close();
  assert(thealgebra->hasElt(hFromE2));
  assert(!(hFromE2.isZero()));
  assert(thealgebra->isIsoGrade(hFromE2, 2));
  std::cerr << hFromE2 << std::endl;

  //--------------------------------------------------------------------
  //
  // NOTE: Because the complexification involves sqrt(2), it would be
  // better to import the equilibrium type and compute the maps here!
  //
  //--------------------------------------------------------------------

  //for now, we import the complexification maps
  std::cerr << "WARNING! for now, we use double precision complexifiers";
  std::cerr << std::endl;
  // rFromC
  std::ifstream rFromCStr(rFromCFileName.c_str());
  assert(rFromCStr.is_open());
  std::vector< Polynomial > rFromC;
  seReadVectorOfPolynomials(rFromCStr, rFromC);
  rFromCStr.close();
  assert(rFromC.size() == thealgebra->getNumVars());
  for (std::vector< Polynomial >::const_iterator iterPoly = rFromC.begin();
       iterPoly != rFromC.end();
       ++iterPoly) {
    assert(thealgebra->hasElt(*iterPoly));
    assert(!(iterPoly->isZero()));
  }
  // cFromR
  std::ifstream cFromRStr(cFromRFileName.c_str());
  assert(cFromRStr.is_open());
  std::vector< Polynomial > cFromR;
  seReadVectorOfPolynomials(cFromRStr, cFromR);
  cFromRStr.close();
  assert(cFromR.size() == thealgebra->getNumVars());
  for (std::vector< Polynomial >::const_iterator iterPoly = cFromR.begin();
       iterPoly != cFromR.end();
       ++iterPoly) {
    assert(thealgebra->hasElt(*iterPoly));
    assert(!(iterPoly->isZero()));
    // SemiClassical transforms get degree 2 from h_bar
    //assert(thealgebra->isIsoGrade(*iterPoly, 1));
  }

  //--------------------------------------------------------------------
  //
  // compose the complexification with the real-diagonal transform
  //
  //--------------------------------------------------------------------

  std::cerr << "Composing the real-diagonal transformation with complex...";
  std::cerr << std::endl;
  //
  std::vector< Polynomial > eFromDC;
  for (std::vector< Polynomial >::const_iterator iterEFromDR = eFromDR.begin();
       iterEFromDR != eFromDR.end();
       ++iterEFromDR) {
    eFromDC.push_back((*iterEFromDR)(rFromC));
  }
  assert(eFromDC.size() == thealgebra->getNumVars());
  std::cerr << "Writing E from DC..." << std::endl;
  std::ofstream eFromDCStr(eFromDCFileName.c_str());
  seWriteVectorOfPolynomials(eFromDCStr, eFromDC);

  //--------------------------------------------------------------------
  //
  // back to dealing with the quadratic part (special)
  //
  //--------------------------------------------------------------------

  //real-diagonalise the term by applying the transformation
  std::cerr << "Diagonalising the quadratic term..." << std::endl;
  const Polynomial hFromDC2(hFromE2(eFromDC));
  assert(thealgebra->isIsoGrade(hFromDC2, 2));

  //remove off-diagonal part and confirm complex diagonal form
  std::pair< Polynomial, Polynomial > dAndN(thealgebra->zero(), thealgebra->zero());
  thealgebra->diagonalAndNonDiagonal(hFromDC2, dAndN);
  assert(thealgebra->isDiagonal(dAndN.first));
  std::cerr << dAndN.first << std::endl;
  std::cerr << dAndN.second << std::endl;

  std::cerr << dAndN.first << std::endl;
  std::cerr << "Infinity norm of off-diagonal complex part: ";
  std::cerr << (dAndN.second.lInfinityNorm()) << std::endl;

  const Polynomial& hFromDC2Diagonal(dAndN.first);

  //--------------------------------------------------------------------
  //
  // Setup the Deprit triangle
  //
  // NOTE:
  //
  // The terms that the method computeNormalFormAndGenerator produces
  // are all _OUTER_ terms.  That is, they have been properly
  // factorially-weighted for direct addition to the resulting
  // polynomial.  This is in contrast to the _INNER_ terms K_{i} and
  // W_{i}, which are the terms inside the Deprit triangle.  Such
  // inner terms cannot be summed directly over i to make a
  // polynomial; they are solely for use in the Deprit triangle
  // algorithm.  Hori's method would avoid such technicalities.
  //
  //--------------------------------------------------------------------

  DepritTriangle normalisationTriangle(*thealgebra, hFromDC2Diagonal);

  // we want to keep the computation in pipeline form as much as
  // possible.  for now, we use a loop here.  ultimately, these parts
  // will reside on separate processors with communication?

  for (Power row = 0; row < ((maxGrade - 2) + 1); ++row) {
    //currentGrade is previous level of computation, we must add one!
    const Power thisGrade(normalisationTriangle.getCurrentGrade() + 1);
    //create placeholder terms
    Polynomial zero(thealgebra->zero());
    Polynomial kTerm(thealgebra->zero());
    Polynomial wTerm(thealgebra->zero());
    //encode the grade as a string for use in file names
    std::stringstream gradeName;
    gradeName << thisGrade;
    //handle the quadratic part specially; we diagonalised it already
    if (row == 0) {
      //file output
      std::cerr << "Writing H from DC..." << std::endl;
      const std::string hName(hFromDCPrefix+"--grade-"+gradeName.str()+".pol");
      std::ofstream hFromDCStr(hName.c_str());
      seWritePolynomial(hFromDCStr, hFromDC2Diagonal);
      hFromDCStr.close();
      //normalisation
      normalisationTriangle.computeNormalFormAndGenerator(hFromDC2Diagonal, kTerm, wTerm);
      //file output
      std::cerr << "Writing K from NC..." << std::endl;
      const std::string kName(kFromNCPrefix+"--grade-"+gradeName.str()+".pol");
      std::ofstream kFromNCStr(kName.c_str());
      seWritePolynomial(kFromNCStr, kTerm);
      kFromNCStr.close();
    }
    else {
      for (int m = 0; m < 72; ++m)
        std::cerr << "=";
      std::cerr << std::endl;
      std::cerr << "Diagonalising input terms..." << std::endl;
      //diagonalisation
      const std::string fName(hFromEPrefix+"--grade-"+gradeName.str()+".pol");
      std::ifstream hFromEGStr(fName.c_str());
      assert(hFromEGStr.is_open());
      const Polynomial hFromEG(seReadPolynomial(hFromEGStr));
      hFromEGStr.close();
      assert(thealgebra->hasElt(hFromEG));
      assert(thealgebra->isIsoGrade(hFromEG, thisGrade));
      std::cerr << "Terms: " << hFromEG.getNumTerms() << std::endl;
      const Polynomial hFromDCG(hFromEG(eFromDC));
      std::cerr << "DiagTerms: " << hFromDCG.getNumTerms() << std::endl;
      //file output
      std::cerr << "Writing H from DC..." << std::endl;
      const std::string hName(hFromDCPrefix+"--grade-"+gradeName.str()+".pol");
      std::ofstream hFromDCStr(hName.c_str());
      seWritePolynomial(hFromDCStr, hFromDCG);
      hFromDCStr.close();
      //normalisation
      normalisationTriangle.computeNormalFormAndGenerator(hFromDCG, kTerm, wTerm);
      //file output
      std::cerr << "Writing K from NC..." << std::endl;
      const std::string kName(kFromNCPrefix+"--grade-"+gradeName.str()+".pol");
      std::ofstream kFromNCStr(kName.c_str());
      seWritePolynomial(kFromNCStr, kTerm);
      kFromNCStr.close();
    }
    //file output
    std::cerr << "Writing W..." << std::endl;
    std::stringstream rowStr;
    rowStr << row;
    const std::string wName(wPrefix+"--inner-"+rowStr.str()+".pol");
    std::ofstream wStr(wName.c_str());
    const Polynomial& wI(normalisationTriangle.getInnerGeneratorTerms().at(row));
    seWritePolynomial(wStr, wI);
    wStr.close();
  }
  //--------------------------------------------------------------------
  //
  // Use the calclated generating function to make the coordinate
  //  transforms using CoordinateChange.* (PRC)
  //
  //--------------------------------------------------------------------


  //--------------------------------------------------------------------
  //
  // generate the complex coordinate transforms
  //
  //--------------------------------------------------------------------

  if ( n_vars % 2 == 0 ) {
      std::cerr << "Calculating coordinate transforms" << std::endl;
  } else {
      std::cerr << "Semi-Classical -no  coordinate transforms" << std::endl;
      return EXIT_SUCCESS;
  }
  const ClassicalLieAlgebra algebra(numDof);
  std::vector< Polynomial > wInnerTerms=
      normalisationTriangle.getInnerGeneratorTerms();
  //
  // we transform each coordinate, both forward and backward
  //
  std::vector< Polynomial > x;
  std::vector< Polynomial > yFromX;
  std::vector< Polynomial > xFromY;
  for (Index i = 0; i < algebra.getNumVars(); ++i) {
      std::cerr << ".";
      const Polynomial xI(algebra.coordinateMonomial(i));
      x.push_back(xI);
      CoordinateChange forwardI(algebra, wInnerTerms, false);
      CoordinateChange reverseI(algebra, wInnerTerms, true);
      Polynomial xIForward(algebra.zero());
      Polynomial xIReverse(algebra.zero());
      Polynomial result(algebra.zero());
      for (Power row = 0; row < ((maxGrade - 2) + 1); ++row) {
          if (row == 0)
	      forwardI.computeNextTerm(xI, result);
          else
	      forwardI.computeNextTerm(algebra.zero(), result);
          xIForward += result;
          if (row == 0)
	      reverseI.computeNextTerm(xI, result);
          else
	      reverseI.computeNextTerm(algebra.zero(), result);
          xIReverse += result;
      }
      yFromX.push_back(xIForward);
      xFromY.push_back(xIReverse);
  }
  //--------------------------------------------------------------------
  //
  // Convert from (a vector of polynomial maps:
  // from vector of complex coordinates to a complex coordinate vector).
  // to (a vector of polynomial maps: from vector of reals to real vector).
  //
  //--------------------------------------------------------------------

  std::cerr << std::endl << "Converting to real" << std::endl;
  // we use the imported complexification maps rFromC cFromR
  std::vector< Polynomial > yFromXreal;
  yFromXreal = polyvecSubstitute(rFromC, polyvecSubstitute(yFromX, cFromR));
  
  const std::string FtName(systemPrefix+"--norm_to_diag.vpol");
  std::ofstream FtStr(FtName.c_str());
  seWriteVectorOfPolynomials(FtStr, yFromXreal);
  FtStr.close();

  // we use the imported complexification maps rFromC cFromR
  std::vector< Polynomial > xFromYreal;
  xFromYreal = polyvecSubstitute(rFromC, polyvecSubstitute(xFromY, cFromR));

  const std::string RtName(systemPrefix+"--diag_to_norm.vpol");
  std::ofstream RtStr(RtName.c_str());
  seWriteVectorOfPolynomials(RtStr, xFromYreal);
  RtStr.close();
  std::cerr << "Completed coordinate transforms" << std::endl;


  //--------------------------------------------------------------------
  //
  // When memory is a problem, we need to look at:-
  //
  // [-2] extending vector< Polynomial >, for example, creates a copy
  //      of the whole structure each time.  therefore, in the generator
  //      (and probably elsewhere), one should either pre-allocate with
  //      dummy values (but then take care with size()) or, better, use
  //      a different structure such as ::list.  Read on the STL types
  //      to determine the best course of action.
  //
  // [-1] recursion will be using the stack.  remove usage or ensure no
  //      non-reference statics within each stack context.
  //
  // [00] use equilibrium type to build C++ complexifier maps;
  //      this will remove some spurious extra terms.
  //
  // [01] more pass-by-reference for results;
  //      avoids temp copying.
  //
  // [02] provide a Polynomial::truncateSmallCoeffs(tolerance) method.
  //      truncation of small terms,
  //      avoids rounding issues (not so important with GMP).
  //
  // [03] store the Deprit triangle elements on disk and only read them
  //      at the point where they are needed.  since poisson brackets
  //      are the time-consuming thing, the read/write reconstruction is
  //      most likely not significant overhead.
  //
  // [04] improve polynomial operations, especially poisson bracket,
  //      replacing them with dedicated reference-result based versions,
  //      especially avoiding the creation of temporary variables, and
  //      avoiding copying in the delivery of results.
  //
  // [05] reduce GMP precision (at first, it was 256; 96 seems okay),
  //      chiefly for speed, although there will be memory impact.
  //
  // [06] pass smart pointers (not so important?),
  //      avoids unnecessary copies.
  //
  // [07] make a system bath with less bath modes! (want this anyway),
  //      starts with a smaller system to begin with; binomial coeffs!
  //
  // [08] do the diagonalisation computation here (using NR< GMP >),
  //      to get maximum accuracy and hence avoid spurious terms.
  //
  // [09] change the power and index type to char (but this leads to
  //      read/write and range complexities).
  //
  // [10] move to a polynomial representation as a map from single
  //      integer to coefficient, where the integer is an index into
  //      a monomial basis table.  The possible problem here is that
  //      one would need an efficient two-way mapping.
  //
  // Consider buying the following references:-
  //
  // [a] Large-scale C++ software design, Lakos.
  //
  // [b] Modern C++ design, Alexandrescu. (I already have this.)
  //
  //--------------------------------------------------------------------

  //finish cleanly
  return EXIT_SUCCESS;
}

