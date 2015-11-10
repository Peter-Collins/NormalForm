//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: IOTest
//
// PURPOSE:
//
// A few experimental routines for parsing of simple s-expressions.
//
// NOTES:
//
// This was designed to get polynomial reading from the python side of
// the code up and running quickly.  This is quite an ugly piece of code
// and it would be a good idea to replace it as soon as time pressures
// allow.
//
// It is important to use enough default precision to read a coefficient
// from a file.  How to unify this is not clear at the moment.  For now,
// we just assume that the client code knows what it is doing and sets
// the required precision.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

//library headers
#include <gmpxx.h>

//project headers
#include "SExprIO.h"

int main(const int argc, const char* argv[]) {
  std::ifstream inFile("SExprIOTest.in", std::ios::in);
  assert(inFile);
  seReadStart(inFile, "hello");
  std::string format(seReadString(inFile));
  std::cerr << format << std::endl;
  std::vector< double > v(seReadVector< double >(inFile, 3));
  for (std::vector< double >::const_iterator p = v.begin(); p != v.end(); ++p) {
    std::cerr << *p << std::endl;
  }
  seReadStart(inFile, "this");
  mpf_set_default_prec(2500);
  mpf_class d;
  inFile >> d;
  std::cerr << std::setprecision(500);
  std::cerr << std::setiosflags(std::ios::scientific);
  std::cerr << d << std::endl;
  seReadEnd(inFile, "this");
  seReadEnd(inFile, "hello");
  inFile.close();
  return EXIT_SUCCESS;
}

