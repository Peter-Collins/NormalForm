//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
//         Peter Collins 2006
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: MakeNormalForm
//
// PURPOSE:
//
// Import diagonalisation, Taylor series files from Python version and
// normalise, etc.
//
// NOTES:
//
// This is taken from SystemBathExample to run on any system
// Pending the C++ version of the diagonalisation code, this module is
// designed to import the Taylor series for the system-bath model from
// the Python code, along with the pre-computed vector of
// diagonalisation polynomials so that the diagonalisation can be
// applied via the C++ code using higher numerical precision, and so
// that we can proceed to perform normalisation, compute the
// coordinate changes, and extract the integrals.
//
// It is important to note that the files imported from Python will be
// at double precision, so one must take care that Polynomial objects
// created within the C++ part of the code have the correct GMP
// precision (higher than double).
//
//----------------------------------------------------------------------

//standard headers

//library headers

//project headers

//function declarations

int main(const int argc, const char* argv[]);
