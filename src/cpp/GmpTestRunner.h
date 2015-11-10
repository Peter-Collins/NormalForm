#ifndef GMP_TEST_RUNNER_H
#define GMP_TEST_RUNNER_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: GmpTestRunner implementation
//
// PURPOSE:
//
// Aggregate basic tests on GNU Multi-Precision functionality.
//
// NOTES:
//
// GMP already has its own test suites.  In particular, after
// compilation of GMP from sources, it is _ESSENTIAL_ to run "make
// check" in the distribution directory.  This will ensure that there
// have been no compiler errors.
//
// The main purpose of the tests here is to check correct inclusion
// and linking to GMP from client code.
//
//----------------------------------------------------------------------

//standard headers
#include <cstdlib>

//library headers
#include <cppunit/ui/text/TestRunner.h>

//project headers
#include "GmpTestRunner.h"
#include "GmpTest.h"

//function declarations
int main( int argc, char **argv);

#endif //GMP_TEST_RUNNER_H
