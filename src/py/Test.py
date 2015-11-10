#!/usr/local/bin/python2.3

"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Test

PURPOSE:

Global test suite to test all modules.

NOTES:

This module imports all the individual unit tests and adds them to a
single master test suite.

I recommend using the pyunit GUI in order to run tests, although it is
just as easy to call this script from the comman line using pythonX.XX
Test.py.

"""

import unittest

import AccuracyTest
import ComplexTest
import CoordinateChangeTest
import DiagonalTest
import EigenSystemTest
import EigenValueVectorPairTest
import EquilibriumTypeTest
import IntegralExtractorTest
import IsogradeInnerTaylorCoeffsTest
import LieAlgebraTest
import LieTriangleTest
import NameLookupSemanticsTest
import PackedPowersTest
import PackedTupleTest
import PolynomialTest
import PolynomialRingTest
import PolynomialRingIOTest
import ShelvedDictTest
import SystemBathTest
import TaylorTest
import TuplePowersTest
import UtilityTest
import UuidTest
import GmpTest

def suite():
    suites = []

    #accumulate the individual test suites
    suites.append(AccuracyTest.suite())
    suites.append(ComplexTest.suite())
    suites.append(CoordinateChangeTest.suite())
    suites.append(DiagonalTest.suite())
    suites.append(EigenSystemTest.suite())
    suites.append(EigenValueVectorPairTest.suite())
    suites.append(EquilibriumTypeTest.suite())
    suites.append(IntegralExtractorTest.suite())
    suites.append(IsogradeInnerTaylorCoeffsTest.suite())
    suites.append(LieTriangleTest.suite())
    suites.append(NameLookupSemanticsTest.suite())
    suites.append(PackedPowersTest.suite())
    suites.append(PackedTupleTest.suite())
    suites.append(PolynomialTest.suite())
    suites.append(PolynomialRingTest.suite())
    suites.append(PolynomialRingIOTest.suite())
    suites.append(ShelvedDictTest.suite())
    suites.append(SystemBathTest.suite())
    suites.append(TaylorTest.suite())
    suites.append(TuplePowersTest.suite())
    suites.append(UtilityTest.suite())
    suites.append(UuidTest.suite())

    #this added at the end; the longest running test
    suites.append(LieAlgebraTest.suite())
    suites.append(GmpTest.suite())

    #return the composite suite
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')

