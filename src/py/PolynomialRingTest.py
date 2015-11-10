"""

Author: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

"""

import unittest
import random
import StringIO

from Powers import Powers
from Polynomial import Polynomial
from PolynomialRing import PolynomialRing
from MLab import fabs

from PolynomialTest import rand_poly
from PolynomialTest import _eq as _poly_eq

def _eq(a, b, tol=1.0e-8):
    """

    We need to loosen tolerance for Moyal bracket calculations.

    """
    return _poly_eq(a, b, tol)

n_cases = 4

class PolynomialRingTest(unittest.TestCase):

    def test_one(self):
        for n_vars in xrange(1, 7):
            ring = PolynomialRing(n_vars)
            one = ring.one()
            self.assert_(isinstance(one, Polynomial))
            self.assert_(one.n_vars() == n_vars)
            self.assert_(len(one) == 1)
            po, co = list(one.powers_and_coefficients())[0]
            self.assert_(len(po) == n_vars)
            for p in po:
                self.assert_(p == 0)
            self.assert_(co == 1.0)

    def test_zero(self):
        for n_vars in xrange(1, 7):
            ring = PolynomialRing(n_vars)
            zero = ring.zero()
            self.assert_(isinstance(zero, Polynomial))
            self.assert_(zero.n_vars() == n_vars)
            self.assert_(len(zero) == 0)
            self.assert_(not zero)

    def test_n_vars(self):
        for n_vars in xrange(1, 25):
            ring = PolynomialRing(n_vars)
            self.assert_(ring.n_vars() == n_vars)

    def test_coordinate_monomial(self):
        for n_vars in xrange(1, 7):
            ring = PolynomialRing(n_vars)
            for var in xrange(n_vars):
                x = ring.coordinate_monomial(var)
                self.assert_(isinstance(x, Polynomial))
                self.assert_(x.n_vars() == n_vars)
                self.assert_(len(x) == 1)
                po, co = list(x.powers_and_coefficients())[0]
                self.assert_(len(po) == n_vars)
                for i, p in enumerate(po):
                    if i == var:
                        self.assert_(p == 1)
                    else:
                        self.assert_(p == 0)
                self.assert_(co == 1.0)

    def test_grad(self):
        poly = Polynomial(3, {Powers((0, 0, 0)): 1.0,
                              Powers((1, 1, 0)): -2.0,
                              Powers((0, 2, 0)): 5.3,
                              Powers((0, 1, 5)): -9})
        ring = PolynomialRing(3)
        grad = ring.grad(poly)
        g0 = Polynomial(3, {Powers((0, 1, 0)): -2.0})
        g1 = Polynomial(3, {Powers((1, 0, 0)): -2.0,
                            Powers((0, 1, 0)): 10.6,
                            Powers((0, 0, 5)): -9})
        g2 = Polynomial(3, {Powers((0, 1, 4)): -45.0})
        for expected, actual in zip((g0, g1, g2), grad):
            self.assert_(expected == actual, (expected, actual))

def suite():
    suites = []
    suites.append(unittest.makeSuite(PolynomialRingTest))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
