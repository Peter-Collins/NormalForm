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
from PolynomialRingIO import PolynomialRingIO
#from numpy.oldnumeric.mlab import fabs
#from MLab import fabs

from PolynomialTest import rand_poly
from PolynomialTest import _eq as _poly_eq

n_cases = 4

class PolynomialRingIOTest(unittest.TestCase):

    def test_read(self):
        ring = PolynomialRing(3)
        poly = Polynomial(3, {Powers((0, 0, 0)): +1,
                              Powers((1, 1, 0)): -2,
                              Powers((0, 2, 0)): -5+3J,
                              Powers((0, 1, 5)): +9J})
        file = StringIO.StringIO('(polynomial\n(num-variables 3)\n(num-monomials 4)\n(powers-format \"dense\")\n(monomials\n((0 0 0) (+1 +0))\n((1 1 0) (-2 +0))\n((0 2 0) (-5 +3))\n((0 1 5) (+0 +9))\n)\n)\n')
        read = PolynomialRingIO(ring).read_sexp_polynomial(file)
        file.close()
        self.assert_(poly == read, read)

    def test_write(self):
        ring = PolynomialRing(3)
        poly = Polynomial(3, {Powers((0, 0, 0)): +1,
                              Powers((1, 1, 0)): -2,
                              Powers((0, 2, 0)): -5+3J,
                              Powers((0, 1, 5)): +9J})
        file = StringIO.StringIO()
        PolynomialRingIO(ring).write_sexp_polynomial(file, poly)
        fstr = file.getvalue()
        file.close()
        file = StringIO.StringIO(fstr)
        read = PolynomialRingIO(ring).read_sexp_polynomial(file)
        file.close()
        self.assert_(poly == read, read)

    def test_read_zero(self):
        ring = PolynomialRing(3)
        poly = Polynomial(3)
        file = StringIO.StringIO('(polynomial\n(num-variables 3)\n(num-monomials 0)\n(powers-format \"dense\")\n(monomials\n)\n)\n')
        read = PolynomialRingIO(ring).read_sexp_polynomial(file)
        file.close()
        self.assert_(poly == read, read)

    def test_write_zero(self):
        ring = PolynomialRing(3)
        poly = Polynomial(3)
        file = StringIO.StringIO()
        PolynomialRingIO(ring).write_sexp_polynomial(file, poly)
        fstr = file.getvalue()
        file.close()
        file = StringIO.StringIO(fstr)
        read = PolynomialRingIO(ring).read_sexp_polynomial(file)
        file.close()
        self.assert_(poly == read, read)

    def test_read_vector_zero(self):
        ring = PolynomialRing(3)
        poly = Polynomial(3)
        polys = [poly]
        file = StringIO.StringIO('(vector-of-polynomials\n(num-polynomials 1)\n(polynomials\n(polynomial\n(num-variables 3)\n(num-monomials 0)\n(powers-format \"dense\")\n(monomials\n)\n)\n)\n)\n')
        read = PolynomialRingIO(ring).read_sexp_vector_of_polynomials(file)
        file.close()
        self.assert_(polys == read, read)

    def test_write_vector(self):
        ring = PolynomialRing(3)
        poly0 = Polynomial(3, {Powers((0, 0, 0)): +1,
                               Powers((1, 1, 0)): -2,
                               Powers((0, 2, 0)): -5+3J,
                               Powers((0, 1, 5)): +9J})
        poly1 = ring.zero()
        poly2 = ring.coordinate_monomial(1)
        polys = [poly0, poly1, poly2]
        file = StringIO.StringIO()
        PolynomialRingIO(ring).write_sexp_vector_of_polynomials(file, polys)
        fstr = file.getvalue()
        file.close()
        file = StringIO.StringIO(fstr)
        read = PolynomialRingIO(ring).read_sexp_vector_of_polynomials(file)
        file.close()
        self.assert_(polys == read, read)

def suite():
    suites = []
    suites.append(unittest.makeSuite(PolynomialRingIOTest))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
