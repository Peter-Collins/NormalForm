# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
import random
from copy import deepcopy

from Polynomial import *
from PolynomialRing import PolynomialRing
from IsogradeInnerTaylorCoeffs import IsogradeInnerTaylorCoeffs
from Powers import Powers

n_cases = 2 #32

tmp_path = '../tmp/'

class ListToPoly(unittest.TestCase):

    def setUp(self):
        ring = PolynomialRing(6)
        self.iso = IsogradeInnerTaylorCoeffs(ring, offset=1)

    def test_empty_list(self):
        p_list = []
        self.assert_(self.iso.list_to_poly(p_list) == Polynomial(6))

    def test_empty_poly(self):
        p_poly = Polynomial(6)
        self.assert_(self.iso.poly_to_list(p_poly) == [])

    def test_inverse_eg1(self):
        p_poly = Polynomial(6, terms={Powers((1, 0, 0, 0, 0, 0)): 1.0})
        p_list = self.iso.poly_to_list(p_poly)
        self.assert_(p_list == [p_poly])
        q_poly = self.iso.list_to_poly(p_list)
        self.assert_(q_poly == p_poly, (q_poly, p_poly))

    def test_inverse_eg2(self):
        p1 = Polynomial(6, terms={Powers((1, 0, 0, 0, 0, 0)): 1.0})
        p3 = Polynomial(6, terms={Powers((0, 1, 1, 1, 0, 0)): 1.0J})
        z = Polynomial(6)
        p_poly = p1 + p3
        p_list = self.iso.poly_to_list(p_poly)
        self.assert_(p_list == [(1)*p1, (1*1)*z, (1*1*2)*p3], p_list)
        q_poly = self.iso.list_to_poly(p_list)
        self.assert_(q_poly == p_poly)


def suite():
    suites = []
    suites.append(unittest.makeSuite(ListToPoly))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
