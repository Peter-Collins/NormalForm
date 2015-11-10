# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from Powers import Powers
from Polynomial import Polynomial
import Taylor

class Sine(unittest.TestCase):

    def test_terms(self):
        s = Taylor.Sine(1, 0)
        self.assert_(not s[0])
        self.assert_(not s[2])
        self.assert_(not s[4])
        self.assert_(s[1] == Polynomial(1, terms={Powers((1,)): +1.0}))
        self.assert_(s[3] == Polynomial(1, terms={Powers((3,)): -1.0/6.0}), s[3])
        self.assert_(s[5] == Polynomial(1, terms={Powers((5,)): +1.0/120.0}), s[5])

    def test_terms_embed(self):
        s = Taylor.Sine(2, 1)
        self.assert_(not s[0])
        self.assert_(not s[2])
        self.assert_(not s[4])
        self.assert_(s[1] == Polynomial(2, terms={Powers((0, 1)): +1.0}))
        self.assert_(s[3] == Polynomial(2, terms={Powers((0, 3)): -1.0/6.0}), s[3])
        self.assert_(s[5] == Polynomial(2, terms={Powers((0, 5)): +1.0/120.0}), s[5])

    def test_terms_cached(self):
        s = Taylor.Cached(Taylor.Sine(2, 1))
        self.assert_(not s[0])
        self.assert_(not s[2])
        self.assert_(not s[4])
        self.assert_(s[1] == Polynomial(2, terms={Powers((0, 1)): +1.0}))
        self.assert_(s[3] == Polynomial(2, terms={Powers((0, 3)): -1.0/6.0}), s[3])
        self.assert_(s[5] == Polynomial(2, terms={Powers((0, 5)): +1.0/120.0}), s[5])

class Cosine(unittest.TestCase):

    def test_terms(self):
        s = Taylor.Cosine(1, 0)
        self.assert_(not s[1])
        self.assert_(not s[3])
        self.assert_(not s[5])
        self.assert_(s[0] == Polynomial(1, terms={Powers((0,)): +1.0}))
        self.assert_(s[2] == Polynomial(1, terms={Powers((2,)): -1.0/2.0}), s[2])
        self.assert_(s[4] == Polynomial(1, terms={Powers((4,)): +1.0/24.0}), s[4])

    def test_terms_embed(self):
        s = Taylor.Cosine(2, 1)
        self.assert_(not s[1])
        self.assert_(not s[3])
        self.assert_(not s[5])
        self.assert_(s[0] == Polynomial(2, terms={Powers((0, 0)): +1.0}))
        self.assert_(s[2] == Polynomial(2, terms={Powers((0, 2)): -1.0/2.0}), s[2])
        self.assert_(s[4] == Polynomial(2, terms={Powers((0, 4)): +1.0/24.0}), s[4])

class Sum(unittest.TestCase):

    def test_sine_plus_cosine(self):
        s = Taylor.Cached(Taylor.Sine(2, 0))
        c = Taylor.Cached(Taylor.Cosine(2, 1))
        r = s+c
        self.assert_(r[0] == Polynomial(2, terms={Powers((0, 0)): +1.0}), r[0])
        self.assert_(r[1] == Polynomial(2, terms={Powers((1, 0)): +1.0}), r[1])
        self.assert_(r[2] == Polynomial(2, terms={Powers((0, 2)): -1.0/2.0}), r[2])
        self.assert_(r[3] == Polynomial(2, terms={Powers((3, 0)): -1.0/6.0}), r[3])
        self.assert_(r[4] == Polynomial(2, terms={Powers((0, 4)): +1.0/24.0}), r[4])
        self.assert_(r[5] == Polynomial(2, terms={Powers((5, 0)): +1.0/120.0}), r[5])

class Product(unittest.TestCase):

    def test_sine_times_cosine(self):
        s = Taylor.Cached(Taylor.Sine(2, 0))
        c = Taylor.Cached(Taylor.Cosine(2, 1))
        r = s*c
        self.assert_(not r[0])
        self.assert_(r[1] == Polynomial(2, terms={Powers((1, 0)): +1.0}), r[1])

    def test_sine_times_sine(self):
        s = Taylor.Cached(Taylor.Sine(2, 0))
        r = s*s
        self.assert_(not r[0])
        self.assert_(not r[1])
        self.assert_(r[2] == Polynomial(2, terms={Powers((2, 0)): +1.0}))
        self.assert_(not r[3])
        self.assert_(r[4] == Polynomial(2, terms={Powers((4, 0)): 2.0*(-1.0/6.0)}), r[4])
        self.assert_(not r[5])
        self.assert_(r[6] == Polynomial(2, terms={Powers((6, 0)): +2.0*1.0/120.0+1.0/36.0}), r[6])

class Bernoulli(unittest.TestCase):

    def test_values(self):
        b = Taylor.bernoulli
        self.assertEquals(b(0), +1.0)
        self.assertEquals(b(1), -1.0/2.0)
        self.assertEquals(b(2), +1.0/6.0)
        self.assertEquals(b(3), +0.0)
        self.assertEquals(b(4), -1.0/30.0)
        self.assertEquals(b(5), +0.0)
        self.assertEquals(b(6), +1.0/42.0)
        self.assertEquals(b(7), +0.0)
        self.assertEquals(b(8), -1.0/30.0)
        self.assertEquals(b(9), +0.0)

class Tanh(unittest.TestCase):

    def test_terms(self):
        t = Taylor.Tanh(1, 0)
        self.assert_(not t[0])
        self.assert_(t[1] == Polynomial(1, terms={Powers((1,)): 1.0}), t[1])

def suite():
    suites = []
    suites.append(unittest.makeSuite(Bernoulli))
    suites.append(unittest.makeSuite(Tanh))
    suites.append(unittest.makeSuite(Product))
    suites.append(unittest.makeSuite(Sum))
    suites.append(unittest.makeSuite(Sine))
    suites.append(unittest.makeSuite(Cosine))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
