#!/usr/local/bin/python2.3

# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
from TuplePowers import TuplePowers
from PackedPowers import *

from random import Random
gen = Random(6543210)

class Default(unittest.TestCase):

    """Test the properties of the default constructor."""

    def test_call(self):
        m = PackedPowers()
        self.assertRaises(TypeError, m.__call__, 1.0)

    def test_degree(self):
        m = PackedPowers()
        self.assertEquals(m.degree(), 1)

    def test_repr(self):
        m = PackedPowers()
        self.assertEquals(repr(m), 'PackedPowers((1,), bits=4)')

    def test_diff_coeff(self):
        m = PackedPowers()
        c, d = m.diff(0)
        self.assertEquals(c, 1)

    def test_diff_powers(self):
        m = PackedPowers()
        c, d = m.diff(0)
        self.assertEquals(repr(d), 'PackedPowers((0,), bits=4)')

    def test_pows(self):
        p = PackedPowers()
        for i in xrange(10):
            self.assertEquals(p**i, PackedPowers((i,)))

def rand_int(max_val):
    return gen.randrange(0, max_val)
    
def rand_int_tuple(length, max_val):
    return tuple([rand_int(max_val) for i in xrange(length)])
        
def min_bits_needed_to_represent_power(p):
    bits = 1
    while(p>0):
        p = p>>1
        bits += 1
    return bits

class PowersLike(unittest.TestCase):

    def test_empty(self):
        t = ()
        m = TuplePowers(t)
        p = PackedPowers(t, bits=4)
        self.assert_(t == p)
        self.assert_(m == p)

    def test_not_eq(self):
        t = ()
        p = PackedPowers(t, bits=4)
        self.assert_(not (TuplePowers((1,)) == p))

    def test_example(self):
        t = (1, 2, 3)
        m = TuplePowers(t)
        p = PackedPowers(t, bits=3)
        self.assert_(t == p)
        self.assert_(m == p)

    def test_dict(self):
        t = (1, 2, 3)
        m = TuplePowers(t)
        p = PackedPowers(t, bits=3)
        d = {}
        d[m] = -99.0
        self.assert_(d.has_key(p))
        self.assert_(d[p] == d[m])

class TupleLike(unittest.TestCase):

    def test_empty(self):
        t = ()
        p = PackedPowers(t, bits=4)
        self.assert_(t == p)

    def test_not_eq(self):
        t = ()
        p = PackedPowers(t, bits=4)
        self.assert_(not ((1,) == p))

    def test_example(self):
        t = (1, 2, 3)
        p = PackedPowers(t, bits=3)
        self.assert_(t == p)

    def test_dict(self):
        t = (1, 2, 3)
        p = PackedPowers(t, bits=3)
        d = {}
        d[t] = -99.0
        self.assert_(d.has_key(p))
        self.assert_(d[p] == d[t])

    def test_to_tuple(self):
        for i in xrange(50):
            t = rand_int_tuple(40, 14)
            m = PackedPowers(t, bits=4)
            tt = m.to_tuple()
            self.assert_(type(tt) == tuple)
            self.assert_(t == tt)

class BitsDoNotAffectEq(unittest.TestCase):

    def test_me(self):
        t = (1, 2, 3, 4)
        p0 = PackedPowers(t, bits=4)
        p1 = PackedPowers(t, bits=5)
        self.assert_(p0 == p1)
        
        
class Examples(unittest.TestCase):

    """Test many of the routines for a set of examples."""

    def setUp(self):
        #would need to encode bits in the examples?!?
        examples = [(),
                    (1,),
                    (2,),
                    (1,2,3),
                    (1,1,0,5)]
        for bits in xrange(5,10):
            max_val = (1<<bits)-1
            for length in xrange(1,10):
                for case in xrange(10):
                    t = rand_int_tuple(length, max_val)
                    examples.append(t)
        self.examples = []
        for eg in examples:
            if eg:
                max_pow = max(eg)
                bits = min_bits_needed_to_represent_power(max_pow)
                self.examples.append((eg, rand_int(3)+bits))
            else:
                self.examples.append((eg, rand_int(4)+3))

    def test_call_wrong_len(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            for l in xrange(0, len(eg)*2):
                x = rand_int_tuple(l, 10)
                if not l==len(eg):
                    self.assertRaises(IndexError, p.__call__, x)

    def test_call_non_seq(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertRaises(TypeError, p.__call__, 1.0)

    def test_call_all_zero(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            zero = (0.0,)*len(eg)
            if p.degree()==0:
                self.assertEquals(p(zero), 1.0)
            else:
                self.assertEquals(p(zero), 0.0)

    def test_call_all_one(self):
        """any monomial (no coeff) evaluated at (1.0,)*len gives 1.0"""
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            one = (1.0,)*len(eg)
            self.assertEquals(p(one), 1.0)

    def test_call_any_zero(self):
        """any monomial (no coeff) evaluated at (1.0,)*len gives 1.0"""
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            for var in range(len(eg)):
                x = [float(i+1) for i in rand_int_tuple(len(eg), 3)]
                x[var] = 0.0
                if p[var]==0:
                    self.assert_(not p(x)==0.0)
                else:
                    self.assertEquals(p(x), 0.0)

    def test_call_negative_counts(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            deg = sum(eg)
            neg = sum([1 for i in eg if i%2])%2
            eno = (-1.0,)*len(eg)
            r = p(eno)
            if deg==0:
                self.assert_(r==1.0)
            if neg:
                self.assert_(r==-1.0)
            else:
                self.assert_(r==1.0)

    def test_number_of_variables(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p.number_of_variables(), len(eg))

    def test_degree(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p.degree(), sum(eg))

    def test_repr(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(repr(p), 'PackedPowers(%s, bits=%d)'%(repr(eg),
                                                                      bits))

    def test_pow_zero_is_one(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p**0, PackedPowers((0,)*len(eg)))

    def test_pow_one_is_same(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p**1, p)

    def test_pow(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            for i in xrange(10):
                eg10 = tuple([i*j for j in eg])
                try:
                    self.assertEquals(p**i, PackedPowers(eg10))
                except InsufficientBitsError:
                    pass

    def test_mul(self):
        lengths = 50
        cases = 20
        bits = 25
        max_val = int((1<<bits)-1)/2
        for length in xrange(lengths):
            for case in xrange(cases):
                pt = rand_int_tuple(length, max_val)
                p = PackedPowers(pt, bits)
                qt = rand_int_tuple(length, max_val)
                q = PackedPowers(qt, bits)
                rt = tuple([a+b for a,b in zip(pt, qt)])
                r = PackedPowers(rt, bits)
                self.assertEquals(p*q, r)
                self.assertEquals(repr(p*q),
                                  'PackedPowers(%s, bits=%d)'%(repr(rt),
                                                                 bits))

    def test_equal_self(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p, p)

    def test_less_equal_self(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assert_(p<=p)

    def test_greater_equal_self(self):
        for eg in self.examples:
            p = TuplePowers(eg)
            self.assert_(p>=p)

    def test_not_greater_self(self):
        for eg in self.examples:
            p = TuplePowers(eg)
            self.assert_(not (p>p))

    def test_not_less_self(self):
        for eg in self.examples:
            p = TuplePowers(eg)
            self.assert_(not (p<p))

    def test_not_not_equal_self(self):
        for eg in self.examples:
            p = TuplePowers(eg)
            self.assert_(not (p!=p))

    def test_equal_other_same(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            self.assertEquals(p, PackedPowers(eg, bits))

    def test_diff(self):
        for eg, bits in self.examples:
            p = PackedPowers(eg, bits)
            for var in xrange(len(eg)):
                coeff, q = p.diff(var)
                self.assertEquals(coeff, eg[var])
                msg = repr(eg)+' -> '+repr(q)
                for j in xrange(len(eg)):
                    if j==var:
                        if eg[j]==0:
                            self.assertEquals(q[j], 0, msg)
                        else:
                            self.assertEquals(q[j], eg[j]-1, msg)
                    else:
                        self.assertEquals(q[j], eg[j], msg)

class Diff(unittest.TestCase):

    """Diff-specific range tests."""

    def test_negative_index(self):
        for i in xrange(-10, 0):
            m = PackedPowers()
            self.assertRaises(IndexError, m.diff, i)

    def test_bad_index(self):
        for i in xrange(1, 10):
            m = PackedPowers()
            self.assertRaises(IndexError, m.diff, i)

class Pow(unittest.TestCase):

    """Pow-specific range tests."""

    def test_negative_pow(self):
        for i in xrange(-10, 0):
            m = PackedPowers()
            self.assertRaises(ValueError, m.__pow__, i)

def suite():
    suites = []
    suites.append(unittest.makeSuite(BitsDoNotAffectEq))
    suites.append(unittest.makeSuite(PowersLike))
    suites.append(unittest.makeSuite(TupleLike))
    suites.append(unittest.makeSuite(Default))
    suites.append(unittest.makeSuite(Examples))
    suites.append(unittest.makeSuite(Diff))
    suites.append(unittest.makeSuite(Pow))
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
