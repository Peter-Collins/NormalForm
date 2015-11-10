#!/usr/local/bin/python2.3

# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
from PackedTuple import *

max_bits = 30
max_term = (1<<max_bits)-1

class Examples(unittest.TestCase):

    examples = ((),
                (1,),
                (2,1),
                tuple(range(3)),
                tuple(range(10)),
                (1,2,11),
                tuple(range((1<<4)-1)))

    def test_len(self):
        bits = 4
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            self.assertEquals(len(packed), len(eg))

    def test_to_list(self):
        bits = 5
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            self.assertEquals(packed.to_list(), list(eg))

    def test_to_tuple(self):
        bits = 5
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            self.assertEquals(packed.to_tuple(), tuple(eg))

    def test_iter(self):
        bits = 5
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            t = tuple([x for x in packed])
            self.assertEquals(t, eg)

    def test_str(self):
        bits = 4
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            self.assertEquals(str(packed), str(eg))

    def test_repr(self):
        bits = 4
        for eg in self.examples:
            packed = PackedTuple(eg, bits=bits)
            self.assertEquals(repr(packed), 'PackedTuple(%s, bits=%d)'%(repr(eg),
                                                                        bits))

    def test_get_item(self):
        """__getitem__ should return relevant values."""
        for eg in self.examples:
            packed = PackedTuple(eg)
            for i, x in enumerate(eg):
                self.assertEquals(x, packed[i])

    def test_get_slice(self):
        for eg in self.examples:
            packed = PackedTuple(eg)
            for i in range(len(eg)):
                for j in range(i, len(eg)):
                    s = packed[i:j]
                    for k, x in enumerate(eg[i:j]):
                        self.assertEquals(s[k], x)

    def test_set_item_same(self):
        """__setitem__ with same item gives no change."""
        for eg in self.examples:
            packed = PackedTuple([0 for x in eg])
            for i, x in enumerate(eg):
                packed[i] = x
                self.assertEquals(x, packed[i])

    def test_set_item_zero(self):
        """__setitem__ with zero gives zero."""
        for eg in self.examples:
            packed = PackedTuple(eg)
            for i, x in enumerate(eg):
                packed[i] = 0
                self.assertEquals(0, packed[i])

    def test_hash(self):
        for eg in self.examples:
            packed = PackedTuple(eg)
            self.assertEquals(hash(packed), hash(eg))

class TupleLike(unittest.TestCase):

    """It is important that packed tuples behave like tuples."""

    def test_eq_to_tuple(self):
        t = (0, 1, 2, 3, 4)
        p = PackedTuple(t, bits=4)
        self.assert_(t == p)

    def test_as_dict_key(self):
        t = (0, 1, 2, 3, 4)
        p = PackedTuple(t, bits=4)
        d = {}
        d[p] = -99.0
        self.assert_(d.has_key(t))
        self.assert_(d[t] == d[p])
        self.assert_(d[t] == -99.0)

class Bits(unittest.TestCase):

    def test_value(self):
        for bits in xrange(1, max_bits):
            p = PackedTuple((), bits=bits)
            m = p.max_size()
            p = PackedTuple((m,), bits=bits)
            self.assertEquals(p._bits, bits)

class EmptyTuple(unittest.TestCase):

    def test_no_zeroth_element(self):
        for i in range(1, max_bits+1):
            empt = PackedTuple((), bits=i)
            self.assertRaises(IndexError, empt.__getitem__, 0)

    def test_len(self):
        for i in range(1, max_bits+1):
            empt = PackedTuple((), bits=i)
            self.assertEquals(len(empt), len(()))

    def test_hash(self):
        for i in range(1, max_bits+1):
            empt = PackedTuple((), bits=i)
            self.assertEquals(hash(empt), hash(()))

class BadTuple(unittest.TestCase):

    examples = ((1,2,3,-1),
                (-3,),
                (-1,10,1001))

    def test_negative_value(self):
        """Ensure that a tuple with negative value throws exception."""
        for eg in self.examples:
            try:
                packed = PackedTuple(eg)
            except NegativeIntegerError:
                pass
            else:
                self.assert_(0, 'Should throw exception at negative item.')

    def test_too_large_for_bits(self):
        """Ensure that exceeding the number of bits throws error."""
        self.assertRaises(InsufficientBitsError, PackedTuple, (max_bits+1,))

    def test_negative_set_item(self):
        """Ensure that we cannot try to set a negative value."""
        p = PackedTuple(range(10))
        try:
            p[0] = -1
        except NegativeIntegerError:
            pass
        else:
            self.assert_(0, 'Cannot set a negative item.')

class BadBits(unittest.TestCase):

    def test_non_positive_bits(self):
        for i in range(17):
            try:
                p = PackedTuple((1,2,3), bits=-i)
            except NonPositiveBitsError:
                pass
            else:
                self.assert_(0, 'Bits must be strictly positive.')

    def test_huge_bits(self):
        try:
            p = PackedTuple(range(500000), bits=100)
        except TooManyBitsError:
            pass
        else:
            self.assert_(0, 'Maximum 31 bits with this implementation.')

class MaxSize(unittest.TestCase):

    def test_value(self):
        for bits in range(1, max_bits+1):
            p = PackedTuple((), bits=bits)
            self.assertEquals(p.max_size(), (1<<bits)-1)

    def test_store(self):
        for bits in range(1, max_bits+1):
            p = PackedTuple((), bits=bits)
            m = p.max_size()
            t = tuple([i%(m+1) for i in range(66)])
            p = PackedTuple(t, bits=bits)
            for i in range(len(p)):
                self.assertEquals(p[i], i%(m+1), '%d'%i)

#test_len

def suite():
    suites = []
    suites.append(unittest.makeSuite(TupleLike))
    suites.append(unittest.makeSuite(Examples))
    suites.append(unittest.makeSuite(EmptyTuple))
    suites.append(unittest.makeSuite(BadTuple))
    suites.append(unittest.makeSuite(BadBits))
    suites.append(unittest.makeSuite(MaxSize))
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
