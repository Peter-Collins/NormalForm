#!/usr/local/bin/python2.3

# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
from SparseTuple import *

class Examples(unittest.TestCase):

    examples = ((),
                (1,),
                (2,1),
                tuple(range(3)),
                tuple(range(10)),
                (1,2,11),
                tuple(range((1<<4)-1)))

    def test_not_eq_same_length_sparse(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            t = tuple(range(99, 99+len(eg)))
            self.assert_(len(s) == len(t), (s, t))
            if len(s) > 0:
                u = SparseTuple(t[:-1]+(999,))
                self.assert_(not s == u, (s, u))

    def test_not_eq_same_length_tuple(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            t = tuple(range(99, 99+len(eg)))
            self.assert_(len(s) == len(t), (s, t))
            if len(s) > 0:
                self.assert_(not s == t, (s, t))
            else:
                self.assert_(s == t, (s, t))

    def test_eq_tuple_self(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            self.assert_(s == eg, (s, eg))

    def test_eq_to_tuple_self(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            t = s.to_tuple()
            self.assert_(s == t, (s, t))

    def test_not_eq_other(self):
        not_in_examples = tuple(range(51))
        for eg in self.examples:
            s = SparseTuple(eg)
            self.assert_(not (s == not_in_examples), (s, not_in_examples))

    def test_eq_self(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            self.assert_(s == s, (s, s))

    def test_eq_copy(self):
        for eg in self.examples:
            s = SparseTuple(eg)
            v = SparseTuple(eg)
            self.assert_(s == v, (s, v))

    def _tuple_to_nonzero_dict(self, t):
        d = {}
        for i, e in enumerate(t):
            if e != 0:
                d[i] = e
        return d

    def test_both_constructors(self):
        for eg in self.examples:
            d = self._tuple_to_nonzero_dict(eg)
            p0 = SparseTuple(eg)
            p1 = SparseTuple(len(eg), d)
            self.assertEquals(p0, p1)

    def test_len(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            self.assertEquals(len(sparse), len(eg))

    def test_to_list(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            self.assertEquals(sparse.to_list(), list(eg))

    def test_to_tuple(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            self.assertEquals(sparse.to_tuple(), tuple(eg))

    def test_iter(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            t = tuple([x for x in sparse])
            self.assertEquals(t, eg)

    def test_str(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            self.assertEquals(str(sparse), str(eg))

    def test_repr(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            d = {}
            for i, e in enumerate(eg):
                if e != 0:
                    d[i] = e
            self.assertEquals(repr(sparse), 'SparseTuple(%d, %s)'%(len(eg),
                                                                   repr(d)))

    def test_construction(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            d = {}
            for i, e in enumerate(eg):
                if e != 0:
                    d[i] = e
            sparse2 = SparseTuple(len(eg), d)
            self.assertEquals(sparse, sparse2)
                
    def test_get_item(self):
        """__getitem__ should return relevant values."""
        for eg in self.examples:
            sparse = SparseTuple(eg)
            for i, x in enumerate(eg):
                self.assertEquals(x, sparse[i])

    def test_get_slice(self):
        for eg in self.examples:
            sparse = SparseTuple(eg)
            for i in range(len(eg)):
                for j in range(i, len(eg)):
                    s = sparse[i:j]
                    for k, x in enumerate(eg[i:j]):
                        self.assertEquals(s[k], x)

    def test_set_item_same(self):
        """__setitem__ with same item gives no change."""
        for eg in self.examples:
            sparse = SparseTuple([0 for x in eg])
            for i, x in enumerate(eg):
                sparse[i] = x
                self.assertEquals(x, sparse[i])

    def test_set_item_zero(self):
        """__setitem__ with zero gives zero."""
        for eg in self.examples:
            sparse = SparseTuple(eg)
            for i, x in enumerate(eg):
                sparse[i] = 0
                self.assertEquals(0, sparse[i])

    def test_hash(self):
        
        """Note that we only insist that sparse tuples can be used in
        the place of the equivalent tuple in a dictionary.  This does
        not mean that they have to have the same hash, only that the
        lookup works correctly (which uses eq as well as hash)."""

        t = (0, 1, 2, 3, 4)
        p = SparseTuple(t)
        self.assert_(hash(t) == hash(p))

class TupleLike(unittest.TestCase):

    """It is important that sparse tuples behave like tuples."""

    def test_eq_to_tuple(self):
        t = (0, 1, 2, 3, 4)
        p = SparseTuple(t)
        self.assert_(t == p)

    def test_as_dict_key(self):
        """

        This test depends on has values being the same.

        We have deprecated this requirement.

        """
        t = (0, 1, 2, 3, 4)
        p = SparseTuple(t)
        d = {}
        d[p] = -99.0
        self.assert_(d.has_key(t))
        self.assert_(d[t] == d[p])
        self.assert_(d[t] == -99.0)

class EmptyTuple(unittest.TestCase):

    def test_no_zeroth_element(self):
        empt = SparseTuple(())
        self.assertRaises(IndexError, empt.__getitem__, 0)

    def test_len(self):
        empt = SparseTuple(())
        self.assertEquals(len(empt), len(()))

    def test_hash(self):
        empt = SparseTuple(())
        self.assertEquals(hash(empt), hash(()))

class BadTuple(unittest.TestCase):

    examples = ((1,2,3,-1),
                (-3,),
                (-1,10,1001))

    def test_negative_value(self):
        """Ensure that a tuple with negative value throws exception."""
        for eg in self.examples:
            try:
                sparse = SparseTuple(eg)
            except NegativeIntegerError:
                pass
            else:
                self.assert_(0, 'Should throw exception at negative item.')

    def test_negative_set_item(self):
        """Ensure that we cannot try to set a negative value."""
        p = SparseTuple(range(10))
        try:
            p[0] = -1
        except NegativeIntegerError:
            pass
        else:
            self.assert_(0, 'Cannot set a negative item.')

def suite():
    suites = []
    suites.append(unittest.makeSuite(TupleLike))
    suites.append(unittest.makeSuite(Examples))
    suites.append(unittest.makeSuite(EmptyTuple))
    suites.append(unittest.makeSuite(BadTuple))
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
