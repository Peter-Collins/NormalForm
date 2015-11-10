# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from EigenSystem import *

class SortGeneral(unittest.TestCase):

    def test_example(self):
        seq = [9, 8, 7, 6, 3, 4, 4, 5, 3, 2, 1, 0]
        sorted = sort_by_ranking_function(seq)
        self.assertEquals(sorted, [0, 1, 2, 3, 3, 4, 4, 5, 6, 7, 8, 9])

class SortByPosition(unittest.TestCase):

    def _first(self, x):
        return x[0]

    def _last(self, x):
        return x[-1]

    def setUp(self):
        self.seq = [(1, 'd', 1), (2, 'b', 10), (3, 'c', 4), (4, 'a', 2)]

    def test_first(self):
        sorted = sort_by_ranking_function(self.seq, self._first)
        self.assertEquals(sorted,
                          [(1, 'd', 1), (2, 'b', 10), (3, 'c', 4), (4, 'a', 2)])

    def test_last(self):
        sorted = sort_by_ranking_function(self.seq, self._last)
        self.assertEquals(sorted, 
                          [(1, 'd', 1), (4, 'a', 2), (3, 'c', 4), (2, 'b', 10)])

    def test_next_to_last(self):
        def next_to_last(x):
            return x[-2]
        sorted = sort_by_ranking_function(self.seq, next_to_last)
        self.assertEquals(sorted, 
                          [(4, 'a', 2), (2, 'b', 10), (3, 'c', 4), (1, 'd', 1)])

class SortByRealAndImag(unittest.TestCase):

    def test_order(self):
        seq = [1+2J, -3J, +3J, 1-2J, -5, -7J, 5]
        res = [5, -5, 1+2J, 1-2J, -7J, 3J, -3J]
        class dummy:
            def __init__(self, val, vec):
                self.val = val
                self.vec = vec
            def __eq__(self, other):
                return (self.val == other.val) and (self.vec == other.vec)
        def add_vec(val_seq):
            return [dummy(val, [0.0]) for val in val_seq]
        sorted = sort_by_ranking_function(add_vec(seq), rank_real_before_imag)
        self.assertEquals(sorted, add_vec(res))

class ValVec:

    """Dummy class to test eigenvalues and eigvectors."""

    def __init__(self, val, vec):
        self.val = val
        self.vec = vec

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return (self.val==other.val and self.vec==other.vec)

class LocatePairs(unittest.TestCase):

    def test_exists(self):
        seq = []
        eig = EigenSystem(seq, 1.0e-15)
        eig._collect_plus_minus_pairs()

    def test_len(self):
        seq = [+1, -1.001, -2.0003, +2.001, 5]
        eig = EigenSystem(seq, 1.0e-15)
        try:
            out = eig._collect_plus_minus_pairs()
        except IndexError:
            pass
        else:
            assert 0, 'odd number of entries should be detected.'

    def test_empty(self):
        seq = []
        eig = EigenSystem(seq, 1.0e-15)
        eig._collect_plus_minus_pairs()
        self.assertEquals(eig.get_plus_minus_pairs(), [])

    def test_example(self):
        seq = [ValVec(+1, 'v0'),
               ValVec(-1.001, 'v1'),
               ValVec(-2.0003, 'v2'),
               ValVec(+2.001, 'v3')]
        eig = EigenSystem(seq, 1.0e-1)
        eig._collect_plus_minus_pairs()
        exp = [(ValVec(+1, 'v0'), ValVec(-1.001, 'v1')),
               (ValVec(+2.001, 'v3'), ValVec(-2.0003, 'v2'))]
        self.assertEquals(eig.get_plus_minus_pairs(), exp)

    def test_exceed_tolerance(self):
        seq = [ValVec(+1, 0),
               ValVec(-1.1, 1),
               ValVec(-2.0003, 2),
               ValVec(+2.001, 3)]
        try:
            eig = EigenSystem(seq, 1.0e-1)
            eig._collect_plus_minus_pairs()
        except ValueError:
            pass
        else:
            assert 0, 'value exceeding tolerance should be detected.'

    def test_too_close_to_zero(self):
        seq = [ValVec(+1, 0),
               ValVec(-1.1, 1),
               ValVec(-0.0003, 2),
               ValVec(+0.001, 3)]
        try:
            eig = EigenSystem(seq, 1.0e-1)
            eig._collect_plus_minus_pairs()
        except ValueError:
            pass
        else:
            assert 0, 'values too close to zero should be detected.'

class Pure(unittest.TestCase):

    def test_zero_gives_both_false(self):
        c = 0.0
        self.assert_(not is_pure_real(c))
        self.assert_(not is_pure_imag(c))

    def test_nonzero_float_is_pure_real(self):
        c = -0.1
        self.assert_(is_pure_real(c))
        self.assert_(not is_pure_imag(c))

    def test_nonzero_imag_is_pure_imag(self):
        c = 0.3J
        self.assert_(not is_pure_real(c))
        self.assert_(is_pure_imag(c))

    def test_mixed_complex_gives_both_false(self):
        c = 54.3-0.3J
        self.assert_(not is_pure_real(c))
        self.assert_(not is_pure_imag(c))

class Type(unittest.TestCase):

    def test_determine_empty(self):
        seq = []
        eig = EigenSystem(seq, 1.0e-1)
        eig._collect_plus_minus_pairs()
        eig._determine_equilibrium_type()
        self.assertEquals(eig.get_equilibrium_type(), '')

    def test_determine_non_seq(self):
        seq = [1,2]
        try:
            eig = EigenSystem(seq, 1.0e-1)
            eig.plus_minus_pairs = seq
            eig._determine_equilibrium_type()
        except TypeError:
            pass
        else:
            assert 0

    def test_determine_non_pairs(self):
        seq = [(1,2,3)]
        try:
            eig = EigenSystem(seq, 1.0e-1)
            eig.plus_minus_pairs = seq
            eig._determine_equilibrium_type()
        except ValueError:
            pass
        else:
            assert 0

    def test_determine_example(self):
        seq = [ValVec(+1J, 'v0'),
               ValVec(-1.001J, 'v1'),
               ValVec(-2.0003, 'v2'),
               ValVec(+2.001, 'v3'),
               ValVec(-2.0003J, 'v4'),
               ValVec(+2.001J, 'v5')]
        eig = EigenSystem(seq, 1.0e-1)
        eig._collect_plus_minus_pairs()
        eig._determine_equilibrium_type()
        self.assertEquals(eig.get_equilibrium_type(), 'csc')

    def test_saddles_come_first(self):
        seq = [ValVec(+1J, (0.1,0.2,0.3)),
               ValVec(-1.001J, (0.4,0.5,0.6)),
               ValVec(-2.0003, (-0.1,0.9,-0.4)),
               ValVec(+2.001, (0.8,0.7,0.6)),
               ValVec(-2.0003J, (0.2,0.3,0.4)),
               ValVec(+2.001J, (-0.2,0.3,-0.4))]
        eig_tolerance = 1.0e-3
        eig = EigenSystem(seq, eig_tolerance)
        eig._truncate_eigen_values_and_vectors()
        eig._sort_eigensystem_reals_first()
        eig._collect_plus_minus_pairs()
        eig._determine_equilibrium_type()
        eq_type = eig.get_equilibrium_type()
        had_centre_yet = False
        for saddle in eq_type:
            if saddle == 's':
                assert not had_centre_yet
            else:
                assert saddle == 'c', saddle
                had_centre_yet = True

def suite():
    suites = []
    suites.append(unittest.makeSuite(Pure))
    suites.append(unittest.makeSuite(Type))
    suites.append(unittest.makeSuite(SortGeneral))
    suites.append(unittest.makeSuite(SortByPosition))
    suites.append(unittest.makeSuite(SortByRealAndImag))
    suites.append(unittest.makeSuite(LocatePairs))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
