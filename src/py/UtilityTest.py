# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

import Utility

class Factorial(unittest.TestCase):

    def test_values(self):
        self.assertEquals(Utility.factorial(0), 1)
        self.assertEquals(Utility.factorial(1), 1)
        self.assertEquals(Utility.factorial(2), 2)
        self.assertEquals(Utility.factorial(5), 120)

class Binomial(unittest.TestCase):

    def setUp(self):
        self._choose = Utility.binomial

    def _choose_old(self, n, k):
        if n < 0: raise ValueError
        if k < 0 or k > n: raise ValueError
        b = [0] * (n + 1)
        b[0] = 1
        for i in xrange(1, n + 1):
            b[i] = 1
            j = i - 1
            while j > 0:
                b[j] += b[j - 1]
                j -= 1
        return b[k]

    def test_neg(self):
        try:
            self._choose(-1, 0)
        except ValueError:
            pass
        else:
            assert 0, 'negative first arg to Binomial should raise'

    def test_zero_choose_zero(self):
        self.assertEquals(self._choose(0, 0), 1) #?check?

    def test_zero_choose_negative(self):
        try:
            self._choose(0, -1)
        except ValueError:
            pass
        else:
            assert 0, 'choosing negative in Binomial should raise'

    def test_zero_choose_positive(self):
        try:
            self._choose(0, 1)
        except ValueError:
            pass
        else:
            assert 0, 'choosing positive in Binomial zero should raise'

    def test_values(self):
        examples = {(0, 0): 1,
                    (1, 0): 1,
                    (1, 1): 1,
                    (2, 0): 1,
                    (2, 1): 2,
                    (2, 2): 1,
                    (4, 2): 6,
                    (9, 8): 9}
        for (n, k), b in examples.iteritems():
            self.assertEquals(self._choose(n, k), b)

    def test_identity0(self):
        for n in xrange(10):
            lhs = self._choose(n, 0)
            rhs = self._choose(n, n)
            self.assert_(lhs == rhs == 1)

    def test_identity1(self):
        for n in xrange(10):
            for k in xrange(n+1):
                lhs = self._choose(n, k)
                rhs = self._choose(n, n-k)
                self.assert_(lhs == rhs)

    def test_identity2_does_not_work_due_to_negatives(self):
        for n in xrange(10):
            for k in xrange(n+1):
                lhs = self._choose(n, k)*(-1 ** k)
                try:
                    rhs = self._choose(k-n-1, k)
                except ValueError:
                    pass
                else:
                    self.assert_(lhs == rhs)

    def test_identity3(self):
        for n in xrange(10):
            for k in xrange(n+1-1):
                lhs = self._choose(n, k+1)
                rhs = self._choose(n, k)*(n-k)/(k+1)
                self.assert_(lhs == rhs)

    def test_identity4(self):
        for n in xrange(10):
            for k in xrange(1, n+1):
                lhs = self._choose(n+1, k)
                rhs = self._choose(n, k)+self._choose(n, k-1)
                self.assert_(lhs == rhs)

    def test_cf_implementations(self):
        for i in xrange(1, 10):
            for j in xrange(1, i+1):
                b0 = self._choose_old(i, j)
                b1 = self._choose(i, j)
                self.assertEquals(b0, b1, (i,j,b0,b1))

    def test_range_downward_misses_final_value(self):
        c = 5
        for i in range(5, 0, -1):
            c -= 1
        self.assertEquals(c, 0)

    def test_xrange_downward_misses_final_value(self):
        c = 5
        for i in xrange(5, 0, -1):
            c -= 1
        self.assertEquals(c, 0)

def suite():
    suites = []
    suites.append(unittest.makeSuite(Factorial))
    suites.append(unittest.makeSuite(Binomial))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
