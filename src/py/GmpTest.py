# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

"""

Test the basic Python wrappers to the GNU Multi-Precision library.

"""

import unittest
import pickle as pickle

try:
    from gmpy import mpf, mpz
except:
    print 'NO GMP INSTALLED!'

def equal_to_bits_prec(a, b, bits):
    return abs(a-b) < 1.0/float(1L<<bits)

class BasicFloats(unittest.TestCase):

    def test_make_one(self):
        bits = 128
        base = 10
        ones = []
        ones.append(mpf('1.0', bits, base))
        ones.append(mpf(1.0, bits))
        ones.append(mpf(1.0))
        for one in ones:
            self.assert_(one == 1.0)

    def test_make_one_plus_epsilon(self):
        bits = 128
        base = 10
        one = mpf('1', bits, base)
        eps = mpf('1e-38', bits, base)
        self.assert_(one != eps)
        self.assert_(eps > 0.0)
        self.assert_(one+eps > 1.0)
        self.assert_(float(one+eps) == 1.0)

    def test_addition(self):
        bits = 512
        base = 10
        one = mpf('1', bits, base)
        eps = mpf('1e-38', bits, base)
        self.assert_(mpz(1.0/((one+eps)-one)) == 1e+38)

    def test_numbers_are_immutable_in_present_version(self):
        one = mpf('1.0', 32)
        one_orig = one
        one += 1.0
        self.assert_(one == 2.0)
        self.assert_(one_orig == 1.0)

    def test_cannot_pickle_and_unpickle_in_raw_form(self):
        bits = 512
        base = 10
        one = mpf('1', bits, base)
        try:
            s = pickle.dumps(one)
        except pickle.PicklingError:
            pass
        else:
            assert 0, 'expected not to be able to pickle in gmpy! new version?'

    def test_can_pickle_string_and_retrieve(self):
        bits = 512
        base = 10
        num = mpf('3.23987239874534573452439582735234e-10', bits, base)
        one = mpf('1', bits, base)
        num = num+one
        def mpf_to_pkl(mpf):
            return (mpf.digits(), mpf.getrprec())
        def pkl_to_mpf((digit_str, prec)):
            return mpf(digit_str, prec, 10)
        s = pickle.dumps(mpf_to_pkl(num))
        num_new = pkl_to_mpf(pickle.loads(s))
        self.assert_(equal_to_bits_prec(num, num_new, bits))

    def test_equal_to_bits_prec(self):
        bits = 4
        num = mpf('1.0', bits)+1.0/mpf(1<<(bits+1), bits)
        for i in xrange(bits*2):
            if i <= bits:
                self.assert_(equal_to_bits_prec(1.0, num, i))
            else:
                self.assert_(not equal_to_bits_prec(1.0, num, i))

def suite():
    suites = []
    suites.append(unittest.makeSuite(BasicFloats))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
