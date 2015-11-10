"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: ComplexTest

PURPOSE:

Unit-tests for utilities for handling complex numbers.

NOTES:

"""

import unittest

import Complex

def _is_bracketed(st):
    """Test that a string is parenthetic."""
    return (st.startswith('(') and st.endswith(')'))

class PrettyPrinting(unittest.TestCase):
    """

    Test the pretty-printing facilities for complex coefficients.

    """
    def setUp(self):
        """Construct some strings representing one and zero."""
        self._prec = 12
        self._zero = ('+0.'+('0'*self._prec)+'e+00')
        self._one = ('+1.'+('0'*self._prec)+'e+00')
    def test_zero(self):
        """Pretty-print zero complex number."""
        st = Complex.pretty_complex(0)
        self.assert_(_is_bracketed(st), st)
        self.assert_(st[1:].startswith(self._zero), st)
    def test_one(self):
        """Pretty-print pure real one."""
        st = Complex.pretty_complex(1)
        self.assert_(_is_bracketed(st), st)
        self.assert_(st[1:].startswith(self._one), st)
    def test_j(self):
        """Pretty-print pure imaginary unit J."""
        st = Complex.pretty_complex(1J)
        self.assert_(_is_bracketed(st), st)
        self.assert_(st[:-1].endswith(self._one+'j'), st)
    def test_one_plus_j(self):
        """Pretty-print mixed complex number one plus J."""
        st = Complex.pretty_complex(1+1J)
        self.assert_(_is_bracketed(st), st)
        self.assert_(st[1:].startswith(self._one+' '), st)
        self.assert_(st[:-1].endswith(self._one+'j'), st)
    def test_example(self):
        """Pretty-print an example in scientific format."""
        st = Complex.pretty_complex(-0.25+1.125J)
        re = '-2.5'+('0'*(self._prec-1))+'e-01'
        im = '+1.125'+('0'*(self._prec-3))+'e+00'
        self.assertEquals(st, '(%s %sj)'%(re, im))

def suite():
    """Automatically collect tests into a test suite."""
    suites = []
    suites.append(unittest.makeSuite(PrettyPrinting))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
