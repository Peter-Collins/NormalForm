"""

AUTHOR: Dr Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: AccuracyTest

PURPOSE:

Unit tests for the Accuracy module.

NOTES:

"""

import unittest

from Accuracy import trunc

class ForcePure(unittest.TestCase):
    """

    Test routines for converting complex scalars and tuples with tiny
    real/imag parts into pure real/imag.

    """

    def setUp(self):
        def closure(x):
            return trunc(x, tolerance=1.0e-15)
        self.trunc = trunc

    def test_float(self):
        c = 0.5
        self.assertEquals(self.trunc(c), 0.5)
        c = -2.0e-15
        self.assertEquals(self.trunc(c), -2.0e-15)
        c = -3.0e-17
        self.assertEquals(self.trunc(c), 0.0)

    def test_real(self):
        c = complex(0.5, 0.0)
        self.assertEquals(self.trunc(c), 0.5)
        c = complex(0.5, 1.0e-16)
        self.assertEquals(self.trunc(c), 0.5)

    def test_imag(self):
        c = complex(0.0, -2.3)
        self.assertEquals(self.trunc(c), -2.3J)
        c = complex(-5.0e-16, -2.3)
        self.assertEquals(self.trunc(c), -2.3J)

    def test_complex(self):
        c = complex(0.5, 0.1)
        self.assertEquals(self.trunc(c), 0.5+0.1J)
        c = complex(1.0e-16, -2.3e-18)
        self.assertEquals(self.trunc(c), 0.0)

def suite():
    suites = []
    suites.append(unittest.makeSuite(ForcePure))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
