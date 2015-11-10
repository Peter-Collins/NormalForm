#!/usr/local/bin/python2.3

# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
#import X

class TestXFeature(unittest.TestCase):

    def test_that_test_fixture_is_called(self):
        self.failUnless(False)

def suite():
    suites = []
    suites.append(unittest.makeSuite(TestXFeature))
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
