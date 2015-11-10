#!/usr/local/bin/python2.3

# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
import Uuid

class Examples(unittest.TestCase):

    def test_make_some_uuids(self):
        ids = []
        for i in xrange(1000):
            id = Uuid.uuidgen()
            ids.append(id)
        for x, y in zip(ids[:-1], ids[1:]):
            assert x != y

def suite():
    suites = []
    suites.append(unittest.makeSuite(Examples))
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
