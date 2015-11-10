# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from EigenValueVectorPair import EigenValueVectorPair

class EigenValueVectorPairTest(unittest.TestCase):

    def test_basic_access(self):
        """
        
        Simply construct and access the components.

        """
        val = 3.0
        vec = (-1.0, 0.5, -2.0)
        val_vec = EigenValueVectorPair(val, vec)
        self.assertEquals(val_vec.val, val)
        self.assertEquals(val_vec.vec, vec)

def suite():
    suites = []
    suites.append(unittest.makeSuite(EigenValueVectorPairTest))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
