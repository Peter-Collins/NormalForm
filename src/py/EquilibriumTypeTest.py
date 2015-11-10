# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from EquilibriumType import *

class Type(unittest.TestCase):

    def test_empty_type_to_str(self):
        eq_type = ''
        out = equilibrium_type_to_str(eq_type)
        self.assertEquals(out, '')

    def test_saddle_type_to_str(self):
        eq_type = 's'
        out = equilibrium_type_to_str(eq_type)
        self.assertEquals(out, 'saddle')

    def test_centre_type_to_str(self):
        eq_type = 'c'
        out = equilibrium_type_to_str(eq_type)
        self.assertEquals(out, 'centre')

    def test_zero_type_to_str(self):
        eq_type = '0'
        out = equilibrium_type_to_str(eq_type)
        self.assertEquals(out, 'zero')

    def test_empty_type_to_tex(self):
        eq_type = ''
        out = equilibrium_type_to_tex(eq_type)
        self.assertEquals(out, '')

    def test_sadle_type_to_tex(self):
        eq_type = 's'
        out = equilibrium_type_to_tex(eq_type)
        self.assertEquals(out, 'saddle')

    def test_centre_type_to_tex(self):
        eq_type = 'c'
        out = equilibrium_type_to_tex(eq_type)
        self.assertEquals(out, 'centre')

    def test_zero_type_to_tex(self):
        eq_type = '0'
        out = equilibrium_type_to_tex(eq_type)
        self.assertEquals(out, 'zero')

    def test_equilibrium_type_to_str(self):
        eq_type = 'scssscc000'
        out = equilibrium_type_to_str(eq_type)
        exp = 'saddle x centre x saddle^(3) x centre^(2) x zero^(3)'
        self.assertEquals(out, exp, '%s\n%s'%(out, exp))

    def test_equilibrium_type_to_tex(self):
        eq_type = 'scssscc000'
        out = equilibrium_type_to_tex(eq_type)
        exp = 'saddle$\\times$centre$\\times$saddle$^{3}$$\\times$centre$^{2}$$\\times$zero$^{3}$'
        self.assertEquals(out, exp, '%s\n%s'%(out, exp))

def suite():
    suites = []
    suites.append(unittest.makeSuite(Type))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
