# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from Powers import Powers
from CoordinateChange import *

class EmptyGeneratingFunction(unittest.TestCase):

    """

    With a zero, or constant, generating function, all functions
    should be preserved under the Lie transform.

    At present, though, our formulation will only preserve terms in a
    polynomial that are of degree at least one.  In other words, the
    constant term is neglected.

    Mathematically, the constant term should survive during a Lie
    transform, so this is probably an artifact of the Lie triangle
    method.

    We could always simply copy any degree zero terms across at the
    beginning of the transform.  We could likewise do the same during
    normalisation?

    """

    def test_diag_in_norm(self):
        dof = 2
        steps = 5
        alg = LieAlgebra(dof)
        w_list = [2.0*alg.one()]*(steps+1)
        for index in xrange(alg.n_vars()):
            f_s = alg.coordinate_monomial(index)+2.3*alg.one()
            x_i_list = []
            x_ij_dict = {}
            changer = CoordinateChange(alg, w_list)
            changer.express_diag_in_norm(f_s, x_i_list, x_ij_dict, steps)
            id_s = alg.coordinate_monomial(index)
            self.assert_(len(x_i_list) >= 1)
            self.assert_(x_i_list[0] == id_s, x_i_list[0])
            for i in xrange(1, len(x_i_list)):
                self.assert_(x_i_list[i] == alg.zero())

    def test_norm_in_diag(self):
        dof = 2
        steps = 5
        alg = LieAlgebra(dof)
        w_list = [3.0*alg.one()]*(steps+1)
        for index in xrange(alg.n_vars()):
            f_s = alg.coordinate_monomial(index)
            x_i_list = []
            x_ij_dict = {}
            changer = CoordinateChange(alg, w_list)
            changer.express_norm_in_diag(f_s, x_i_list, x_ij_dict, steps)
            id_s = alg.coordinate_monomial(index)
            self.assert_(len(x_i_list) >= 1)
            self.assert_(x_i_list[0] == id_s, x_i_list[0])
            for i in xrange(1, len(x_i_list)):
                self.assert_(x_i_list[i] == alg.zero())

def suite():
    suites = []
    suites.append(unittest.makeSuite(EmptyGeneratingFunction))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
