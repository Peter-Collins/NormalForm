"""

AUTHOR: Peter Collins, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: NfExample

PURPOSE:

Provide various example configuration and Hamiltonians for the
normalization algorithms.

NOTES:

These are copied from  projects/nf--unified--0.5/src/py/NormalFormTest.py.

"""

from Powers import Powers
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from NormalFormIO import read_ascii_polynomial
import logging
logger = logging.getLogger()
ma_path = '../../test/ma-files/'

def read_ma_poly(name):
    global ma_path
    import os
    import sys
    logger.info('reading')
    logger.info(name)
    ma_path = os.path.join(sys.path[0], ma_path)
    in_file = open(os.path.join(ma_path, name), 'r')
    p = read_ascii_polynomial(in_file,
                              is_xxpp_format=True)
    in_file.close()
    logger.info('done')
    return p

class NfExample:
    """

    An example system on which to perform the normal form computation.
    Sublass so as to provide self.dof and self.h_er or otherwise
    override the accessors get_dof() and get_h_er().

    """
    def get_lie_algebra(self):
        return self.lie
    def get_dof(self):
        return self.get_lie_algebra().dof()
    def get_h_er(self):
        return self.h_er

class Hill(NfExample):
    """

    Hill's equations for 3-DoF.

    """
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'hill_l1_18'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')
        terms = {Powers((2, 0, 0, 0, 0, 0)): -4.0, #x^2
                 Powers((0, 0, 2, 0, 0, 0)): +2.0, #y^2
                 Powers((0, 0, 0, 0, 2, 0)): +2.0, #z^2
                 Powers((0, 2, 0, 0, 0, 0)): +0.5, #px^2
                 Powers((0, 0, 0, 2, 0, 0)): +0.5, #py^2
                 Powers((0, 0, 0, 0, 0, 2)): +0.5, #pz^2
                 Powers((1, 0, 0, 1, 0, 0)): -1.0, #xpy
                 Powers((0, 1, 1, 0, 0, 0)): +1.0} #ypx
        assert len(terms) == 8
        h_2 = self.lie.polynomial(terms)
        assert (h_2 == self.lie.isograde(self.h_er, 2))
        del h_2

class RydbergCircular(NfExample):
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'cp_saddle_18'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class Hcn(NfExample):
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'hcn_saddle_10'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class RydbergCrossed(NfExample):
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'rydberg_saddle_16'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class SingleSaddle(NfExample):
    def __init__(self):
        dof = 1
        self.lie = LieAlgebra(dof)
        self.prefix = 'saddle_1dof'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class SingleCentre(NfExample):
    def __init__(self):
        dof = 1
        self.lie = LieAlgebra(dof)
        self.prefix = 'centre_1dof'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class SingleSaddleBad(NfExample):
    def __init__(self):
        dof = 1
        self.lie = LieAlgebra(dof)
        self.prefix = 'bad_1dof'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class AnotherSaddle(NfExample):
    def __init__(self):
        dof = 1
        self.lie = LieAlgebra(dof)
        self.h_er = self.lie.polynomial({Powers((1, 1)): +1.0,
                                         Powers((3, 0)): -1.0,
                                         Powers((4, 0)): -1.0,
                                         Powers((5, 0)): -1.0})

class Rftbp(NfExample):
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'rftbp_eros_eq1_10'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class Crtbp(NfExample):
    def __init__(self):
        dof = 3
        self.lie = LieAlgebra(dof)
        self.prefix = 'crtbp_sun_neptune_l1_14'
        self.h_er = read_ma_poly(self.prefix+'_equi_to_tham.pol')

class SystemBath(NfExample):
    def __init__(self, n_bath_modes):
        from SystemBath import new_random_system_bath

        #input system bath parameters
        #n_bath_modes = 8 #17
        system_mass = 1.0
        imaginary_harmonic_frequency_at_barrier = -0.5
        reciprocal_barrier_height_above_well_bottom = 0.6
        random_seed = 54321

        #prefix for logging etc
        self.prefix = 'system-bath--dof-%d' % (1 + n_bath_modes)
        #file_prefix = '%s--grade-%d' % (self.prefix, max_degree)

        #build system bath and retrieve algebra
        sb = new_random_system_bath(n_bath_modes,
                                    system_mass,
                                    imaginary_harmonic_frequency_at_barrier,
                                    reciprocal_barrier_height_above_well_bottom,
                                    random_seed)
        self.lie = sb.lie_algebra()
        self.h_er = sb.hamiltonian_real()
        assert not (self.h_er.n_vars()%2)

