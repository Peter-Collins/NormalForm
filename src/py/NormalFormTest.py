# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import os
from sys import stderr

import logging

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

#clear any existing log file:
logfile_name = 'NormalFormTest.log'
dummy = open(logfile_name, 'w')
dummy.close()

file_handler = logging.FileHandler(logfile_name)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler(stderr)
stream_handler.setFormatter(formatter)

logger = logging.getLogger() #('NormalFormTest')
logger.setLevel(logging.INFO)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)

using_hotshot = False

if using_hotshot:
    import hotshot
    import hotshot.stats
else:
    import profile
    import pstats

from Powers import Powers
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial, read_ascii_polynomial
use_generator = 1
if use_generator:
    from NormalForm import NormalForm
else:
    from NormalFormOld import NormalForm
from SystemBath import SystemBath, new_random_system_bath

ma_path = '../../test/ma-files/'

def read_ma_poly(name):
    logger.info('reading')
    logger.info(name)
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
        sb = new_random_system_bath(n_bath_modes,
                                    1.0,
                                    -0.5,
                                    0.6,
                                    random_seed=54321)
        self.lie = sb.lie_algebra()
        self.h_er = sb.hamiltonian_real()
        assert not (self.h_er.n_vars()%2)

def main():
    if 0:
        example = Hcn()
        example = RydbergCircular()
        example = RydbergCrossed()
        example = SingleSaddleBad()
        example = SingleCentre()
        example = SingleSaddle()
        example = AnotherSaddle()
        example = Rftbp()
        example = Crtbp()
        example = Hill()
    if 1:
        example = Hill()
        #example = SystemBath(9) #SystemBath(49) #49=50dof, 100vars
        
    #for hill, degree greater than 5 we seem to get huge errors in coord changes.
    lie = example.get_lie_algebra()
    degree = 6
    h_er_trunc = lie.isograde(example.get_h_er(), 0, degree+1)

    #tolerance is not really unified in the code.
    tolerance = 5.0e-14
    nf = NormalForm(lie, h_er_trunc)
    nf.set_tolerance(tolerance)
    nf.perform_all_computations(degree)

def create_empty_file(name):
    file = open(name, 'w')
    file.close()

if __name__ == '__main__':
    profile_name = 'NormalFormTest.prof'
    proportion_worst_results = 0.125

    create_empty_file(profile_name)
    
    if using_hotshot:
        profiler = hotshot.Profile(profile_name)
        benchtime = profiler.runcall(main) #main call
        profiler.close()
        stats = hotshot.stats.load(profile_name)
    else:
        profile.run('main()', profile_name) #main call
        stats = pstats.Stats(profile_name)

    stats.strip_dirs()
    stats.sort_stats('time', 'cum', 'calls')
    print '[1] Statistics:'
    stats.print_stats(proportion_worst_results)
    print '[2] Callers for the above:'
    stats.print_callers(proportion_worst_results)

