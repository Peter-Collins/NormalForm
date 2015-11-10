This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

#setup logging
import os
import sys
import logging

#profiling
import profile
import pstats

#main imports
sys.path = [os.path.join(os.path.split(sys.path[0])[0], "py")] + sys.path
from Powers import Powers
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from NormalForm import NormalForm
from Diagonal import Complexifier
from SystemBath import SystemBath, new_random_system_bath
from PolynomialRingIO import PolynomialRingIO

stderr = sys.stderr

def create_empty_file(name):
    file = open(name, 'w')
    file.close()

def setup_logging(logfile_name, do_stream=False):
    global logger
    create_empty_file(logfile_name)
    file_handler = logging.FileHandler(logfile_name)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(stderr)
    stream_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    if do_stream:
        logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

def compute_diagonalisation(lie_algebra, h_er_trunc, max_degree, tolerance, system_prefix):
    """

    Here, just go as far as we can with double precision, and dump the
    output files.

    """
    
    nf = NormalForm(lie_algebra, h_er_trunc)
    nf.set_tolerance(tolerance)
    desired_grade = max_degree
    nf.desired_grade = max_degree
    logger.info('main computation started...')
    #setup:
    nf.h_er = nf.alg.isograde(nf.h_er, 0, desired_grade+1)
    nf.confirm_equilibrium_point()
    #output max grade
    print desired_grade
    #output h_er
    print nf.h_er
    #diagonalize
    nf.diagonalize()
    #make a polynomial ring io
    ring_io = PolynomialRingIO(lie_algebra)
    #output vector of polynomials for er in terms of dr
    er_from_dr = nf.er_in_terms_of_dr
    file_name = '%s--e-from-d.vpol' % system_prefix
    file_ostr = open(file_name, 'w')
    ring_io.write_sexp_vector_of_polynomials(file_ostr, er_from_dr)
    file_ostr.close()
    #output equilibrium type
    print nf.eig.get_equilibrium_type()
    #compute the complexification maps
    eq_type = nf.eig.get_equilibrium_type()
    com = Complexifier(lie_algebra, eq_type)
    r_from_c = com.calc_sub_complex_into_real()
    c_from_r = com.calc_sub_real_into_complex()
    #output the complexification maps
    file_name = '%s--r-from-c.vpol' % system_prefix
    file_ostr = open(file_name, 'w')
    ring_io.write_sexp_vector_of_polynomials(file_ostr, r_from_c)
    file_ostr.close()

def dump_isogrades(graded_algebra, poly, max_grade, system_prefix, poly_name):
    logger.info('dumping taylor series...')
    file_prefix = '%s--%s' % (system_prefix, poly_name)
    ring_io = PolynomialRingIO(graded_algebra)
    for grade in xrange(max_grade + 1):
        file_name = '%s--grade-%d.pol' % (file_prefix, grade)
        file_ostr = open(file_name, 'w')
        p_grade = graded_algebra.isograde(poly, grade)
        ring_io.write_sexp_polynomial(file_ostr, p_grade)
        file_ostr.close()
    logger.info('...done')

def compute_normal_form(lie_algebra, h_er_trunc, max_degree, tolerance):
    #perform computations
    nf = NormalForm(lie_algebra, h_er_trunc)
    nf.set_tolerance(tolerance)
    nf.perform_all_computations(max_degree)

def main():
    """

    13/02/2005: pending a theoretical understanding of gauge
    transformations, etc., for the Lennard-Jones clusters, in
    particular Argon-6 (for which we have minimum energy path data
    from Wales), I will look at a system-bath model of the same size.
    Argon-6 (original coordinate, not centre-of-mass coordinates) has
    6 particles and thus the original (pre-gauge) Hamiltonian has 6 x
    6 = 36 variables, grouped into 18 degrees of freedom.  An 18 dof
    system bath model has 17 bath modes.

    """

    #input system bath parameters
    n_bath_modes = 8 #17
    system_mass = 1.0
    imaginary_harmonic_frequency_at_barrier = -0.5
    reciprocal_barrier_height_above_well_bottom = 0.6
    random_seed = 54321

    #input normal form computation parameters
    max_degree = 20
    tolerance = 5.0e-12

    #start logging
    system_prefix = 'system-bath--dof-%d' % (1 + n_bath_modes)
    file_prefix = '%s--grade-%d' % (system_prefix, max_degree)
    logfile_name = file_prefix + '.log'
    setup_logging(logfile_name, do_stream=True)

    #build system bath and retrieve algebra
    sb = new_random_system_bath(n_bath_modes,
                                system_mass,
                                imaginary_harmonic_frequency_at_barrier,
                                reciprocal_barrier_height_above_well_bottom,
                                random_seed)
    lie = sb.lie_algebra()
    h_er = sb.hamiltonian_real()

    #dump the Taylor series
    dump_isogrades(lie, h_er, max_degree, system_prefix, 'h-from-e')

    #truncate at requied grade
    h_er_trunc = lie.isograde(h_er, 0, max_degree + 1)

    #perform diagonalisation, normal form, coordinate changes, etc.
    compute_diagonalisation_only = True #False
    if compute_diagonalisation_only:
        diagonalisation_degree = 2
        compute_diagonalisation(lie, h_er_trunc, diagonalisation_degree + 1, tolerance, system_prefix)
    else:
        compute_normal_form(lie, h_er_trunc, max_degree, tolerance)

if __name__ == '__main__':
    #setup profiling
    profile_name = 'system-bath.prof'
    proportion_worst_results = 0.125
    create_empty_file(profile_name)
    
    #perform computations
    profile.run('main()', profile_name)

    #output profile
    stats = pstats.Stats(profile_name)
    stats.strip_dirs()
    stats.sort_stats('time', 'cum', 'calls')
    print '[1] Statistics:'
    stats.print_stats(proportion_worst_results)
    print '[2] Callers for the above:'
    stats.print_callers(proportion_worst_results)


