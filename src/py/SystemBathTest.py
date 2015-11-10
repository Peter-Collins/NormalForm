# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from math import sqrt

from Polynomial import Polynomial
from EigenSystem import EigenSystem
from Diagonal import Diagonalizer, Complexifier
from SystemBath import *

imag_harm_freq_barr = -0.5

class EquilibriumPoint(unittest.TestCase):
    """

    Ensure that the generated Hamiltonian has an equilibrium
    point at the origin.

    """

    def test_no_constant_term(self):
        """

        Ensure that there is no constant term.

        """
        sb = new_random_system_bath(0, 1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        h_er = sb.hamiltonian_real()
        self.assert_(not h_er.homogeneous(0))

    def test_no_linear_term(self):
        """

        Ensure that there is no linear term.

        """
        sb = new_random_system_bath(0, 1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        h_er = sb.hamiltonian_real()
        self.assert_(not h_er.homogeneous(1))

class Saddle(unittest.TestCase):
    """

    Test the system bath with no bath modes.  This should be a simple saddle.

    """

    def test_zero_and_first_degree(self):
        sb = new_random_system_bath(0, 1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        h = sb.hamiltonian_real()
        self.assert_(h.homogeneous(0) == Polynomial(2))
        self.assert_(h.homogeneous(1) == Polynomial(2))

    def test_quadratic_terms(self):
        tolerance = 1.0e-15
        n = 0 #bath modes
        
        sb = new_random_system_bath(n, 1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        lie = sb.lie_algebra()
        h = sb.hamiltonian_real()

        #Sanity checks:
        assert h.degree() == 4
        assert h.n_vars() == 2*n+2
        assert len(h) == (3) + (2*n) + (n) #system+bath+coupling

        h_2 = h.homogeneous(2)
        self.assert_(h_2) #non-zero
        diag = Diagonalizer(lie)
        eig = diag.compute_eigen_system(h_2, tolerance)
        diag.compute_diagonal_change()

        eq_type = eig.get_equilibrium_type()
        self.assertEquals(eq_type, 's')

        eigs = [pair.val for pair in eig.get_raw_eigen_value_vector_pairs()]

        mat = diag.get_matrix_diag_to_equi()
        assert diag.matrix_is_symplectic(mat)

        sub_diag_into_equi = diag.matrix_as_vector_of_row_polynomials(mat)
        comp = Complexifier(lie, eq_type)
        sub_complex_into_real = comp.calc_sub_complex_into_real()

        h_2_diag = h_2.substitute(sub_diag_into_equi)
        h_2_comp = h_2_diag.substitute(sub_complex_into_real)
        h_2_comp = h_2_comp.with_small_coeffs_removed(tolerance)
        self.assert_(lie.is_diagonal_polynomial(h_2_comp))

        h_diag = h.substitute(sub_diag_into_equi)
        h_comp = h_diag.substitute(sub_complex_into_real)
        h_comp = h_comp.with_small_coeffs_removed(tolerance)

        h_comp_2 = h_comp.homogeneous(2)
        self.assert_((h_comp_2-h_2_comp).l1_norm() < tolerance,
                     '%s != %s'%(h_comp_2, h_2_comp))
        self.assert_(not h_comp.homogeneous(1))

class SaddleCentre(unittest.TestCase):
    """

    Test the system bath with a single bath modes.  This should have
    saddle-centre structure.

    """

    def test_zero_and_first_degree(self):
        n_bath_modes = 1

        lie = LieAlgebra(n_bath_modes+1)
        sb = new_random_system_bath(n_bath_modes,
                                    1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        h = sb.hamiltonian_real()
        self.assert_(h.homogeneous(0) == Polynomial(n_bath_modes*2+2))
        self.assert_(h.homogeneous(1) == Polynomial(n_bath_modes*2+2))

    def test_quadratic_terms(self):
        n_bath_modes = 1
        tolerance = 1.0e-14
        
        sb = new_random_system_bath(n_bath_modes,
                                    1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        lie = sb.lie_algebra()
        h = sb.hamiltonian_real()

        #Sanity checks:
        assert h.degree() == 4
        assert h.n_vars() == 2*n_bath_modes+2
        assert len(h) == (3) + (2*n_bath_modes) + (n_bath_modes)

        h_2 = h.homogeneous(2)
        self.assert_(h_2) #non-zero
        diag = Diagonalizer(lie)
        eig = diag.compute_eigen_system(h_2, tolerance)
        diag.compute_diagonal_change()

        eq_type = eig.get_equilibrium_type()
        self.assertEquals(eq_type, 'sc')

        eigs = [pair.val for pair in eig.get_raw_eigen_value_vector_pairs()]

        mat = diag.get_matrix_diag_to_equi()
        assert diag.matrix_is_symplectic(mat)

        sub_diag_into_equi = diag.matrix_as_vector_of_row_polynomials(mat)
        comp = Complexifier(lie, eq_type)
        sub_complex_into_real = comp.calc_sub_complex_into_real()

        h_2_diag = h_2.substitute(sub_diag_into_equi)
        h_2_comp = h_2_diag.substitute(sub_complex_into_real)
        h_2_comp = h_2_comp.with_small_coeffs_removed(tolerance)
        self.assert_(lie.is_diagonal_polynomial(h_2_comp))

        h_diag = h.substitute(sub_diag_into_equi)
        h_comp = h_diag.substitute(sub_complex_into_real)
        h_comp = h_comp.with_small_coeffs_removed(tolerance)

        h_comp_2 = h_comp.homogeneous(2)
        self.assert_((h_comp_2-h_2_comp).l1_norm() < tolerance,
                     '%s != %s'%(h_comp_2, h_2_comp))
        self.assert_(not h_comp.homogeneous(1))

class SaddleCentreCentre(unittest.TestCase):
    """

    Test the system bath with two bath modes; saddle-centre-centre structure.

    """

    def test_zero_and_first_degree(self):
        n_bath_modes = 2

        lie = LieAlgebra(n_bath_modes+1)
        sb = new_random_system_bath(n_bath_modes,
                                    1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        h = sb.hamiltonian_real()
        self.assert_(h.homogeneous(0) == Polynomial(n_bath_modes*2+2))
        self.assert_(h.homogeneous(1) == Polynomial(n_bath_modes*2+2))

    def test_quadratic_terms(self):
        n_bath_modes = 2
        tolerance = 1.0e-14
        
        sb = new_random_system_bath(n_bath_modes,
                                    1.0, imag_harm_freq_barr, 0.6, random_seed=54321)
        lie = sb.lie_algebra()
        h = sb.hamiltonian_real()

        #Sanity checks:
        assert h.degree() == 4
        assert h.n_vars() == 2*n_bath_modes+2
        assert len(h) == (3) + (2*n_bath_modes) + (n_bath_modes)

        h_2 = h.homogeneous(2)
        self.assert_(h_2) #non-zero
        diag = Diagonalizer(lie)
        eig = diag.compute_eigen_system(h_2, tolerance)
        diag.compute_diagonal_change()

        eq_type = eig.get_equilibrium_type()
        self.assertEquals(eq_type, 'scc')

        eigs = [pair.val for pair in eig.get_raw_eigen_value_vector_pairs()]

        mat = diag.get_matrix_diag_to_equi()
        assert diag.matrix_is_symplectic(mat)

        sub_diag_into_equi = diag.matrix_as_vector_of_row_polynomials(mat)
        comp = Complexifier(lie, eq_type)
        sub_complex_into_real = comp.calc_sub_complex_into_real()

        h_2_diag = h_2.substitute(sub_diag_into_equi)
        h_2_comp = h_2_diag.substitute(sub_complex_into_real)
        h_2_comp = h_2_comp.with_small_coeffs_removed(tolerance)
        self.assert_(lie.is_diagonal_polynomial(h_2_comp))

        h_diag = h.substitute(sub_diag_into_equi)
        h_comp = h_diag.substitute(sub_complex_into_real)
        h_comp = h_comp.with_small_coeffs_removed(tolerance)

        h_comp_2 = h_comp.homogeneous(2)
        self.assert_((h_comp_2-h_2_comp).l1_norm() < tolerance,
                     '%s != %s'%(h_comp_2, h_2_comp))
        self.assert_(not h_comp.homogeneous(1))

def suite():
    suites = []
    suites.append(unittest.makeSuite(EquilibriumPoint))
    suites.append(unittest.makeSuite(SaddleCentreCentre))
    suites.append(unittest.makeSuite(SaddleCentre))
    suites.append(unittest.makeSuite(Saddle))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')

