# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
from Powers import Powers
from Polynomial import Polynomial
from NormalFormIO import read_ascii_polynomial
from LieAlgebra import LieAlgebra
from EigenSystem import EigenSystem
from Diagonal import Diagonalizer, Complexifier
from LieTriangle import LieTriangle
from IntegralExtractor import IntegralExtractor

class Hill:
    
    def __init__(self, max_grade):
        dof = 3
        self.alg = LieAlgebra(dof)
        self.h = self.alg.isograde(self.hill_about_l1(), 0, max_grade+1)
        self.h_dc = self.diagonalize(self.alg, self.h)
        self.h_dr = self.h_dc.substitute(self.sub_r_into_c)
        self.k_dc = self.normalize_to_grade(max_grade)
        self.k_dr = self.k_dc.substitute(self.sub_r_into_c)

    def hill_about_l1(self):
        #print 'reading hill hamiltonian for test...'
        in_name = '../../test/ma-files/hill_l1_18_equi_to_tham.pol'
        in_file = open(in_name, 'r')
        h = read_ascii_polynomial(in_file, is_xxpp_format=1)
        in_file.close()
        #print 'done'
        terms = {Powers((2, 0, 0, 0, 0, 0)): -4.0, #x^2
                 Powers((0, 0, 2, 0, 0, 0)): +2.0, #y^2
                 Powers((0, 0, 0, 0, 2, 0)): +2.0, #z^2
                 Powers((0, 2, 0, 0, 0, 0)): +0.5, #px^2
                 Powers((0, 0, 0, 2, 0, 0)): +0.5, #py^2
                 Powers((0, 0, 0, 0, 0, 2)): +0.5, #pz^2
                 Powers((1, 0, 0, 1, 0, 0)): -1.0, #xpy
                 Powers((0, 1, 1, 0, 0, 0)): +1.0} #ypx
        assert len(terms) == 8
        h_2 = self.alg.polynomial(terms)
        assert (h_2 == self.alg.isograde(h, 2))
        return h

    def diagonalize(self, dof, h):
        tolerance = 5.0e-15

        self.diag = Diagonalizer(self.alg)
        self.eig = self.diag.compute_eigen_system(h, tolerance)
        self.diag.compute_diagonal_change()

        mat = self.diag.get_matrix_diag_to_equi()
        assert self.diag.matrix_is_symplectic(mat)
        
        sub_diag_into_equi = self.diag.matrix_as_vector_of_row_polynomials(mat)
        h_diag = h.substitute(sub_diag_into_equi)

        self.eq_type = self.eig.get_equilibrium_type()
        self.comp = Complexifier(self.alg, self.eq_type)
        self.sub_c_into_r = self.comp.calc_sub_complex_into_real()
        self.sub_r_into_c = self.comp.calc_sub_real_into_complex()
        h_comp = h_diag.substitute(self.sub_c_into_r)

        h_comp = h_comp.with_small_coeffs_removed(tolerance)
        assert (self.alg.is_diagonal_polynomial(self.alg.isograde(h_comp, 2)))

        return h_comp

    def normalize_to_grade(self, max_grade):
        steps = max_grade
        h_dc = self.alg.isograde(self.h_dc, 2, max_grade+1)
        h_dc_2 = self.alg.isograde(h_dc, 2)
        k_comp = self.alg.zero()
        w_comp = self.alg.zero()
        #the following is not very efficient for separating grades!
        h_list = [self.alg.isograde(h_dc, i + 2) for i in xrange(0, steps + 1)]
        nf = LieTriangle(self.alg, h_dc_2)
        for i in xrange(0, steps + 1):
            h_term = h_list[i]
            k_term, w_term = nf.compute_normal_form_and_generating_function(h_term)
            k_comp += k_term
            w_comp += w_term

        #print k_comp
        k_real = k_comp.substitute(self.sub_r_into_c)
        assert (k_real.imag().l_infinity_norm() < 1.0e-14)
        k_real = k_real.real()

        return k_comp

class NormalizeHill(unittest.TestCase):

    def setUp(self):
        self.hill = Hill(4) #6

    def tearDown(self):
        del self.hill

    def test_hill_basics(self):
        isograde = self.hill.alg.isograde
        #print 'diagonalized(R)[gra2]', isograde(self.hill.h_dr, 2)
        #print 'diagonalized(C)[gra2]', isograde(self.hill.h_dc, 2)
        self.assert_(len(isograde(self.hill.h_dr, 2)) == 5)
        self.assert_(len(isograde(self.hill.h_dc, 2)) == 3)
        #print 'normalized(R)[gra2]', isograde(self.hill.k_dr, 2)
        #print 'normalized(C)[gra2]', isograde(self.hill.k_dc, 2)
        self.assert_(isograde(self.hill.k_dr, 2) == isograde(self.hill.h_dr, 2))
        self.assert_(isograde(self.hill.k_dc, 2) == isograde(self.hill.h_dc, 2))

    def test_set_hamiltonian(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)

    def test_find_lost_integrals_and_non_simple_one_lost(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()

        lost = extractor._lost_integrals
        self.assert_(not lost) #none lost for hills

    def test_find_lost_integrals_and_non_simple_one_simple(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()

        simple_part = extractor._simple_part
        self.assert_(simple_part == self.hill.k_dc)

    def test_find_lost_integrals_and_non_simple_one_non_simple(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()

        non_simple = extractor._non_simple_integral
        self.assert_(not non_simple)

    def test_list_the_simple_integrals(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()

        kept = extractor._kept_integrals
        self.assert_(kept == [0, 1, 2])

    def test_express_simple_part_as_polynomial_over_all_integrals(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()

        simple_in_ints = extractor._simple_in_integrals
        #print 'simple part in all ints', simple_in_ints
        self.assert_(len(simple_in_ints) == len(self.hill.k_dc))

    def test_express_all_as_polynomial_over_all_integrals(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()
        extractor.express_both_parts_as_polynomial_over_all_integrals()

        total_in_ints = extractor._total_in_integrals
        #print 'whole polynomial in all ints', total_in_ints
        self.assert_(len(total_in_ints) == len(self.hill.k_dc))

    def test_express_all_integrals_in_normal_form_coords(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()
        extractor.express_all_integrals_in_normal_form_coords()

        ints_in_coords = extractor._integrals_in_coords
        self.assert_(len(ints_in_coords) == 3)
        #print 'all integrals in coords', ints_in_coords
        for i in ints_in_coords:
            self.assert_(len(i) == 1) #1-term

    def test_integral_expressions_against_original_hamiltonian(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()
        extractor.express_both_parts_as_polynomial_over_all_integrals()
        extractor.express_all_integrals_in_normal_form_coords()

        total_in_ints = extractor._total_in_integrals
        ints_in_coords = extractor._integrals_in_coords

        eval_at_ints = total_in_ints.substitute(ints_in_coords)
        diff = eval_at_ints - self.hill.k_dc
        self.assert_(diff.l_infinity_norm() < 1.0e-15, diff)

    def test_integral_expressions_against_original_hamiltonian_tmp(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()
        extractor.express_both_parts_as_polynomial_over_all_integrals()
        extractor.express_all_integrals_in_normal_form_coords()

        total_in_ints = extractor._total_in_integrals
        ints_in_coords = extractor._integrals_in_coords

        res = Polynomial(6)
        for po, co in total_in_ints.powers_and_coefficients():
            tmp = Polynomial.One(6)
            for i, po_i in enumerate(po):
                if po_i > 0:
                    tmp *= (ints_in_coords[i]**po_i)
            res += co * tmp
        eval_at_ints = res

        diff = eval_at_ints - self.hill.k_dc
        self.assert_(diff.l_infinity_norm() < 1.0e-10, diff)

    def test_integral_expressions_against_real_hamiltonian(self):
        extractor = IntegralExtractor(self.hill.alg)
        extractor.set_complex_normal_hamiltonian(self.hill.k_dc)
        extractor.find_lost_simple_integrals_and_non_simple_integral()
        extractor.list_the_simple_integrals()
        extractor.express_simple_part_as_polynomial_over_all_integrals()
        extractor.express_both_parts_as_polynomial_over_all_integrals()
        extractor.express_all_integrals_in_normal_form_coords()
        k_dc_in_ints = extractor._total_in_integrals
        ints_in_coords = extractor._integrals_in_coords
        r_ints_in_coords = [one_int.substitute(self.hill.sub_r_into_c) for one_int in ints_in_coords]
        res = Polynomial(6)
        for po, co in k_dc_in_ints.powers_and_coefficients():
            tmp = Polynomial.One(6)
            for i, po_i in enumerate(po):
                if po_i > 0:
                    tmp *= (r_ints_in_coords[i]**po_i)
            res += co * tmp
        eval_at_ints = res

        diff = eval_at_ints - self.hill.k_dr
        self.assert_(diff.l_infinity_norm() < 1.0e-15, diff)

        eval_at_ints = k_dc_in_ints.substitute(r_ints_in_coords)
        diff = eval_at_ints - self.hill.k_dr
        self.assert_(diff.l_infinity_norm() < 1.0e-15, diff)
        

def suite():
    suites = []
    suites.append(unittest.makeSuite(NormalizeHill))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
