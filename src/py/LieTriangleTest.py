"""

Author: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

"""

import unittest
import os.path

import LieTriangle

from PolynomialRing import PolynomialRing
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from NormalFormIO import read_ascii_polynomial
from IsogradeInnerTaylorCoeffs import IsogradeInnerTaylorCoeffs
from Powers import Powers
from EigenSystem import EigenSystem
from Diagonal import Diagonalizer, Complexifier

test_dir = '../../test/'

class SumPolynomialList(unittest.TestCase):

    """Ensure that factorials appear correctly."""

    def test_example(self):
        ring = PolynomialRing(3)
        p0 = Polynomial(3)
        p0 += Polynomial.Monomial((0,0,0))
        assert ring.grade(p0) == 0
        p1 = Polynomial(3)
        p1 += Polynomial.Monomial((0,1,0))
        p1 += Polynomial.Monomial((0,0,1))
        assert ring.grade(p1) == 1
        p2 = Polynomial(3)
        p2 += Polynomial.Monomial((0,2,0))
        p2 += Polynomial.Monomial((0,1,1))
        p2 += Polynomial.Monomial((1,0,1))
        assert ring.grade(p2) == 2
        p3 = Polynomial(3)
        p3 += Polynomial.Monomial((3,0,0))
        p3 += Polynomial.Monomial((0,1,2))
        assert ring.grade(p3) == 3
        iso = IsogradeInnerTaylorCoeffs(ring, 0)
        p_sum = iso.list_to_poly([1*p0, 1*p1, 2*p2, 6*p3])
        self.assert_(ring.isograde(p_sum, 0) == p0)
        self.assert_(ring.isograde(p_sum, 1) == p1)
        self.assert_(ring.isograde(p_sum, 2) == p2)
        self.assert_(ring.isograde(p_sum, 3) == p3)

class DbKeys(unittest.TestCase):

    """Test that database keys for the $H_{ij}$ are constructed
    correctly."""

    def test_one(self):
        self.assertEquals(LieTriangle._make_db_key(0, 1), '(0,1)')
        
    def test_another(self):
        self.assertEquals(LieTriangle._make_db_key(100,99), '(100,99)')

class Binomial(unittest.TestCase):

    def setUp(self):
        self._choose = LieTriangle.binomial

    def _choose_old(self, n, k):
        if n < 0: raise ValueError
        if k < 0 or k > n: raise ValueError
        b = [0] * (n + 1)
        b[0] = 1
        for i in xrange(1, n + 1):
            b[i] = 1
            j = i - 1
            while j > 0:
                b[j] += b[j - 1]
                j -= 1
        return b[k]

    def test_neg(self):
        try:
            self._choose(-1, 0)
        except ValueError:
            pass
        else:
            assert 0, 'negative first arg to Binomial should raise'

    def test_zero_choose_zero(self):
        self.assertEquals(self._choose(0, 0), 1) #?check?

    def test_zero_choose_negative(self):
        try:
            self._choose(0, -1)
        except ValueError:
            pass
        else:
            assert 0, 'choosing negative in Binomial should raise'

    def test_zero_choose_positive(self):
        try:
            self._choose(0, 1)
        except ValueError:
            pass
        else:
            assert 0, 'choosing positive in Binomial zero should raise'

    def test_values(self):
        examples = {(0, 0): 1,
                    (1, 0): 1,
                    (1, 1): 1,
                    (2, 0): 1,
                    (2, 1): 2,
                    (2, 2): 1,
                    (4, 2): 6,
                    (9, 8): 9}
        for (n, k), b in examples.iteritems():
            self.assertEquals(self._choose(n, k), b)

    def test_identity0(self):
        for n in xrange(10):
            lhs = self._choose(n, 0)
            rhs = self._choose(n, n)
            self.assert_(lhs == rhs == 1)

    def test_identity1(self):
        for n in xrange(10):
            for k in xrange(n+1):
                lhs = self._choose(n, k)
                rhs = self._choose(n, n-k)
                self.assert_(lhs == rhs)

    def test_identity2_does_not_work_due_to_negatives(self):
        for n in xrange(10):
            for k in xrange(n+1):
                lhs = self._choose(n, k)*(-1 ** k)
                try:
                    rhs = self._choose(k-n-1, k)
                except ValueError:
                    pass
                else:
                    self.assert_(lhs == rhs)

    def test_identity3(self):
        for n in xrange(10):
            for k in xrange(n+1-1):
                lhs = self._choose(n, k+1)
                rhs = self._choose(n, k)*(n-k)/(k+1)
                self.assert_(lhs == rhs)

    def test_identity4(self):
        for n in xrange(10):
            for k in xrange(1, n+1):
                lhs = self._choose(n+1, k)
                rhs = self._choose(n, k)+self._choose(n, k-1)
                self.assert_(lhs == rhs)

    def test_cf_implementations(self):
        for i in xrange(1, 10):
            for j in xrange(1, i+1):
                b0 = self._choose_old(i, j)
                b1 = self._choose(i, j)
                self.assertEquals(b0, b1, (i,j,b0,b1))

    def test_range_downward_misses_final_value(self):
        c = 5
        for i in range(5, 0, -1):
            c -= 1
        self.assertEquals(c, 0)

    def test_xrange_downward_misses_final_value(self):
        c = 5
        for i in xrange(5, 0, -1):
            c -= 1
        self.assertEquals(c, 0)

class BlankLieTriangle(unittest.TestCase):

    """Test what is necessary to make a blank normal form, so that we
    can do a basic test of the methods."""

    def _create_nf(self):
        dof = 3
        steps = 2
        alg = LieAlgebra(dof)
        h = alg.zero()
        w = []
        h_db_name = os.path.join(test_dir, 'tmp/_test_h.db')
        k_db_name = os.path.join(test_dir, 'tmp/_test_k.db')
        self.freqs = [-1, +3J, +0.123456]
        h += alg.polynomial({Powers((1,1,0,0,0,0)): self.freqs[0],
                             Powers((0,0,1,1,0,0)): self.freqs[1],
                             Powers((0,0,0,0,1,1)): self.freqs[2]})
        h_2 = h
        state = (h_db_name, w, -1)
        normalizer = LieTriangle.LieTriangle(alg, h_2, state)
        return normalizer

    def test_create_nf(self):
        normalizer = self._create_nf()

class SolveHomologicalEquation(BlankLieTriangle):

    """Solution of the homological equation from a given residue
    (non-normalised terms)."""

    def setUp(self):
        self.normalizer = self._create_nf()

    def test_equal_powers_raises_zero_division_error(self):

        """Equal powers on any canonically-conjugate pair should at
        the very least raise a zero-division error."""
        
        rem = 0.5*Polynomial.Monomial((1, 1, 0, 0, 0, 0))
        try:
            gen = self.normalizer.solve_homological_eqn(rem)
        except ZeroDivisionError:
            pass
        else:
            assert 0, 'equal powers in homological equation should raise'

    def test_powers_differing_by_unity_induce_sum(self):
        norm = self.normalizer
        c = 0.5
        for mon in [Polynomial.Monomial((1, 2, 3, 4, 5, 6)),
                    Polynomial.Monomial((1, 2, 1, 2, 1, 2)),
                    Polynomial.Monomial((0, 1, 0, 1, 0, 1))]:
            rem = Polynomial(6)
            rem += c*mon
            gen = norm.solve_homological_eqn(rem)
            coeff = -c/sum(self.freqs)
            expected = Polynomial(6)+coeff*mon
            self.assert_((gen - expected).l1_norm() < 1.0e-15, (gen, expected))

class HillBase:
    
    def hill_about_l1(self):
        in_name = os.path.join(test_dir, 'ma-files/hill_l1_18_equi_to_tham.pol')
        in_file = open(in_name, 'r')
        h = read_ascii_polynomial(in_file, is_xxpp_format=1)
        in_file.close()
        terms = {Powers((2, 0, 0, 0, 0, 0)): -4.0, #x^2
                 Powers((0, 0, 2, 0, 0, 0)): +2.0, #y^2
                 Powers((0, 0, 0, 0, 2, 0)): +2.0, #z^2
                 Powers((0, 2, 0, 0, 0, 0)): +0.5, #px^2
                 Powers((0, 0, 0, 2, 0, 0)): +0.5, #py^2
                 Powers((0, 0, 0, 0, 0, 2)): +0.5, #pz^2
                 Powers((1, 0, 0, 1, 0, 0)): -1.0, #xpy
                 Powers((0, 1, 1, 0, 0, 0)): +1.0} #ypx
        assert len(terms) == 8
        dof = 3
        alg = LieAlgebra(dof)
        h_2 = alg.polynomial(terms)
        self.assert_(h_2 == alg.isograde(h, 2))
        return h

    def diagonalize(self, dof, h):
        tolerance = 5.0e-15

        lie = self.alg
        diag = Diagonalizer(lie)
        eig = diag.compute_eigen_system(h, tolerance=tolerance)
        diag.compute_diagonal_change()

        mat = diag.get_matrix_diag_to_equi()
        assert diag.matrix_is_symplectic(mat)
        
        sub_diag_into_equi = diag.matrix_as_vector_of_row_polynomials(mat)
        h_diag = h.substitute(sub_diag_into_equi)

        eq_type = eig.get_equilibrium_type()
        comp = Complexifier(lie, eq_type)
        self.sub_c_into_r = comp.calc_sub_complex_into_real()
        self.sub_r_into_c = comp.calc_sub_real_into_complex()
        h_comp = h_diag.substitute(self.sub_c_into_r)

        h_comp = h_comp.with_small_coeffs_removed(tolerance)
        self.assert_(lie.is_diagonal_polynomial(h_comp.homogeneous(2)))

        return h_comp

    def truncate_and_print_list_polynomials(self, hl, tolerance=1.0e-12):
        for p in hl:
            print p.with_small_coeffs_removed(tolerance)

    def setUp(self):
        self.dof = 3
        self.alg = LieAlgebra(self.dof) #classical
        self.h = self.alg.isograde(self.hill_about_l1(), 0, 6)
        self.h_dc = self.diagonalize(self.dof, self.h)

class NormalizeAlreadyNormalized(HillBase, unittest.TestCase):

    def test_hill_quadratic_invariant(self):
        steps = 10
        k_poly = self.alg.zero()
        h_dc = self.alg.isograde(self.h_dc, 2)
        nf = LieTriangle.LieTriangle(self.alg, h_dc)
        for i in xrange(steps):
            h_term = self.alg.isograde(h_dc, i + 2)
            k_term, w_term = nf.compute_normal_form_and_generating_function(h_term)
            if i == 0:
                self.assert_(k_term == self.alg.isograde(h_dc, 2))
            else:
                self.assert_(k_term == self.alg.zero())
            h_dict, w_list, row = nf.get_state()
            self.assert_(w_list[i] == self.alg.zero())

    def test_hill_quadratic_invariant_in_steps(self):
        steps = 10
        k_poly = self.alg.zero()
        h_dc = self.alg.isograde(self.h_dc, 2)
        state = None
        for i in xrange(steps):
            previous_row = i - 1
            nf = LieTriangle.LieTriangle(self.alg, h_dc, state)
            h_term = self.alg.isograde(h_dc, i + 2)
            k_term, w_term = nf.compute_normal_form_and_generating_function(h_term)
            if i == 0:
                self.assert_(k_term == h_dc)
            else:
                self.assert_(k_term == self.alg.zero())
            state = nf.get_state()
            h_dict, w_list, row = state
            self.assert_(w_list[i] == self.alg.zero())

class NormalizeHillThreeDof(HillBase, unittest.TestCase):
    
    def test_hill_cubic(self):
        steps = 5
        k_list = []
        h_dc = self.alg.isograde(self.h_dc, 2, 4)
        h_dc_2 = self.alg.isograde(h_dc, 2)
        k_comp = self.alg.zero()
        w_comp = self.alg.zero()
        nf = LieTriangle.LieTriangle(self.alg, h_dc_2)
        for i in xrange(0, steps + 1):
            h_term = self.alg.isograde(h_dc, i + 2)
            k_term, w_term = nf.compute_normal_form_and_generating_function(h_term)
            k_comp += k_term
            w_comp += w_term
        k_real = k_comp.substitute(self.sub_r_into_c)
        self.assert_(k_real.imag().l_infinity_norm() < 1.0e-14)
        k_real = k_real.real()
        #print k_real

    def test_hill_cubic_cf_in_steps(self):
        last_step = 4
        h_dc = self.alg.isograde(self.h_dc, 2, 4)
        h_dc_2 = self.alg.isograde(h_dc, 2)

        #perform the computation in one call of the Lie triangle:
        k_poly = self.alg.zero()
        nf = LieTriangle.LieTriangle(self.alg, h_dc_2)
        for i in xrange(0, last_step + 1):
            h_term = self.alg.isograde(h_dc, i + 2)
            k_term, w_term = nf.compute_normal_form_and_generating_function(h_term)
            self.assert_(self.alg.is_isograde(k_term, i + 2))
            self.assert_(self.alg.is_isograde(w_term, i + 2))
            if i == 0:
                self.assert_(k_term == h_dc_2)
            k_poly += k_term
        h_dict, w_list, row = nf.get_state()
        self.assert_(row == last_step)

        # now, perform the computation in steps, i.e., repeatedly make
        # new triangle loaded with dict and list data and attempt to
        # continue:

        steps_k_poly = self.alg.zero()
        state = None
        for i in xrange(0, last_step + 1):
            steps_nf = LieTriangle.LieTriangle(self.alg, h_dc_2, state)
            h_term = self.alg.isograde(h_dc, i + 2)
            k_term, w_term = steps_nf.compute_normal_form_and_generating_function(h_term)
            self.assert_(self.alg.is_isograde(k_term, i + 2))
            self.assert_(self.alg.is_isograde(w_term, i + 2))
            if i == 0:
                self.assert_(k_term == h_dc_2)
            steps_k_poly += k_term
            state = steps_nf.get_state()
        steps_h_dict, steps_w_list, steps_row = state
        self.assert_(steps_row == last_step)

        self.assert_(h_dict == steps_h_dict)
        self.assert_(w_list == steps_w_list)

        #compare the two methods:
        l1_tol = 1.0e-15
        for i in xrange(last_step + 1):
            x, y = (w_list[i], steps_w_list[i])
            self.assert_(self.alg.grade(x) == self.alg.grade(y))
            err = (x - y).l1_norm()
            self.assert_(err < l1_tol, err)
            err = (self.alg.isograde(k_poly, i+2) - self.alg.isograde(steps_k_poly, i+2)).l1_norm()
            self.assert_(err < l1_tol, err)

def suite():
    suites = []
    suites.append(unittest.makeSuite(NormalizeHillThreeDof))
    suites.append(unittest.makeSuite(SumPolynomialList))
    suites.append(unittest.makeSuite(NormalizeAlreadyNormalized))
    suites.append(unittest.makeSuite(Binomial))
    suites.append(unittest.makeSuite(SolveHomologicalEquation))
    suites.append(unittest.makeSuite(BlankLieTriangle))
    suites.append(unittest.makeSuite(DbKeys))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
