# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
## Automatically adapted for numpy.oldnumeric Dec 16, 2008 by alter_code1.py
try:
    import numpy.oldnumeric.mlab as MLab
    import numpy.oldnumeric.linear_algebra as LinearAlgebra
    from numpy.oldnumeric import matrixmultiply
except:
    import MLab, LinearAlgebra
    from Numeric import matrixmultiply
from math import sqrt
from Diagonal import Diagonalizer, Complexifier
from Powers import Powers
from Polynomial import Polynomial
from LieAlgebra import LieAlgebra

class HamiltonsEquations(unittest.TestCase):

    def test_skew_symmetric_diagonal_2x_2(self):

        """Test the basic structure of the skew-symmetric matrix."""
        
        dof = 3
        lie = LieAlgebra(dof)
        diag = Diagonalizer(lie)
        J = diag.skew_symmetric_matrix()
        for qi in xrange(0, dof, 2):
            pi = qi+1
            self.assertEquals(J[qi,qi],  0.0)
            self.assertEquals(J[qi,pi], +1.0)
            self.assertEquals(J[pi,qi], -1.0)
            self.assertEquals(J[pi,pi],  0.0)

    def test_skew_symmetric_off_diagonal_2x_2(self):
        
        """Test the basic structure of the skew-symmetric matrix."""
        
        dof = 3
        lie = LieAlgebra(dof)
        diag = Diagonalizer(lie)
        J = diag.skew_symmetric_matrix()
        for qi in xrange(0, dof, 2):
            pi = qi+1
            for j in xrange(0, dof):
                if j!=qi and j!=pi:
                    self.assertEquals(J[qi,j], 0.0)
                    self.assertEquals(J[pi,j], 0.0)
                    self.assertEquals(J[j,qi], 0.0)
                    self.assertEquals(J[j,pi], 0.0)

    def test_example_skew(self):

        """Multiply the skew-symmetric matrix by a specific vector."""
        
        dof = 3
        lie = LieAlgebra(dof)
        diag = Diagonalizer(lie)
        J = diag.skew_symmetric_matrix()
        x = (1,2,3,4,5,6)
        y = matrixmultiply(J, x)
        # self.assertEquals(y, (2,-1,4,-3,6,-5))
        z = y == (2,-1,4,-3,6,-5)
        self.assertTrue(z.all())

    def test_cf_hamilton_linear_to_quadratic_matrix(self):

        """Compare the results of (1) asking directly for the linear
        matrix for a given Hamiltonian, and (2) forming Hamilton's
        equations and then linearizing them."""
        
        dof = 3
        terms = {Powers((2, 0, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 1, 0, 1, 1, 0)): 7.2,
                 Powers((0, 1, 0, 0, 1, 0)): 7.12,
                 Powers((0, 0, 3, 0, 1, 0)): -4.0 }
        h = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        diag = Diagonalizer(alg)
        lin = diag.linear_matrix(h)
        rhs = diag.hamiltons_equations_rhs(h)
        rhs_lin = tuple([alg.isograde(h_term, 1) for h_term in rhs])
        for i, pol in enumerate(rhs_lin):
            for m, c in pol.powers_and_coefficients():
                mon = alg.monomial(m, c)
                self.assertEquals(alg.grade(mon), 1)
                for j, f in enumerate(m):
                    if f==1:
                        self.assertEquals(lin[i, j], c)

    def test_eigensystem(self):

        """Multiply eigenvectors and the original matrix by the
        eigenvalues in order to check the integrity of the
        eigensystem."""
        
        dof = 3
        terms = {Powers((2, 0, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 0, 1, 0, 1, 0)): 7.2,
                 Powers((0, 0, 0, 0, 0, 2)): 7.2,
                 Powers((0, 0, 0, 1, 1, 0)): 7.12 }
        g = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        diag = Diagonalizer(alg)
        lin = diag.linear_matrix(g)
        val_vec_pairs = diag.eigenvalue_eigenvector_pairs(g)
        for p in val_vec_pairs:
            prod_s = p.vec*p.val
            prod_v = matrixmultiply(lin, p.vec)
            for x in prod_s-prod_v:
                self.assert_(abs(x)<1.0e-15)

    def test_matrix_as_polynomials(self):
        dof = 1
        mat = ((1, 2), (3, 4))
        alg = LieAlgebra(dof)
        diag = Diagonalizer(alg)
        vp = diag.matrix_as_vector_of_row_polynomials(mat)
        self.assert_(vp[0]((0.0, 0.0)) == 0.0)
        self.assert_(vp[0]((1.0, 0.0)) == 1.0)
        self.assert_(vp[0]((1.0, 1.0)) == 3.0)
        self.assert_(vp[0]((0.0, 1.0)) == 2.0)
        self.assert_(vp[1]((0.0, 0.0)) == 0.0)
        self.assert_(vp[1]((1.0, 0.0)) == 3.0)
        self.assert_(vp[1]((1.0, 1.0)) == 7.0)
        self.assert_(vp[1]((0.0, 1.0)) == 4.0)

class CanonicalChange(unittest.TestCase):

    def setUp(self):

        """Set up an example from Hill's equations."""

        x_2 = Powers((2, 0, 0, 0, 0, 0))
        px2 = Powers((0, 2, 0, 0, 0, 0))
        y_2 = Powers((0, 0, 2, 0, 0, 0))
        py2 = Powers((0, 0, 0, 2, 0, 0))
        z_2 = Powers((0, 0, 0, 0, 2, 0))
        pz2 = Powers((0, 0, 0, 0, 0, 2))
        xpy = Powers((1, 0, 0, 1, 0, 0))
        ypx = Powers((0, 1, 1, 0, 0, 0))
        terms = {px2: 0.5, py2: 0.5, pz2: 0.5,
                 xpy: -1.0, ypx: 1.0,
                 x_2: -4.0, y_2: 2.0, z_2: 2.0}
        assert len(terms) == 8

        self.h_2 = Polynomial(6, terms=terms)
        self.lie = LieAlgebra(3)
        self.diag = Diagonalizer(self.lie)
        self.eq_type = 'scc'
        e = []
        e.append(+sqrt(2.0*sqrt(7.0)+1.0))
        e.append(-e[-1])
        e.append(complex(0.0, +sqrt(2.0*sqrt(7.0)-1.0)))
        e.append(-e[-1])
        e.append(complex(0.0, +2.0))
        e.append(-e[-1])
        self.eig_vals = e

    def test_basic_matrix(self):

        """This test is rather monolithic, but at least it implements
        a concrete example that we can compare with our earlier
        computations.  It also tests the mutual-inverse character of
        the equi-to-diag and diag-to-equi transformations."""

        tolerance = 5.0e-15
        eig = self.diag.compute_eigen_system(self.h_2, tolerance)
        self.diag.compute_diagonal_change()

        eq_type = eig.get_equilibrium_type()
        self.assertEquals(eq_type, self.eq_type)

        eigs = [pair.val for pair in eig.get_raw_eigen_value_vector_pairs()]
        for actual, expected in zip(eigs, self.eig_vals):
            self.assert_(abs(actual-expected) < tolerance, (actual, expected))

        mat = self.diag.get_matrix_diag_to_equi()
        assert self.diag.matrix_is_symplectic(mat)

        sub_diag_into_equi = self.diag.matrix_as_vector_of_row_polynomials(mat)
        mat_inv = LinearAlgebra.inverse(MLab.array(mat))
        sub_equi_into_diag = self.diag.matrix_as_vector_of_row_polynomials(mat_inv)
        h_diag_2 = self.h_2.substitute(sub_diag_into_equi)
        h_2_inv = h_diag_2.substitute(sub_equi_into_diag)
        self.assert_(h_2_inv) #non-zero
        self.assert_(not h_2_inv.is_constant())
        self.assert_(self.lie.is_isograde(h_2_inv, 2))
        self.assert_((self.h_2-h_2_inv).l1_norm() < 1.0e-14)
        comp = Complexifier(self.diag.get_lie_algebra(), eq_type)
        sub_complex_into_real = comp.calc_sub_complex_into_real()
        h_comp_2 = h_diag_2.substitute(sub_complex_into_real)
        h_comp_2 = h_comp_2.with_small_coeffs_removed(tolerance)
        self.assert_(self.lie.is_diagonal_polynomial(h_comp_2))

class RealVsComplex(unittest.TestCase):

    def test_saddle(self):
        lie = LieAlgebra(1)
        diag = Diagonalizer(lie)
        q = lie.q
        p = lie.p
        orig = q(0)*p(0)
        eq_type = 's'
        comp = Complexifier(lie, eq_type)
        r2c = comp.calc_sub_complex_into_real()
        c2r = comp.calc_sub_real_into_complex()
        comp = orig.substitute(r2c)
        real = comp.substitute(c2r)
        diff = orig-real
        for m, c in diff.powers_and_coefficients():
            self.assert_(abs(c)<1.0e-15)

    def test_centre(self):
        lie = LieAlgebra(1)
        diag = Diagonalizer(lie)
        q = lie.q
        p = lie.p
        orig = (0.5*q(0)**2)+(0.5*p(0)**2)
        self.assertEquals(len(orig), 2)
        eq_type = 'c'
        comp = Complexifier(lie, eq_type)
        r2c = comp.calc_sub_complex_into_real()
        c2r = comp.calc_sub_real_into_complex()
        comp = orig.substitute(r2c)
        self.assertEquals(len(comp), 1)
        for m, c in comp.powers_and_coefficients():
            self.assert_(c.real == 0.0)
            self.assert_(m[0] == 1)
            self.assert_(m[1] == 1)
        real = comp.substitute(c2r)
        diff = orig-real
        for m, c in diff.powers_and_coefficients():
            self.assert_(abs(c)<1.0e-15)

    def test_centre_saddle(self):
        lie = LieAlgebra(2)
        diag = Diagonalizer(lie)
        q = lie.q
        p = lie.p
        cent = (0.5*q(0)**2)+(0.5*p(0)**2)
        sadd = (1.0*q(1))*(1.0*p(1))
        orig = cent + sadd
        self.assertEquals(len(orig), 3)
        eq_type = 'cs'
        comp = Complexifier(lie, eq_type)
        r2c = comp.calc_sub_complex_into_real()
        c2r = comp.calc_sub_real_into_complex()
        comp = orig.substitute(r2c)
        self.assertEquals(len(comp), 2)
        real = comp.substitute(c2r)
        diff = orig-real
        for m, c in diff.powers_and_coefficients():
            self.assert_(abs(c)<1.0e-15)

    def test_mutual_inverses(self):
        lie = LieAlgebra(6)
        diag = Diagonalizer(lie)
        pol = lie.zero()
        for i in xrange(6):
            pol += lie.q(i)
            pol += lie.p(i)
        eq_type = 'scccsc'
        comp = Complexifier(lie, eq_type)
        sub_c_into_r = comp.calc_sub_complex_into_real()
        sub_r_into_c = comp.calc_sub_real_into_complex()
        pol_c = pol.substitute(sub_c_into_r)
        pol_r = pol_c.substitute(sub_r_into_c)
        self.assert_(len(pol_r) == len(pol))
        for i in xrange(6):
            self.assert_((pol_r-pol).l_infinity_norm() < 1.0e-15)
            self.assert_((pol_r-pol).l_infinity_norm() < 1.0e-15)

def suite():
    suites = []
    suites.append(unittest.makeSuite(HamiltonsEquations))
    suites.append(unittest.makeSuite(CanonicalChange))
    suites.append(unittest.makeSuite(RealVsComplex))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
