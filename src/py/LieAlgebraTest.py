#Author: Dr. Andrew David Burbanks, 2005.
# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
import random

from LieAlgebra import *
from Powers import Powers
from Polynomial import Polynomial
## Automatically adapted for numpy.oldnumeric Dec 16, 2008 by alter_code1.py
try:
    from numpy.oldnumeric.mlab import fabs
except:
    from MLab import fabs

from PolynomialTest import rand_poly
from PolynomialTest import _eq as _poly_eq

def _eq(a, b, tol=1.0e-8):
    """

    We need to loosen tolerance for Moyal bracket calculations.

    """
    return _poly_eq(a, b, tol)

n_cases = 4

class Diagonal(unittest.TestCase):

    def test_not_diagonal(self):

        """Ensure that the is_diagonal test can return false by giving
        a non-diagonal example."""
        
        dof = 3
        terms = {Powers((2, 0, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 1, 0, 1, 1, 0)): 7.2,
                 Powers((0, 1, 0, 0, 1, 0)): 7.12,
                 Powers((0, 0, 3, 0, 1, 0)): -4.0 }
        h = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        self.assert_(not alg.is_diagonal_polynomial(h))

    def test_off_diagonal_but_equal_powers(self):

        """Test the is_diagonal check on a non-diagonal polynomial in
        which each coordinate and its associated momentum are raised
        to the same power."""
        
        dof = 3
        terms = {Powers((2, 2, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 0, 1, 1, 1, 1)): 7.2,
                 Powers((0, 0, 0, 0, 1, 1)): 7.12 }
        g = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        self.assert_(not alg.is_diagonal_polynomial(g))

    def test_diagonal_part(self):
        dof = 3
        terms = {Powers((2, 2, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 0, 3, 3, 0, 0)): 7.2,
                 Powers((0, 0, 0, 0, 1, 1)): 7.12 }
        g = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        self.assert_(alg.diagonal_part_of_polynomial(g) == g)
        h = 0.54*Polynomial.Monomial((2, 0, 1, 0, 0, 0))
        g += h
        self.assert_(alg.diagonal_part_of_polynomial(g)+h == g)

    def test_is_diagonal(self):

        """Test the is_diagonal check on a diagonal polynomial."""
        
        dof = 3
        terms = {Powers((2, 2, 0, 0, 0, 0)): -0.3,
                 Powers((1, 1, 0, 0, 0, 0)): 0.33,
                 Powers((0, 0, 3, 3, 0, 0)): 7.2,
                 Powers((0, 0, 0, 0, 1, 1)): 7.12 }
        g = Polynomial(2*dof, terms=terms)
        alg = LieAlgebra(dof)
        self.assert_(alg.is_diagonal_polynomial(g))

class ClassicalToSemiclassicalTest(unittest.TestCase):

    def test_example_conversion_via_substitution(self):
        #make a classical algebra
        alg = LieAlgebra(2)
        
        #two classical polynomials
        q0, p0 = alg.q(0), alg.p(0)
        q1, p1 = alg.q(1), alg.p(1)
        a = q0**2 + 2.0*q0*p0 + q1*p1
        b = q0 + q0*p0**3 + p1
        alg.check_elt(a)
        alg.check_elt(b)

        #form the poisson bracket and ensure non-zero
        pb_a_b = alg.poisson_bracket(a, b)
        self.assert_(pb_a_b)

        #make the corresponding semiclassical algebra
        sem = SemiclassicalLieAlgebra(alg.dof())

        #make a conversion substitution
        destination = [sem.coordinate_monomial(i) for i in xrange(alg.n_vars())]

        #conver the two polynomials
        sa = a.substitute(destination)
        sb = b.substitute(destination)
        sem.check_elt(sa)
        sem.check_elt(sb)
        self.assert_(len(sa) == len(a))
        self.assert_(len(sb) == len(b))

        #form the moyal bracket
        mb_a_b = sem.moyal_bracket(sa, sb)

        #convert the earlier poisson bracket to moyal and check
        pb_a_b_convert = pb_a_b.substitute(destination)
        self.assertEquals(mb_a_b, pb_a_b_convert)

    def test_example_conversion_via_class(self):
        #make a classical algebra
        cla = LieAlgebra(2)
        
        #two classical polynomials
        q0, p0 = cla.q(0), cla.p(0)
        q1, p1 = cla.q(1), cla.p(1)
        a = q0**2 + 2.0*q0*p0 + q1*p1
        b = q0 + q0*p0**3 + p1
        cla.check_elt(a)
        cla.check_elt(b)

        #form the poisson bracket and ensure non-zero
        pb_a_b = cla.poisson_bracket(a, b)
        self.assert_(pb_a_b)

        #make converter
        cla_to_sem = ClassicalToSemiclassical(cla)
        sem = cla_to_sem.semi_classical_algebra()

        #conver the two polynomials
        sa = cla_to_sem(a)
        sb = cla_to_sem(b)
        sem.check_elt(sa)
        sem.check_elt(sb)
        self.assert_(len(sa) == len(a))
        self.assert_(len(sb) == len(b))

        #form the moyal bracket
        mb_a_b = sem.moyal_bracket(sa, sb)

        #convert the earlier poisson bracket to moyal and check
        pb_a_b_convert = cla_to_sem(pb_a_b)
        self.assertEquals(mb_a_b, pb_a_b_convert)

class PoissonMoyal(unittest.TestCase):

    def test_poisson_bracket(self):
        alg = LieAlgebra(2)
        q0, p0 = alg.q(0), alg.p(0)
        q1, p1 = alg.q(1), alg.p(1)
        a = q0**2 + 2.0*q0*p0 + q1*p1
        b = q0 + q0*p0**3 + p1
        print alg.bracket(a, b), 'test here please'

    def test_moyal_bracket(self):
        alg = SemiclassicalLieAlgebra(2)
        q0, p0 = alg.q(0), alg.p(0)
        q1, p1 = alg.q(1), alg.p(1)
        hbar = alg.h_bar()
        a = q0**2 + 2.0*q0*p0 + q1*p1
        b = q0 + q0*p0**3 + p1*hbar
        print alg.bracket(a, b), 'test here please'

    def test_moyal_bracket_via_product(self):
        alg = SemiclassicalLieAlgebra(2)
        q0, p0 = alg.q(0), alg.p(0)
        q1, p1 = alg.q(1), alg.p(1)
        hbar = alg.h_bar()
        a = q0**2 + 2.0*q0*p0 + q1*p1
        b = q0 + q0*p0**3 + p1*hbar
        moy_ab = alg.moyal_bracket(a, b)
        moy_ab_prod = alg.moyal_bracket_via_product(a, b)
        print moy_ab
        print moy_ab_prod
        print 'test here please'
        self.assertEquals(moy_ab_prod, moy_ab)

class MoyalBracketIsLieBracket(unittest.TestCase):

    def setUp(self):
        self.n_cases = n_cases

    def _gen_arguments(self, alg):
        vars = alg.n_vars()
        a = rand_poly(vars = vars, max_terms=5, max_power=4)
        b = rand_poly(vars = vars, max_terms=5, max_power=4)
        c = rand_poly(vars = vars, max_terms=5, max_power=4)
        yield (a, b, c)

    def _test_identity(self, lhs, rhs):
        for case in xrange(self.n_cases):
            dof = 2*random.randint(1, 4)
            alg = SemiclassicalLieAlgebra(dof)
            vars = alg.n_vars()
            for a, b, c in self._gen_arguments(alg):
                left = lhs(alg, a, b, c)
                right = rhs(alg, a, b, c)
                self.assert_(_eq(left, right), 
                             '\nLHS:%s\n\nRHS:%s\n'%(left, right))

    def test_moyal_linear_in_first_mul(self):
        s = 2.0
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(s*a, b)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return s*mb(a, b)
        self._test_identity(lhs, rhs)

    def test_moyal_linear_in_first_mul_complex(self):
        s = 0.5+2.0J
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(s*a, b)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return s*mb(a, b)
        self._test_identity(lhs, rhs)

    def test_moyal_linear_in_first_add(self):
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a+b, c)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, c)+mb(b, c)
        self._test_identity(lhs, rhs)

    def test_moyal_linear_in_second_mul(self):
        s = 0.5
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, s*b)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return s*mb(a, b)
        self._test_identity(lhs, rhs)

    def test_moyal_linear_in_second_mul_complex(self):
        s = 0.5+2.0J
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, s*b)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return s*mb(a, b)
        self._test_identity(lhs, rhs)

    def test_moyal_linear_in_second_add(self):
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, b+c)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, b)+mb(a, c)
        self._test_identity(lhs, rhs)

    def test_moyal_reverse(self):
        def lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return  mb(a, b)
        def rhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return -mb(b, a)
        self._test_identity(lhs, rhs)

    def test_zzz_jacobi_identity(self):
        
        """Likely to fail spectacularly with many variables!"""
        
        def jacobi_lhs(alg, a, b, c):
            mb = alg.moyal_bracket
            return mb(a, mb(b, c))+mb(b, mb(c, a))+mb(c, mb(a, b))
        for case in xrange(1+int(n_cases/2)):
            dof = 2*random.randint(1, 4)
            alg = SemiclassicalLieAlgebra(dof)
            vars = alg.n_vars()
            a = rand_poly(vars, max_power=4, max_terms=5, max_coeff=0.5)
            b = rand_poly(vars, max_power=4, max_terms=5, max_coeff=0.5)
            c = rand_poly(vars, max_power=4, max_terms=5, max_coeff=0.5)
            self.assert_(_eq(jacobi_lhs(alg, a, b, c), alg.zero(), tol=5.0e-8))

    def test_moyal_self(self):
        for case in xrange(self.n_cases):
            dof = random.randint(1, 5)
            alg = SemiclassicalLieAlgebra(dof)
            vars = alg.n_vars()
            a = rand_poly(vars)
            self.assert_(_eq(alg.moyal_bracket(a, a), alg.zero(), tol=1.0e-8))

def suite():
    suites = []
    suites.append(unittest.makeSuite(ClassicalToSemiclassicalTest))
    do_semiclassical = False
    if do_semiclassical:
        suites.append(unittest.makeSuite(PoissonMoyal))
        suites.append(unittest.makeSuite(MoyalBracketIsLieBracket))
    else:
        print 'Semiclassical unit tests supressed for speed [see LieAlgebraTest.suite]'
    suites.append(unittest.makeSuite(Diagonal))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
