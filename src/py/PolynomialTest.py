# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest
import random
from copy import deepcopy

from Polynomial import *
from NormalFormIO import read_ascii_polynomial, write_ascii_polynomial, xxpp_to_xpxp, xpxp_to_xxpp
from Powers import Powers

n_cases = 2 #32

tmp_path = '../../test/tmp/'

random.seed(54321)

class Comparison(unittest.TestCase):

    def test_self_equality(self):
        for vars in xrange(1, 50):
            p = Polynomial(vars)
            self.assert_(p == p)
            self.assertEquals(p, p)

    def test_equality_with_complex_coeffs(self):
        common_powers = Powers((1, 2, 3))
        coeff = 1.0-0.5j
        p = Polynomial(3)
        q = Polynomial(3)
        for pol in q, p:
            pol[common_powers] = coeff
        self.assert_(q == p)

class FromPowers(unittest.TestCase):

    def test_type(self):
        p = Polynomial.Monomial((1, 2, 3))
        self.assert_(isinstance(p, Polynomial))

    def test_get_item_exists(self):
        p = Polynomial.Monomial((1, 2, 3))
        self.assertEquals(p[Powers((1, 2, 3))], 1.0)

    def test_get_item_only_exists(self):
        p = Polynomial.Monomial((1, 2, 3))
        self.assertEquals(p[Powers((2, 3, 4))], 0.0)

    def test_len(self):
        p = Polynomial.Monomial((1, 2, 3))
        self.assertEquals(len(p), 1)

    def test_has_term(self):
        t = (1, 2, 3)
        p = Polynomial.Monomial(t)
        self.assert_(p.has_term(Powers(t)))
        self.assert_(not p.has_term(Powers((2, 2, 3))))

    def test_call(self):
        p = -1.5*Polynomial.Monomial((1, 2, 3))
        self.assertEquals(p((1.0, 1.0, 1.0)), -1.5)

    def test_call_negative_count(self):
        for v in xrange(1, 10):
            one = (1, )*v
            p = 2.0*Polynomial.Monomial(one)
            mone = (-1, )*v
            if v%2:
                self.assertEquals(p(mone), -2.0)
            else:
                self.assertEquals(p(mone), +2.0)

class Default(unittest.TestCase):

    """

    Test a default-constructed Polynomial.

    """

    def test_repr(self):
        t = Polynomial(3)
        s = 'Polynomial(3)'
        self.assertEquals(repr(t), s)

    def test_len(self):
        t = Polynomial(9)
        self.assertEquals(len(t), 0)

class Zero(unittest.TestCase):

    def test_poly_with_no_terms_is_constant(self):
        p = Polynomial(6)
        self.assert_(p.is_constant())

    def test_poly_with_no_terms_is_zero_degree(self):
        p = Polynomial(6)
        self.assertEquals(p.degree(), 0)

    def test_poly_with_no_terms_is_false(self):
        t = Polynomial(52)
        if t:
            assert 0, 'zero polynomial should act like False'
        self._assert_zero_poly(t)

    def test_poly_with_no_nonzero_terms_is_false(self):
        t = Polynomial(3, terms={Powers((3, 6, 4)): 0.0})
        if t:
            assert 0, 'polynomial should not store zero terms'
        self._assert_zero_poly(t)

    def test_poly_with_some_terms_is_true(self):
        t = Polynomial(3, terms={Powers((1, 2, 3)):0.1})
        if t:
            pass
        else:
            assert 0, 'non-zero polynomial should act like True'
        try:
            self._assert_zero_poly(t)
        except AssertionError:
            pass
        else:
            assert 0, 'non-zero polynomial not detected.'

    def _assert_zero_poly(self, p):
        self.assert_(isinstance(p, Polynomial))

        #we can use boolean test.  since Polynomial does not implement
        #__nonzero__, instead __len__ will be called, which gives the
        #number of non-zero terms and hence behaves correctly.

        self.assert_(not p)

    def test_mul_float(self):

        """Multiplication by zero should remove all monomials."""

        t = Polynomial.Monomial((1, 2, 3, 4, 5))
        self.assert_(isinstance(t, Polynomial))
        u = 0.0*t
        self._assert_zero_poly(u)

    def test_mul_other(self):
        t = Polynomial.Monomial((2, 2, 3, 4, 5))
        self.assert_(isinstance(t, Polynomial))
        self.assertEquals(len(t), 1)
        tt = Polynomial.Monomial((1, 2, 3, 4, 5))
        t += tt
        self.assert_(isinstance(t, Polynomial))
        self.assertEquals(len(t), 2, repr(t._terms))
        zero = 0.0*Polynomial.Monomial((5, 2, 8, 4, 5))
        self.assert_(isinstance(zero, Polynomial))
        t0 = t*zero
        self._assert_zero_poly(t0)
        t1 = zero*t
        self._assert_zero_poly(t1)
        zero = 0.0*tt
        self.assert_(isinstance(zero, Polynomial))
        t0 = t*zero
        self._assert_zero_poly(t0)
        t1 = zero*t
        self._assert_zero_poly(t1)

class Norms(unittest.TestCase):

    def test_l1_zero(self):
        p = Polynomial(3)
        self.assertEquals(p.l1_norm(), 0.0)

    def test_l1_example(self):
        p = Polynomial(3)
        p += 1.0*Polynomial.Monomial((1, 2, 3))
        self.assertEquals(p.l1_norm(), 1.0)
        p -= 2.0*Polynomial.Monomial((0, 2, 0))
        self.assertEquals(p.l1_norm(), 3.0)
        p += (0.0+5.0J)*Polynomial.Monomial((0, 1, 1))
        self.assertEquals(p.l1_norm(), 8.0)

class Degree(unittest.TestCase):

    def test_zero(self):
        p = Polynomial(6)
        self.assertEquals(p.degree(), 0)

    def test_one(self):
        for j in xrange(5):
            for i in xrange(6):
                pows = [0,]*6
                pows[i] = j
                p = Polynomial(6)
                p = 1.0*Polynomial.Monomial(tuple(pows))
                self.assertEquals(p.degree(), j)

    def test_example(self):
        p = Polynomial(5)
        p += 0.1*Polynomial.Monomial((1, 0, 1, 1, 0))
        p -= 0.2*Polynomial.Monomial((1, 7, 1, 1, 0))
        p += 0.3*Polynomial.Monomial((0, 0, 1, 1, 2))
        p -= 0.4*Polynomial.Monomial((1, 0, 9, 1, 0))
        self.assertEquals(p.degree(), 11)

class Call(unittest.TestCase):

    def setUp(self):
        self.powers = tuple(xrange(10))
        self.tol = 1.0e-16

    def test_mul(self):
        p = Polynomial.Monomial(self.powers)
        p_accumulate = Polynomial(len(self.powers))
        p_accumulate[Powers((0, )*len(self.powers))] = 1.0
        x = (0.5, -0.2, 0.3, 1.67, -2.001, 0.1, 0.2, 0.3, 0.4, -0.9)
        prod = 1.0
        for i in range(len(self.powers)):
            one = [0,]*len(self.powers)
            one[i] = self.powers[i]
            m_one = Powers(tuple(one))
            p_one = Polynomial.Monomial(tuple(one))
            prod *= p_one(x)
            p_accumulate = p_accumulate * p_one
        self.assert_(abs(p(x)-prod)<self.tol)
        self.assert_(abs(p_accumulate(x)-prod)<self.tol)

    def test_add(self):
        p_as_one = Polynomial.Monomial(self.powers)
        x = (0.5, -0.2, 0.3, 1.67, -2.001, 0.1, 0.2, 0.3, 0.4, -0.9)
        coeff = 1.0/float(len(self.powers))
        total = 0.0
        p_as_many = Polynomial(len(self.powers))
        p_accumulate = Polynomial(len(self.powers))
        p_accumulate[Powers((0, )*len(self.powers))] = 0.0
        for i in range(len(self.powers)):
            m_one = Powers(self.powers)
            p_one = coeff*Polynomial.Monomial(self.powers)
            total += p_one(x)
            p_as_many[m_one] = coeff
            p_accumulate = p_accumulate + p_one
        self.assert_(abs(p_as_one(x)-total)<self.tol)
        self.assert_(abs(p_as_many(x)-total)<self.tol)
        self.assert_(abs(p_accumulate(x)-total)<self.tol)

class Diff(unittest.TestCase):

    def test_diff_sum(self):
        for vars in xrange(1, 100):
            triangle = Powers(tuple(range(1, vars+1)))
            t = Polynomial(vars, terms={triangle: -1.0})
            u = Polynomial(vars)
            for i in xrange(vars):
                u += t.diff(i)
            self.assertEquals(len(u), vars)
            ones = (1.0, )*vars
            self.assertEquals(u(ones), -float((vars*(vars+1))/2))

    def test_diff_pow0(self):
        for vars in xrange(1, 100):
            triangle = Powers(tuple(range(1, vars+1)))
            t = Polynomial(vars, terms={triangle: -1.0})
            for i in xrange(vars):
                u = t.diff_pow(i, 0)
                self.assertEquals(u, t)

    def test_diff_pow1_sum(self):
        for vars in xrange(1, 100):
            triangle = Powers(tuple(range(1, vars+1)))
            t = Polynomial(vars, terms={triangle: -1.0})
            u = Polynomial(vars)
            for i in xrange(vars):
                u += t.diff_pow(i, 1)
            self.assertEquals(len(u), vars)
            ones = (1.0, )*vars
            self.assertEquals(u(ones), -float((vars*(vars+1))/2))

class Constant(unittest.TestCase):

    def test_example(self):
        p = Polynomial(6, terms={Powers((0,0,0,0,0,0)): -0.5})
        self.assert_(p.is_constant())

    def test_too_many_terms(self):
        p = Polynomial(6, terms={Powers((0,0,0,0,0,0)): -0.5,
                                 Powers((1,0,0,0,0,0)): 0.2})
        self.assert_(not p.is_constant())

    def test_non_constant_single_term(self):
        p = Polynomial(6, terms={Powers((0,0,0,2,0,0)): -0.5})
        self.assert_(not p.is_constant())

def rand_powers(vars, max_power):
    return tuple([random.randint(0, max_power) for i in xrange(vars)])

def rand_monomial(vars, max_power):
    return Powers(rand_powers(vars, max_power))

def rand_float(max_size=1.0):
    return max_size*(2.0*random.random()-1.0)

def rand_coeff(max_size=1.0):
    return complex(rand_float(max_size), rand_float(max_size))

def rand_poly(vars=-1, max_vars=40, max_power=10, max_terms=25, max_coeff=1.0):
    if vars <= 0:
        vars = random.randint(1, max_vars)
    t = Polynomial(vars)
    n_terms = random.randint(0, max_terms)
    for i in xrange(n_terms):
        m = rand_monomial(vars, max_power)
        t[m] = rand_coeff(max_coeff)
    return t

def zero(a, b, c):
    vars = a.n_vars()
    return Polynomial(vars)

def one(a, b, c):
    vars = a.n_vars()
    m = Powers((0, )*vars)
    t = Polynomial(vars)
    t[m] = 1.0
    return t

def _eq(a, b, tol=None):
    if tol==None:
        tol = 5.0e-14
    vars = a.n_vars()
    diff = a-b
    if not diff:
        return True
    for m, c in diff.powers_and_coefficients():
        if abs(c)>tol:
            print 'error coeff size:', abs(c)
            return False
    return True

class Ring(unittest.TestCase):

    def setUp(self):
        self.n_cases = n_cases

    def test_sub_for_use_in_eq(self):

        """Important that this is true, or all bets are off!"""

        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                z = zero(a, b, c)
                neg_a = z-a
                for m, c in a.powers_and_coefficients():
                    self.assertEquals(neg_a[m], -c)
                for m, c in neg_a.powers_and_coefficients():
                    self.assertEquals(a[m], -c)

    def _gen_arguments(self, vars):
        a = rand_poly(vars)
        b = rand_poly(vars)
        c = rand_poly(vars)
        yield (a, b, c)

    def _test_identity(self, lhs, rhs):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                left = lhs(a, b, c)
                right = rhs(a, b, c)
                self.assert_(_eq(left, right), 
                             '\nLHS:%s\n\nRHS:%s\n'%(left, right))

    def test_associative_add(self):
        def lhs(a, b, c): return (a+b)+c
        def rhs(a, b, c): return a+(b+c)
        self._test_identity(lhs, rhs)

    def test_associative_mul(self):
        def lhs(a, b, c): return (a*b)*c
        def rhs(a, b, c): return a*(b*c)
        self._test_identity(lhs, rhs)

    def test_commutative_add(self):
        def lhs(a, b, c): return a+b
        def rhs(a, b, c): return b+a
        self._test_identity(lhs, rhs)

    def test_commutative_mul(self):
        def lhs(a, b, c): return a*b
        def rhs(a, b, c): return b*a
        self._test_identity(lhs, rhs)

    def test_identity_add(self):
        def lhs(a, b, c): return a+zero(a, b, c)
        def rhs(a, b, c): return a
        self._test_identity(lhs, rhs)

    def test_identity_mul(self):
        def lhs(a, b, c): return a*one(a, b, c)
        def rhs(a, b, c): return a
        self._test_identity(lhs, rhs)

    def test_absorb_mul(self):
        def lhs(a, b, c): return a*zero(a, b, c)
        def rhs(a, b, c): return zero(a, b, c)
        self._test_identity(lhs, rhs)

    def test_inverse_add(self):
        def lhs(a, b, c): return a+(-a)
        def rhs(a, b, c): return zero(a, b, c)
        self._test_identity(lhs, rhs)

    def test_distrib_left(self):
        def lhs(a, b, c): return a*(b+c)
        def rhs(a, b, c): return (a*b)+(a*c)
        self._test_identity(lhs, rhs)

    def test_distrib_right(self):
        def lhs(a, b, c): return (a+b)*c
        def rhs(a, b, c): return (a*c)+(b*c)
        self._test_identity(lhs, rhs)

    def test_scalar_mult(self):
        def lhs(a, b, c): return 0.5*a+0.5*a
        def rhs(a, b, c): return a
        self._test_identity(lhs, rhs)

    def test_neg_mult(self):
        def lhs(a, b, c): return (-a)
        def rhs(a, b, c): return (-1.0)*a
        self._test_identity(lhs, rhs)

    def test_zero_neg(self):

        """Should ensure that terms which exist in the right operand
        of a subtraction but not in the left are handled correctly."""
        
        def lhs(a, b, c): return zero(a, b, c)-a
        def rhs(a, b, c): return (-1.0)*a
        self._test_identity(lhs, rhs)

    def test_add_sub(self):
        def lhs(a, b, c): return a-(b+c)
        def rhs(a, b, c): return (a-b)-c
        self._test_identity(lhs, rhs)

    def test_zero_add_sub(self):
        def lhs(a, b, c): return (zero(a, b, c)+a)-a
        def rhs(a, b, c): return zero(a, b, c)
        self._test_identity(lhs, rhs)

    def test_add_neg_sub(self):
        def lhs(a, b, c): return -(a+b+c)
        def rhs(a, b, c): return (-a)+(-b)+(-c)
        self._test_identity(lhs, rhs)

    def test_in_place_add(self):
        def lhs(a, b, c):
            d = a.copy()
            d += b
            return d
        def rhs(a, b, c):
            d = a+b
            return d
        self._test_identity(lhs, rhs)

    def test_in_place_add_self(self):
        def lhs(a, b, c):
            d = a.copy()
            d += d
            return d
        def rhs(a, b, c):
            d = 2.0*a
            return d
        self._test_identity(lhs, rhs)

    def test_in_place_sub(self):
        def lhs(a, b, c):
            d = a.copy()
            d -= b
            return d
        def rhs(a, b, c):
            d = a-b
            return d
        self._test_identity(lhs, rhs)

    def test_in_place_sub_self(self):
        def lhs(a, b, c):
            d = a.copy()
            d -= d
            return d
        def rhs(a, b, c):
            return zero(a, b, c)
        self._test_identity(lhs, rhs)

    def test_in_place_mul(self):
        def lhs(a, b, c):
            d = a.copy()
            d *= b
            return d
        def rhs(a, b, c):
            d = a*b
            return d
        self._test_identity(lhs, rhs)

    def test_in_place_mul_self(self):
        def lhs(a, b, c):
            d = a.copy()
            d *= d
            return d
        def rhs(a, b, c):
            d = a*a
            return d
        self._test_identity(lhs, rhs)

    def test_monomials_add(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                s = Polynomial(vars)
                for po, co in a.powers_and_coefficients():
                    mon = co*Polynomial.Monomial(tuple(po))
                    s += mon
                self.assert_(_eq(a, s))

    def test_monomials_set(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                s = Polynomial(vars)
                for po, co in a.powers_and_coefficients():
                    s[po] = co
                self.assert_(_eq(a, s))

    def test_homogeneous_add(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                s = Polynomial(vars)
                for deg in xrange(a.degree()+1):
                    h = a.homogeneous(deg)
                    s += h
                self.assert_(_eq(a, s))

    def test_call_sum(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                s = 0.0
                for m, c in a.powers_and_coefficients():
                    s += c
                ones = (1.0, )*vars
                self.assertEquals(a(ones), s)

    def test_has_terms(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 100)
            for a, b, c in self._gen_arguments(vars):
                for m, c in a.powers_and_coefficients():
                    self.assert_(a.has_term(m))

    def test_powers_and_coefficients(self):
        def lhs(a, b, c):
            d = Polynomial(a.n_vars())
            for m, c in a.powers_and_coefficients():
                self.assertEquals(len(m), a.n_vars())
                d[m] = c
            return d
        def rhs(a, b, c):
            return a
        self._test_identity(lhs, rhs)

    def test_partition_odds_and_evens(self):
        def is_odd_degree(m):
            return (m.degree())%2
        def lhs(a, b, c):
            odds, evens = a.partition(isOdd)
            for m, c in odds.powers_and_coefficients():
                self.assert_((m.degree())%2)
            for m, c in odds.powers_and_coefficients():
                self.assert_(not (m.degree())%2)
            d = odds+evens
            return d
        def rhs(a, b, c):
            return a

class Poisson(unittest.TestCase):

    def setUp(self):
        self.n_cases = max(1, int(n_cases/2))

    def _gen_arguments(self, vars):
        a = rand_poly(vars = vars, max_power=4)
        b = rand_poly(vars = vars, max_power=4)
        c = rand_poly(vars = vars, max_power=4)
        yield (a, b, c)

    def _test_identity(self, lhs, rhs):
        for case in xrange(self.n_cases):
            vars = 2*random.randint(1, 5)
            for a, b, c in self._gen_arguments(vars):
                left = lhs(a, b, c)
                right = rhs(a, b, c)
                self.assert_(_eq(left, right), 
                             '\nLHS:%s\n\nRHS:%s\n'%(left, right))

    def test_poisson_linear_in_first_mul(self):
        pb = Polynomial.poisson_bracket
        s = 2.0
        def lhs(a, b, c): return pb(s*a, b)
        def rhs(a, b, c): return s*pb(a, b)
        self._test_identity(lhs, rhs)

    def test_poisson_linear_in_first_mul_complex(self):
        pb = Polynomial.poisson_bracket
        s = 0.5+2.0J
        def lhs(a, b, c): return pb(s*a, b)
        def rhs(a, b, c): return s*pb(a, b)
        self._test_identity(lhs, rhs)

    def test_poisson_linear_in_first_add(self):
        pb = Polynomial.poisson_bracket
        def lhs(a, b, c): return pb(a+b, c)
        def rhs(a, b, c): return pb(a, c)+pb(b, c)
        self._test_identity(lhs, rhs)

    def test_poisson_linear_in_second_mul(self):
        pb = Polynomial.poisson_bracket
        s = 0.5
        def lhs(a, b, c): return pb(a, s*b)
        def rhs(a, b, c): return s*pb(a, b)
        self._test_identity(lhs, rhs)

    def test_poisson_linear_in_second_mul_complex(self):
        pb = Polynomial.poisson_bracket
        s = 0.5+2.0J
        def lhs(a, b, c): return pb(a, s*b)
        def rhs(a, b, c): return s*pb(a, b)
        self._test_identity(lhs, rhs)

    def test_poisson_linear_in_second_add(self):
        pb = Polynomial.poisson_bracket
        def lhs(a, b, c): return pb(a, b+c)
        def rhs(a, b, c): return pb(a, b)+pb(a, c)
        self._test_identity(lhs, rhs)

    def test_poisson_reverse(self):
        pb = Polynomial.poisson_bracket
        def lhs(a, b, c): return  pb(a, b)
        def rhs(a, b, c): return -pb(b, a)
        self._test_identity(lhs, rhs)

    def test_zzz_jacobi_identity(self):
        
        """Likely to fail spectacularly with many variables!"""
        
        pb = Polynomial.poisson_bracket
        def jacobi_lhs(a, b, c):
            return pb(a, pb(b, c))+pb(b, pb(c, a))+pb(c, pb(a, b))
        for case in xrange(8):
            vars = 2*random.randint(1, 8)
            a = rand_poly(vars, max_power=4, max_terms=25, max_coeff=0.5)
            b = rand_poly(vars, max_power=4, max_terms=25, max_coeff=0.5)
            c = rand_poly(vars, max_power=4, max_terms=25, max_coeff=0.5)
            self.assert_(_eq(jacobi_lhs(a, b, c), zero(a, b, c), tol=5.0e-8))

    def test_poisson_checks_even(self):
        for case in xrange(self.n_cases):
            vars = 1+2*random.randint(0, 25)
            a = rand_poly(vars)
            try:
                a.poisson_bracket(a)
            except IndexError:
                pass
            else:
                self.assert_(False)

    def test_poisson_self(self):
        for case in xrange(self.n_cases):
            vars = 2*random.randint(1, 50)
            a = rand_poly(vars)
            self.assert_(_eq(a.poisson_bracket(a), zero(a, a, a), tol=1.0e-12))

class PowersAndCoeff(unittest.TestCase):

    def setUp(self):
        self.n_cases = max(1, int(n_cases/2))

    def test_monomials_never_repeated(self):
        for case in xrange(self.n_cases):
            vars = 1+2*random.randint(0, 25)
            p = rand_poly(vars)
            found = []
            for m, c in p.powers_and_coefficients():
                for f in found:
                    assert not (f == m) #handles tuple/monomial/packed.
                found.append(m)

    def test_coeffs_never_zero(self):
        
        """This test might not be so convincing given our method for
        generating random polynomials.  Therefore, we set deliberately
        each term to zero and test again."""
        
        for case in xrange(self.n_cases):
            vars = 1+2*random.randint(0, 25)
            p = rand_poly(vars)
            for m, c in p.powers_and_coefficients():
                assert c != 0.0
                q = p.copy()
                q[m] = 0.0
                assert not (q.has_term(m))

class Pow(unittest.TestCase):

    def setUp(self):
        self.n_cases = max(1, int(n_cases/2))

    def test_power_zero_gives_one(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 5)
            p = rand_poly(vars)
            q = p**0
            one = Polynomial(p.n_vars())
            one[Powers((0,)*p.n_vars())] = 1.0
            self.assert_(q == one)

    def test_power_one_gives_same(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 5)
            p = rand_poly(vars)
            q = p**1
            self.assert_(q == p)

    def test_power_two_gives_mul_self(self):
        for case in xrange(self.n_cases):
            vars = random.randint(1, 5)
            p = rand_poly(vars)
            q = p**2
            self.assert_((q-p*p).l_infinity_norm() < 1.0e-14, (q, p))

    def test_power_by_doubling_vs_usual_method(self):
        tol = 1.0e-12
        for e in xrange(5):
            for case in xrange(self.n_cases):
                vars = random.randint(1, 5)
                p = rand_poly(vars)
                q0 = p**e
                q1 = p.pow_doubling(e)
                err = (q0-q1).l_infinity_norm()
                if err > tol:
                    print p
                    print q0
                    print q1
                    print err
                    print
                self.assert_(err < tol)

    def test_pow_associative_up_to_numerical_error(self):
        tol = 1.0e-12
        for e in xrange(1, 5):
            for case in xrange(self.n_cases):
                vars = random.randint(1, 5)
                p = rand_poly(vars)
                q0 = p**(e-1)
                q1 = p**e
                self.assert_((p*q0-q1).l_infinity_norm() < tol)

class Substitute(unittest.TestCase):

    def _identity(self, vars, i):
        p = Polynomial(vars)
        pows = [0, ]*vars
        pows[i] = 1
        p[Powers(tuple(pows))] = 1.0
        return p

    def _square(self, vars, i):
        p = Polynomial(vars)
        pows = [0, ]*vars
        pows[i] = 2
        p[Powers(tuple(pows))] = 1.0
        return p

    def test_sub_identity(self):
        vars = 3
        p = Polynomial(vars)
        p[Powers((2, 1, 0))] = -0.3+0.2J
        p[Powers((1, 1, 0))] = 0.995
        ids = [self._identity(vars, i) for i in xrange(vars)]
        q = p.substitute(ids)
        self.assert_(p==q)

    def test_sub_squares(self):
        vars = 3
        p = Polynomial(vars)
        p[Powers((2, 1, 0))] = -0.3+0.2J
        p[Powers((1, 1, 0))] = 0.995
        sqs = [self._square(vars, i) for i in xrange(vars)]
        q = p.substitute(sqs)
        self.assertEquals(len(q), len(p))
        for m, c in p.powers_and_coefficients():
            doubled = Powers(tuple([2*i for i in m]))
            self.assertEquals(q[doubled], c)

    def test_sub_sum(self):
        vars = 3
        p = Polynomial(vars)
        a = -0.3+0.2J
        b = 0.995
        p += a*Polynomial.Monomial((2, 1, 3))
        p += b*Polynomial.Monomial((1, 1, 0))
        ids = [self._identity(vars, i) for i in xrange(vars)]
        ids[0] = self._identity(vars, 0)+2.0*self._identity(vars, 2)
        q = p.substitute(ids)
        self.assertEquals(len(q), 5)
        self.assertEquals(q[Powers((2, 1, 3))], a)
        self.assertEquals(q[Powers((1, 1, 4))], 4.0*a)
        self.assertEquals(q[Powers((0, 1, 5))], 4.0*a)
        self.assertEquals(q[Powers((1, 1, 0))], b)
        self.assertEquals(q[Powers((0, 1, 1))], 2.0*b)

    def test_sub_neg(self):
        vars = 3
        p = Polynomial(vars)
        a = -0.3+0.2J
        b = 0.995
        p += a*Polynomial.Monomial((2, 2, 3))
        p += b*Polynomial.Monomial((1, 1, 0))
        neg_ids = [-self._identity(vars, i) for i in xrange(vars)]
        q = p.substitute(neg_ids)
        self.assertEquals(len(q), len(p))
        self.assertEquals(q[Powers((2, 2, 3))], -a)
        self.assertEquals(q[Powers((1, 1, 0))], b)
        
class NumericalTruncation(unittest.TestCase):

    def setUp(self):
        p = Polynomial(3, terms={Powers((2, 1, 3)): 0.1,
                                 Powers((1, 1, 4)): 0.2,
                                 Powers((0, 1, 5)): -0.001,
                                 Powers((1, 1, 0)): 0.0000001+0.01J,
                                 Powers((0, 1, 1)): 0.0000001-0.0003J})
        self.p = p

    def test_truncate_smaller_than_exists_returns_same(self):
        p = self.p
        q = p.with_small_coeffs_removed(0.0)
        self.assert_(q == p)

    def test_truncate_removes_proper_terms(self):
        p = self.p
        q0 = p.with_small_coeffs_removed(0.001)
        q1 = Polynomial(3, terms={Powers((2, 1, 3)): 0.1,
                                  Powers((1, 1, 4)): 0.2,
                                  Powers((1, 1, 0)): 0.0000001+0.01J})
        self.assert_(q0 == q1)

    def test_truncate_all_gives_zero(self):
        p = self.p
        q0 = p.with_small_coeffs_removed(999.0)
        q1 = Polynomial(3)
        self.assert_(q0 == q1)
        self.assert_(not q0) #redundancy
        self.assert_(not q1) #redundancy

class IsHomogeneous(unittest.TestCase):

    def test_one_is(self):
        p = Polynomial(2, terms={Powers((1, 2)): 0.5,
                                 Powers((1, 1)): 0.0})
        self.assert_(p.is_homogeneous())
        self.assert_(not p.is_homogeneous(0))
        self.assert_(not p.is_homogeneous(1))
        self.assert_(not p.is_homogeneous(2))
        self.assert_(p.is_homogeneous(3))
        self.assert_(not p.is_homogeneous(4))
        self.assert_(not p.is_homogeneous(5))
        self.assert_(not p.is_homogeneous(6))

    def test_one_is_not(self):
        p = Polynomial(2, terms={Powers((1, 2)): 0.5,
                                 Powers((1, 1)): 0.1})
        self.assert_(not p.is_homogeneous())
        self.assert_(not p.is_homogeneous(0))
        self.assert_(not p.is_homogeneous(1))
        self.assert_(not p.is_homogeneous(2))
        self.assert_(not p.is_homogeneous(3))
        self.assert_(not p.is_homogeneous(4))
        self.assert_(not p.is_homogeneous(5))
        self.assert_(not p.is_homogeneous(6))

class RealAndImag(unittest.TestCase):

    def check_polynomial(self, p, re=None, im=None):
        if re is None:
            re = p.real()
        if im is None:
            im = p.imag()
        for m, c in re.powers_and_coefficients():
            assert c != 0.0
            assert isinstance(c, float)
        for m, c in im.powers_and_coefficients():
            assert c != 0.0
            assert isinstance(c, float)
        q = re + (1.0J)*im
        self.assert_(p == q)

    def test_example(self):
        p = Polynomial(1, terms={Powers((0,)): 1+2J,
                                 Powers((1,)): 1.5,
                                 Powers((5,)): -3.6,
                                 Powers((6,)): 2.3J,
                                 Powers((7,)): -0.2J})
        self.check_polynomial(p)

    def test_real(self):
        p = Polynomial(1, terms={Powers((0,)): 1.5,
                                 Powers((5,)): -3.6})
        re = p.real()
        im = p.imag()
        self.check_polynomial(p, re, im)
        self.assert_(not im)
        
    def test_imag(self):
        p = Polynomial(1, terms={Powers((2,)): 2.3J,
                                 Powers((7,)): -0.2J})
        re = p.real()
        im = p.imag()
        self.check_polynomial(p, re, im)
        self.assert_(not re)

class CoordinateMonomial(unittest.TestCase):

    def test_one_correct(self):
        n_vars = 10
        for i in xrange(n_vars):
            m = Polynomial.CoordinateMonomial(n_vars, i)
            self.assert_(len(m) == 1)
            po, co = list(m.powers_and_coefficients())[0]
            self.assert_(co == 1.0)
            for j in xrange(n_vars):
                if i == j:
                    self.assert_(po[j] == 1)
                else:
                    self.assert_(po[j] == 0)

class InputOutput(unittest.TestCase):

    def setUp(self):
        self.f0 = ['2',
                   '3',
                   '0 0 0.35j',
                   '1 0 -0.2',
                   '1 1 0.1-0.7j']
        self.p0 = Polynomial(2, terms={Powers((0, 0)): 0.35J,
                                       Powers((1, 0)): -0.2,
                                       Powers((1, 1)): (0.1-0.7J)})
        self.f1 = ['4',
                   '2',
                   '0 0 1 0 0.35j',
                   '1 0 1 1 -0.2']
        self.p1 = Polynomial(4, terms={Powers((0, 1, 0, 0)): 0.35J,
                                       Powers((1, 1, 0, 1)): -0.2})

    def test_read_single_old_format(self):
        p0 = read_ascii_polynomial(iter(self.f0), is_xxpp_format=1)
        self.assert_(p0 == self.p0)

    def test_read_two_old_format(self):
        p1 = read_ascii_polynomial(iter(self.f1), is_xxpp_format=1)
        self.assert_(p1 == self.p1)

    def test_read_single_new_format(self):
        p0 = read_ascii_polynomial(iter(self.f0), is_xxpp_format=0)
        self.assert_(p0 == self.p0)

    def test_write_read_single_old_format(self):
        f0 = open(tmp_path+'_tmp.pol', 'w')
        write_ascii_polynomial(f0, self.p0, is_xxpp_format=1)
        f0.close()
        f0 = open(tmp_path+'_tmp.pol', 'r')
        p0 = read_ascii_polynomial(f0, is_xxpp_format=1)
        f0.close()
        self.assert_(p0 == self.p0)

    def test_write_read_two_old_format(self):
        f1 = open(tmp_path+'_tmp.pol', 'w')
        write_ascii_polynomial(f1, self.p1, is_xxpp_format=1)
        f1.close()
        f1 = open(tmp_path+'_tmp.pol', 'r')
        p1 = read_ascii_polynomial(f1, is_xxpp_format=1)
        f1.close()
        self.assert_(p1 == self.p1, (p1, self.p1))

    def test_write_read_single_new_format(self):
        f0 = open(tmp_path+'_tmp.pol', 'w')
        write_ascii_polynomial(f0, self.p0, is_xxpp_format=0)
        f0.close()
        f0 = open(tmp_path+'_tmp.pol', 'r')
        p0 = read_ascii_polynomial(f0, is_xxpp_format=0)
        f0.close()
        self.assert_(p0 == self.p0)

    def test_xxpp_to_xpxp(self):
        xxpp = [1,2,3,4,5,6]
        xpxp = xxpp_to_xpxp(xxpp)
        self.assert_(xpxp == [1, 4, 2, 5, 3, 6], xpxp)

    def test_xpxp_to_xxpp(self):
        xpxp = [1,2,3,4,5,6]
        xxpp = xpxp_to_xxpp(xpxp)
        self.assert_(xxpp == [1, 3, 5, 2, 4, 6], xxpp)

def suite():
    suites = []
    suites.append(unittest.makeSuite(Comparison))
    suites.append(unittest.makeSuite(InputOutput))
    suites.append(unittest.makeSuite(CoordinateMonomial))
    suites.append(unittest.makeSuite(Pow))
    suites.append(unittest.makeSuite(PowersAndCoeff))
    suites.append(unittest.makeSuite(RealAndImag))
    suites.append(unittest.makeSuite(IsHomogeneous))
    suites.append(unittest.makeSuite(Norms))
    suites.append(unittest.makeSuite(Degree))
    suites.append(unittest.makeSuite(NumericalTruncation))
    suites.append(unittest.makeSuite(Constant))
    suites.append(unittest.makeSuite(Substitute))
    suites.append(unittest.makeSuite(FromPowers))
    suites.append(unittest.makeSuite(Default))
    suites.append(unittest.makeSuite(Call))
    suites.append(unittest.makeSuite(Ring))
    suites.append(unittest.makeSuite(Zero))
    suites.append(unittest.makeSuite(Diff))
    suites.append(unittest.makeSuite(Poisson))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
