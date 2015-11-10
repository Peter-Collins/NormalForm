"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Polynomial

PURPOSE:

New class representing multivariate polynomials with fully tested Ring
logic.

NOTES:

This improves on the previous version:-

 - We do not store zero monomials and all monomial powers appear
   exactly once (that is, two monomials with the same powers are
   always represented by a single monomial with coefficient given by
   the sum of the coefficients of the original monomials).

 - [OBSOLETE] The class does not make any assumptions on which Powers
   representation is being used (it only uses the public interface
   from Monomial.py).

@attention: [OBSOLETE] one should be careful that the different Powers
representations recognise equality to each other in terms of the tuple
representation of their powers, otherwise the semantics of using
Powers as keys in a map from powers to coefficients will be faulty.
Unit tests are included to confirm that the logic is correct.

@attention: now, we enforce that exactly one Powers representation is
being used.  This enables the different Power representations to
implement, for example, their own hash functions (for efficiency)
without interfering with the actions of dictionary lookup within the
Polynomial.

Instead of doing:

  poly.add_monomial(powers, coeff)

do:

  poly[powers] += coeff

"""

from Powers import PowersBase, Powers
from Utility import factorial
from ShelvedDict import ShelvedDict

class Polynomial:
    """

    Multivariate polynomials in sparse/distributed representation.

    """
    _new_dict = dict #for memory-based
    #_new_dict = ShelvedDict #for disk-based

    ####################################################################
    # Constructor and copy constructor
    ####################################################################

    def __init__(self, n_vars, terms=None):
        """

        @param n_vars: the number of variables (a non-neg int).

        @param terms: (optional) terms of the polynomial in the form
        of a mapping (e.g., dict) from powers (key) to coefficient
        (value).

        """
        assert n_vars > 0
        self._n_vars = n_vars
        self._terms = Polynomial._new_dict()
        if terms:
            for po, co in terms.iteritems():
                self[po] = co #secure; cannot set zero

    def copy(self):
        """

        Copy constructor.

        """
        res = Polynomial(self._n_vars)
        for powers, coeff in self._terms.iteritems():
            assert coeff != 0.0 #security
            res[powers] = coeff
        return res

    def __del__(self):
        del self._terms

    ####################################################################
    # Constructors for Monomials, IdMonomials, One
    # These should be renamed make_monomial, etc.
    ####################################################################

    def Monomial(powers, coeff=1.0):
        """

        A factory method for making single-term polynomials.

        @param powers: a tuple of non-neg integer powers.

        @param coeff: (optional) coefficient; 1.0 is assumed.

        @return: a polynomial containing a single term (or no terms, in
        the case where coeff is zero).

        """
        if isinstance(powers, tuple):
            powers = Powers(powers)
        assert isinstance(powers, Powers)
        return Polynomial(len(powers), terms={powers: coeff})
    Monomial = staticmethod(Monomial)

    def CoordinateMonomial(n_vars, index):
        """

        Factory for making certain of the basis monomials, namely those
        which contain a single coordinate.

        @param n_vars: the number of variables.
        
        @param index: the index of the single pure power.

        @attention: these monomials alone do **not**, of course, form the
        whole basis, one would have to take *all* monomials with unit
        coefficient, i.e., including those that are products of single
        pure powers.

        """
        assert index >= 0
        assert index < n_vars
        powers = [int(i == index) for i in xrange(n_vars)]
        return Polynomial.Monomial(tuple(powers))
    CoordinateMonomial = staticmethod(CoordinateMonomial)

    def One(n_vars):
        """

        Factory for making the single-term constant polynomial that we can
        identify with the scalar 1.0.

        @param n_vars: the number of variables.

        """
        assert n_vars > 0
        return Polynomial.Monomial((0,)*n_vars)
    One = staticmethod(One)

    ####################################################################
    # String representations
    ####################################################################

    def __repr__(self):
        body = '%d' % (self._n_vars)
        if self._terms:
            body += ', terms=%s' % repr(self._terms)
        return 'Polynomial(%s)' % body

    def __str__(self):
        """

        As repr, but with nicer layout.

        """
        res = 'Polynomial(%s' % self._n_vars
        if self._terms:
            res += ', terms={ \\\n'
            lis = []
            for powers, coeff in self._terms.iteritems():
                lis.append('    %s: %s' % (repr(powers), repr(coeff)))
            res += ', \\\n'.join(lis)
            res += '}'
        res += ')'
        return res

    ####################################################################
    # Inspectors
    ####################################################################

    def n_vars(self):
        """Return the number of variables."""
        return self._n_vars

    def n_terms(self):
        """Return the number of terms in the polynomial."""
        return len(self)

    def powers_and_coefficients(self):
        """Return an iterator over the (powers, coefficient)-pairs."""
        return self._terms.iteritems()

    def coeffs_as_list(self):
        """This could be done as return list(self._terms.itervalues())."""
        return list(self._terms.itervalues())

    def __len__(self):
        """

        Return the number of (non-zero) terms in the Polynomial.
        Note: overriding __len__ (without overriding __nonzero__)
        means that we can use a polynomial in a boolean test; the test
        will produce True for non-zero polynomials (those which have
        terms) and False for those which do not (i.e., those which
        have no terms and thus represent the zero polynomial over the
        corresponding number of variables.

        """
        return len(self._terms)

    def degree(self):
        """Homogeneous degree of the polynomial; maximal degree of powers."""
        res = 0
        for powers, coeff in self._terms.iteritems():
            assert coeff != 0.0 #security
            deg = powers.degree()
            if deg > res:
                res = deg
        return res

    def l1_norm(self):
        """The $l_1$-norm; sum of absolute values of coefficients."""
        res = 0.0
        for coeff in self._terms.itervalues():
            assert coeff != 0.0 #security
            res += abs(coeff)
        return res

    def l_infinity_norm(self):
        """The $l_\\infty$-norm; maximum absolute value of coefficient."""
        res = 0.0
        for coeff in self._terms.itervalues():
            assert coeff != 0.0 #security
            if abs(coeff) > res:
                res = abs(coeff)
        return res

    ####################################################################
    # Comparators
    ####################################################################

    def __cmp__(self, other):
        """

        We prefer direct (rich) equality checks for now.

        """
        raise NotImplementedError

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            if self._n_vars == other._n_vars:
                return self._terms == other._terms
            else:
                return False
        else:
            raise TypeError, 'Must compare Polynomial with similar.'

    ####################################################################
    # Item access
    ####################################################################

    def has_term(self, term):
        """

        Take care with this; one must match the monomial type.

        """
        assert isinstance(term, Powers)
        if len(term)==self._n_vars:
            return self._terms.has_key(term)
        else:
            raise IndexError, "Keys must have correct number of variables"

    def __getitem__(self, key):
        """

        This overrides getitem to return a default zero.

        """
        assert isinstance(key, Powers)
        if len(key) == self._n_vars:
            return self._terms.get(key, 0.0) #self._terms[key]
        else:
            raise IndexError, "Keys must have correct number of variables"

    def __setitem__(self, key, item):
        """

        This overrides setitem so that we can never make
        zero-coefficient monomials, even via self[m] += c.

        """
        assert isinstance(key, Powers), 'Polynomials are Powers -> coeff'
        if len(key) == self._n_vars:
            if item == 0.0:
                if self._terms.has_key(key):
                    del self._terms[key]
                else:
                    pass #IMPORTANT! WE NEVER STORE ZERO ENTRIES.
            else:
                self._terms[key] = item
        else:
            raise IndexError, "Keys must have correct number of variables"

    def __delitem__(self, key):
        assert isinstance(key, Powers)
        if len(key) == self._n_vars:
            self._terms.__delitem__(key)
        else:
            raise IndexError, "Keys must have correct number of variables"

    ####################################################################
    # Property tests
    ####################################################################

    def is_constant(self):
        if len(self._terms) > 1:
            return False
        if len(self._terms):
            monomial = self._terms.keys()[0]
            for i in xrange(len(monomial)):
                if monomial[i] != 0:
                    return False
        return True

    def is_homogeneous(self, degree=-1):
        """

        Check whether the polynomial is homogeneous (of degree).

        """
        if not self._terms:
            return True
        deg = degree
        for m, c in self._terms.iteritems():
            assert c != 0.0
            if deg == -1:
                deg = m.degree()
            else:
                if deg != m.degree():
                    return False
        return True

    ####################################################################
    # Functionoid
    ####################################################################

    def __call__(self, args):
        """

        Evaluate the polynomial with numerical values.  In order to
        evaluate at a vector of other polynomials, use the substitute
        method instead.

        """
        if len(args) == self._n_vars:
            result = 0.0
            for m, c in self._terms.iteritems():
                assert c != 0.0 #security
                result += c * m(args)
            return result
        else:
            raise IndexError, "Wrong number of arguments"

    ####################################################################
    # Arithmetic operations
    ####################################################################

    def __neg__(self):
        res = Polynomial(self._n_vars)
        for m, c, in self._terms.iteritems():
            assert c != 0.0 #security
            res[m] = -c
        return res

    def __mul__(self, other):
        assert isinstance(other, Polynomial)
        res = Polynomial(self._n_vars)
        for ms, cs in self._terms.iteritems():
            assert cs != 0.0 #security
            for mo, co in other._terms.iteritems():
                assert co != 0.0 #security
                t = ms * mo
                if res.has_term(t):
                    res[t] += cs * co
                else:
                    res[t] = cs * co
        return res

    def __rmul__(self, other):
        """

        Used to multiply by float, int, monomial, etc.

        """
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            res[m] = other * c
        return res

    def __pow__(self, other):
        """

        This should be made more efficient by the doubling method.  I
        should do that later.  This should also be made generic.

        """
        assert isinstance(other, int)
        assert other >= 0
        if other == 0:
            return Polynomial.One(self._n_vars) #NOT VERY CLEAN, BUT IT WORKS.
        if not self._terms:
            return Polynomial(self._n_vars) #zero
        res = 1.0 * self
        for i in range(other - 1):
            res = res*self
        return res

    def pow_doubling(self, e):
        """

        The usual" left-to-right efficient exponentiation algorithm.
        e must be vanilla integer, not BigDec.

        """
        assert isinstance(e, int)
        assert e >= 0
        mask = 1L
        while mask <= e:
            mask <<= 1
        mask >>= 1  # position of e's most-significant bit
        result = Polynomial.One(self._n_vars)
        while mask:
            result *= result
            if e & mask:
                result *= self
            mask >>= 1
        return result

    def __iadd__(self, other):
        """

        Implements in-place addition of another polynomial.
        For safety, we handle the case of other is self;
        we ensure that we do not iter over self.

        """
        assert isinstance(other, Polynomial), other.__class__
        if other is self:
            for m in self._terms.keys():
                self._terms[m] *= 2.0
            return self
        if self._n_vars == other._n_vars:
            for m, c in other._terms.iteritems():
                assert c != 0.0 #security
                self[m] += c
            return self
        else:
            raise IndexError, "Number of variables must match"

    def __add__(self, other):
        """

        Addition; not in-place.
        This could perhaps be made more efficient by detecting the
        case when other is self and making a new, doubled, copy of self.

        """
        assert isinstance(other, Polynomial)
        if self._n_vars == other._n_vars:
            if len(self) >= len(other):
                a = self
                b = other
            else:
                a = other
                b = self
            res = a.copy()
            res += b
            return res
        else:
            raise IndexError, "Number of variables must match"

    def __isub__(self, other):
        """

        Implements in-place subtraction of another polynomial.
        For safety, we handle the case of other is self;
        in that case, we delete all terms of the current polynomial.

        """
        assert isinstance(other, Polynomial)
        if other is self:
            for m in self._terms.keys():
                del self._terms[m]
            return self
        if self._n_vars == other._n_vars:
            for m, c in other._terms.iteritems():
                assert c != 0.0 #security
                self[m] -= c
            return self
        else:
            raise IndexError, "Number of variables must match"

    def __sub__(self, other):
        """

        Subtraction; not in-place.
        This could perhaps be made more efficient by detecting the
        case when other is self and returning a new zero polynomial.

        """
        res = self.copy()
        res -= other
        return res

    ####################################################################
    # Other operations
    ####################################################################

    def substitute(self, args):
        """

        Substitute a vector of polynomials for the coordinates in this
        polynomial; this is basically the same as __call__ but uses a
        different zero, namely a zero polynomial.

        Note also that one
        must in general use the one polynomial to take care of the
        case where powers evaluate to a constant.

        """
        if len(args) == self._n_vars:
            assert isinstance(args[0], Polynomial)
            result_space = args[0].n_vars()
            for a in args[1:]:
                assert isinstance(a, Polynomial)
                if a.n_vars() != result_space:
                    raise IndexError, 'incompatible substitution'
            result = Polynomial(result_space)
            for po, co in self._terms.iteritems():
                res = Polynomial.One(result_space) #ensure poly
                assert co != 0.0 #security
                for i in range(len(args)):
                    if po[i] > 0:
                        res *= args[i] ** po[i]
                result += co * res
            assert isinstance(result, Polynomial)
            return result
        else:
            raise IndexError, "Wrong number of arguments"

    def diff(self, var):
        """

        Partially differentiate self with respect to the variable indexed.

        """
        assert isinstance(var, int)
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            dc, dm = m.diff(var)
            if dc == 0.0:
                pass #IMPORTANT! WE DO NOT STORE ZERO ENTRIES.
            else:
                if res.has_term(dm):
                    res[dm] += c * dc
                else:
                    res[dm] = c * dc
        return res

    def diff_pow(self, var, pow):
        """

        Partially differentiate self with respect to the variable indexed.

        """
        assert isinstance(var, int)
        assert isinstance(pow, int)
        if var<0 or var>=self._n_vars:
            raise IndexError(var)
        if pow == 0:
            return self
        if pow == 1:
            return self.diff(var)
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            dc, dm = m.diff_pow(var, pow)
            if dc == 0.0:
                pass #IMPORTANT! WE DO NOT STORE ZERO ENTRIES.
            else:
                if res.has_term(dm):
                    res[dm] += c * dc
                else:
                    res[dm] = c * dc
        return res

    def poisson_bracket(self, other):
        """

        Assumes that the variables are ordered x1, p1, x2, p2, ...,
        xn, pn The two polynomials must have the same number of
        variables and the number is assumed to be even.  This really
        belongs in the LieAlgebra module; I should shift it to
        there.  It sits here pending an update.

        """
        assert isinstance(other, Polynomial)
        if (self._n_vars)%2:
            raise IndexError, "Poisson bracket not defined for odd num coords"
        if self._n_vars == other._n_vars:
            res = Polynomial(self._n_vars)
            for i in range(0, self._n_vars, 2):
                p = (self.diff(i)*other.diff(i + 1))
                q = (self.diff(i + 1)*other.diff(i))
                res += (p-q)
            return res
        else:
            raise IndexError, "Number of variables must match"

    ####################################################################
    # Partitions
    ####################################################################

    def homogeneous(self, degree, up_to=None):
        """

        returns the homogeneous sub polynomial of given degree.

        """
        if up_to==None:
            up_to = degree+1
        assert degree < up_to
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            if degree <= m.degree() < up_to:
                res[m] = c
        return res

    def partition(self, pred):
        """

        Partitions the polynomial into two polys returns polys t, f
        for true and false under the predicate pred.

        """
        t = Polynomial(self._n_vars)
        f = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            if pred(m):
                t[m] = c
            else:
                f[m] = c
        return t, f

    def real(self):
        """

        Return the real polynomial forming the real part of a poly.

        """
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0
            if isinstance(c, float):
                res[m] = c
            else:
                res[m] = c.real
        return res

    def imag(self):
        """

        Return the real polynomial forming the imag part of a poly.

        """
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0
            if isinstance(c, complex):
                res[m] = c.imag
        return res

    def with_small_coeffs_removed(self, tolerance):
        """

        A copy of the polynomial, containing only terms whose
        coefficients are larger in size than the specified numerical
        tolerance.

        """
        res = Polynomial(self._n_vars)
        for m, c in self._terms.iteritems():
            assert c != 0.0 #security
            if abs(c) > tolerance:
                res[m] = c
        return res

