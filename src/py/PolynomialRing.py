"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: PolynomialRing

PURPOSE:

An interface and straightforward implementation for polynomial rings.

NOTES:

In order to provide a more flexible (and re-usable) design, the
Polynomial class has no knowledge of polynomial rings at all.  This
makes Polynomial into a useful tool without the overhead of always
having to manage Polynomials through PolynomialRings.

Our main use for PolynomialRing in the normal form code is to derive
LieAlgebras, in particular the ClassicalLieAlgebra and
SemiClassicalLieAlgebra, which are both types of polynomial ring.

Also provided here is a pretty-printer for complex numbers, which
produces a fixed-width representation, making it easier to read tables
of complex coefficients, for example.

Reading and writing of polynomials is done via very basic
s-expression-like syntax.  Either a single bracketed expression must
fit on one line, or only the head should appear on that line, followed
by the elements (following this same rule), followed by the
terminating bracket on its own line (followed optionally directly by a
comment matching the head).  For example (leading spaces can be
omitted):-

(polynomial
(num-variables 3)
(num-monomials 2)
(powers-format \"dense\")
(monomials
((0 0 0) 1.2+0.5j)
((1 1 0) 0.546)
);;monomials
);;polynomial

"""

from Polynomial import Polynomial
from Powers import Powers
from math import sqrt
from GradedInterface import GradedInterface
from Complex import pretty_complex
import SExprIO

class PolynomialRingInterface(GradedInterface):
    """

    Polynomial rings are graded objects.

    The default (usual) concept of grade is just the degree of the
    polynomial, and thus being isograde means being homogeneous.

    """
    def check_elt(self, elt):
        assert isinstance(elt, Polynomial)
        assert elt.n_vars() == self.n_vars()
    #graded interface:
    def grade(self, elt):
        self.check_elt(elt)
        return elt.degree()
    def is_isograde(self, elt, grade=-1):
        self.check_elt(elt)
        return elt.is_homogeneous(degree=grade)
    def isograde(self, elt, grade, up_to=None):
        self.check_elt(elt)
        return elt.homogeneous(grade, up_to=up_to)
    #finite dimensional
    def n_vars(self):
        raise NotImplementedError
    #ring
    def one(self):
        return Polynomial.One(self.n_vars())
    def zero(self):
        return Polynomial(self.n_vars())
    #default basis
    def monomial(self, powers, coeff=1.0):
        assert len(powers) == self.n_vars()
        return Polynomial.Monomial(powers, coeff)
    def coordinate_monomial(self, i, pow=1):
        pows = [0,]*self.n_vars()
        pows[i] = pow
        return Polynomial.Monomial(tuple(pows))
    def polynomial(self, terms=None):
        return Polynomial(self.n_vars(), terms)
    #i/o
    def display(self, pol_p, names=None):
        """

        Pretty-print

        """
        self.check_elt(pol_p)
        if not names:
            names = ['x_{%d}'%i for i in xrange(0, self.n_vars())]
        assert len(names) == self.n_vars()
        terms = []
        desc = '[%d terms]:\n= '%len(pol_p)
        for po, co in pol_p.powers_and_coefficients():
            term = []
            term.append(pretty_complex(co))
            for i in xrange(self.n_vars()):
                name, power = names[i], po[i]
                if power:
                    if power == 1:
                        term.append('%s'%name)
                    else:
                        term.append('%s^{%d}'%(name, power))
            terms.append(' '.join(term))
        return desc + '\n+ '.join(terms)
    #miscellaneous
    def grad(self, poly):
        self.check_elt(poly)
        return tuple([poly.diff(i) for i in xrange(0, self.n_vars())])

class PolynomialRing(PolynomialRingInterface):
    """

    Polynomial ring with explicit dimension.

    """
    def __init__(self, n_vars):
        self._n_vars = n_vars
    #polynomial ring interface:
    def n_vars(self):
        return self._n_vars

