"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: PolynomialRingIO

PURPOSE:

Provide output routines for polynomials and polynomial rings in bare
ASCII and S-expression formats.

NOTES:

I would rather keep I/O, other than extremely simple I/O and
pretty-printing, out of the main classes themselves.  Otherwise, the
classes and their interfaces start to become more complicated; their
core functionality becomes less clear.

Hence this IO class, which does little more than hold a reference to a
polynomial ring and provide the functionality.

"""

import SExprIO
from Powers import Powers
from Polynomial import Polynomial

class PolynomialRingIO:
    def __init__(self, polynomial_ring):
        self.ring = polynomial_ring
    def write_sexp_vector_of_polynomials(self, ostream, polys):
        ostream.write('(vector-of-polynomials\n')
        ostream.write('(num-polynomials %d)\n' % len(polys))
        ostream.write('(polynomials\n')
        for poly in polys:
            self.write_sexp_polynomial(ostream, poly)
        ostream.write(')\n')
        ostream.write(')\n')
    def read_sexp_vector_of_polynomials(self, istream):
        assert not SExprIO.strip_begin(istream.next(), 'vector-of-polynomials')
        #
        elts = SExprIO.elts(SExprIO.strip_braces(istream.next(), 'num-polynomials'))
        assert len(elts) == 1
        n_polys = int(elts[0])
        assert n_polys >= 1
        #
        assert not SExprIO.strip_begin(istream.next(), 'polynomials')
        #
        polys = []
        for i in xrange(n_polys):
            polys.append(self.read_sexp_polynomial(istream))
        assert not SExprIO.strip_end(istream.next(), 'polynomials')
        assert not SExprIO.strip_end(istream.next(), 'vector-of-polynomials')
        return polys
    def write_sexp_polynomial(self, ostream, poly):
        """

        Write a polynomial to an output stream in SEXP format.

        @param ostream: an iterable of input lines (e.g. file-like),

        @param poly: a polynomial,

        """
        assert poly.n_vars() == self.ring.n_vars(), 'Wrong number of variables in ring.write_sexp_polynomial.'
        ostream.write('(polynomial\n')
        ostream.write('(num-variables %d)\n'%poly.n_vars())
        ostream.write('(num-monomials %d)\n'%len(poly))
        ostream.write('(powers-format \"dense\")\n')
        ostream.write('(monomials\n')
        self._write_dense_sexp_polynomial(ostream, poly)
        ostream.write(')\n') #;;monomials\n')
        ostream.write(')\n') #;;polynomial\n')
    def read_sexp_polynomial(self, istream):
        """

        Read a polynomial from an input stream in SEXP format.

        @param istream: an iterable of input lines (e.g. file-like),

        @return: polynomial.

        """
        assert not SExprIO.strip_begin(istream.next(), 'polynomial')
        #
        elts = SExprIO.elts(SExprIO.strip_braces(istream.next(), 'num-variables'))
        assert len(elts) == 1
        n_vars = int(elts[0])
        assert n_vars >= 0
        assert n_vars == self.ring.n_vars(), 'Wrong number of variables in ring.read_sexp_polynomial.'
        #
        elts = SExprIO.elts(SExprIO.strip_braces(istream.next(), 'num-monomials'))
        assert len(elts) == 1
        n_monomials = int(elts[0])
        assert n_monomials >= 0
        #
        elts = SExprIO.elts(SExprIO.strip_braces(istream.next(), 'powers-format'))
        assert len(elts) == 1
        #
        assert not SExprIO.strip_begin(istream.next(), 'monomials')
        #
        if elts[0] == '\"dense\"':
            p = self._read_dense_sexp_polynomial(istream, n_vars, n_monomials)
        else:
            assert elts[0] == '\"sparse\"'
            p = self._read_sparse_sexp_polynomial(istream, n_vars, n_monomials)
        assert len(p) == n_monomials, 'Wrong number of monomials read.'
        #
        assert not SExprIO.strip_end(istream.next(), 'monomials')
        assert not SExprIO.strip_end(istream.next(), 'polynomial')
        return p
    def _write_dense_sexp_polynomial(self, ostream, poly):
        """

        Write a polynomial in the dense powers format, where each
        monomial is listed as n_vars integers representing
        non-negative powers, followed by a coefficient.

        """
        count = len(poly)
        for pows, coeff in poly.powers_and_coefficients():
            powers_list = pows
            ostream.write('((')
            for i, p in enumerate(powers_list):
                if i == 0:
                    ostream.write('%d'%p)
                else:
                    ostream.write(' %d'%p)
            ostream.write(') ')
            if isinstance(coeff, complex):
                ostream.write('(%+1.15e %+1.15e)'%(coeff.real, coeff.imag))
            else:
                ostream.write('(%+1.15e %+1.15e)'%(coeff, 0.0))
            count -= 1
            ostream.write(')\n')
        assert count == 0, 'Wrong number of monomials written in ring.write_sexp_polynomial.'
    def _read_dense_sexp_polynomial(self, istream, n_vars, n_monomials):
        """

        Read a polynomial in the dense powers format, where each
        monomial is listed as n_vars integers representing
        non-negative powers, followed by a coefficient.

        """
        p = self.ring.zero()
        for i in xrange(n_monomials):
            line = istream.next()
            sinner = SExprIO.strip_braces(line)
            elts = SExprIO.elts(sinner)
            spowers = SExprIO.strip_braces(' '.join(elts[:n_vars]))
            scoeff = SExprIO.strip_braces(''.join(elts[n_vars:]))+'j'
            epowers = SExprIO.elts(spowers)
            assert len(epowers) == n_vars
            powers = tuple([int(e) for e in epowers])
            coeff = complex(scoeff)
            m = coeff*Polynomial.Monomial(powers)
            assert not p.has_term(Powers(powers))
            p += m
        return p
    def _write_sparse_sexp_polynomial(self, istream, n_vars, n_monomials):
        raise NotImplementedError
    def _read_sparse_sexp_polynomial(self, istream, n_vars, n_monomials):
        raise NotImplementedError
