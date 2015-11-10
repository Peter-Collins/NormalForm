"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: IsogradeInnerTaylorCoeffs

PURPOSE:

Implement the relationship between a list of "inner" Taylor series
coefficients and the corresponding polynomial, in which these
coefficients have leading factorial multipliers.

NOTES:

"""

from Utility import factorial
from PolynomialRing import PolynomialRingInterface

class IsogradeInnerTaylorCoeffs:
    """

    Clumsily-named class that represents an identification between the
    isograde parts of a polynomial and the list of *inner* Taylor
    coefficients, i.e., those **not including the factorial
    denominator**, that one uses in the Lie Triangle.

    """

    def __init__(self, graded_algebra, offset):
        """

        @param n_vars: the number of variables for the polynomial ring.

        @param offset: an offset which specifies the grade
        associated with the zeroth element of a list.

        @attention: in the normalization Lie Triangle, one uses an
        offset of 2, i.e., the list begins with quadratic terms, where
        as in the computation of the coordinate-change maps, one uses
        an offset of 1, so that the relevant lists begin with linear
        terms.  It is the purpose of this class to put the logic of
        converting from one form to another in a single place, so that
        erroneous indices do not creep in.

        """
        assert isinstance(graded_algebra, PolynomialRingInterface)
        self.alg = graded_algebra
        self.offset = offset
        assert self.alg.n_vars() >= 0
        assert self.offset >= 0

    def lie_algebra(self):
        return self.alg
    
    def poly_to_inner_taylor(self, poly, i):
        """

        Convert a polynomial to an inner Taylor coefficient by
        taking the grade (i+offset) part and multiplying by
        factorial(i).

        """
        self.alg.check_elt(poly)
        return (factorial(i)*self.alg.isograde(poly, (i+self.offset)))

    def inner_taylor_to_poly(self, poly, i):
        """

        Convert an inner Taylor coefficient of grade (i+offset) to
        a polynomial by dividing by factorial(i).

        """
        self.alg.check_elt(poly)
        gra = self.alg.grade(poly)
        assert (gra == 0) or (gra == i+self.offset)
        return (1.0/factorial(i))*poly

    def list_to_poly(self, p_list, poly=None):
        """

        Dof, list of polys, offset so that input[i] has grade
        i+offset.

        """
        if poly == None:
            poly = self.alg.zero()
        else:
            assert poly == self.alg.zero()
        for i, poly_i in enumerate(p_list):
            poly += self.inner_taylor_to_poly(poly_i, i)
        return poly

    def poly_to_list(self, poly, p_list=None):
        """

        Dof, poly, offset so that output list[i] has grade
        i+offset.  Mutates the input list-like.

        """
        if p_list == None:
            p_list = []
        else:
            assert p_list == []
        for i in xrange(0, self.alg.grade(poly)+1-self.offset):
            p_list.append(self.poly_to_inner_taylor(poly, i))
        return p_list

