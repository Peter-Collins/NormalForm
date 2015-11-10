"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: IntegralExtractor

PURPOSE:

Extract the complex integrals of motion from a complex normal form.

NOTES:

This will also handle the case of partial normalizations, collecting
the relevant simple integrals (those that are a straightforward
product of a simple complex normal form coordinate with its respective
conjugate momentum) and expressing the remainder of the normal form
Hamiltonian as a single non-simple integral.

"""

from sets import Set
from PolynomialRing import PolynomialRing
from Polynomial import Polynomial

class IntegralExtractor:

    """

    Extract integrals from a complex normalized Hamiltonian, including
    the partially-normalized case.

    We define a simple integral to be an integral that may be written
    as the product of a single normal form coordinate and its
    associated conjugate momentum.

    The canonically-conjugate variable pairs corresponding to simple
    integrals are automatically indentified, and hence those
    corresponding to lost integrals are also found.

    The class also finds the remainder part of the Hamiltonian, i.e.,
    the part which cannot be expressed as a polynomial of the simple
    integrals.  This non-simple part of the Hamiltonian is itself an
    integral of the motion, which we call the non-simple integral.

    The simple part of the Hamiltonian is then expressed in a compact
    form as a polynomial over the simple integrals.

    Expressions for all the integrals in terms of the normal form
    coordinates are also computed.  We always arrange that the
    non-simple integral is the last one in the list.

    """

    def __init__(self, lie_algebra):
        self._alg = lie_algebra

    def set_complex_normal_hamiltonian(self, h_comp):
        self._alg.check_elt(h_comp)
        self._h = h_comp

    def find_lost_simple_integrals_and_non_simple_integral(self):

        """Determine which simple integrals (i.e., those that are
        simply products of complex normal form coordinates and their
        associated conjugate momentum), are lost due to partial
        normalization terms in a complex normal form Hamiltonian.  At
        the same time, we can determine the non-simple integral which
        results from such terms."""

        #bind globals
        _alg = self._alg
        _h = self._h
        _xrange = xrange

        lost_integrals = Set()
        simple_part = _alg.zero()
        non_simple_integral = _alg.zero()
        for m, c in _h.powers_and_coefficients():
            powers = m.to_tuple()
            evens, odds = powers[0::2], powers[1::2]
            if evens == odds:
                simple_part[m] = c
            else:
                non_simple_integral[m] = c
                for i in _xrange(_alg.dof()):
                    if powers[2*i+1] != powers[2*i]:                    
                        lost_integrals.add(i)

        self._lost_integrals = lost_integrals
        self._simple_part = simple_part
        self._non_simple_integral = non_simple_integral

    def list_the_simple_integrals(self):

        """Given which simple integrals are lost, figure out which
        ones remain and make a sorted list of them to use later to
        compress powers into a more compact representation which only
        lists the simple integral coordinates."""

        _dof = self._alg.dof()
        _lost = self._lost_integrals
        
        kept = list(Set(xrange(_dof)).difference(_lost))
        kept.sort()
        assert len(kept) <= _dof
        if self._non_simple_integral:
            assert len(kept) == _dof

        self._kept_integrals = kept

    def express_simple_part_as_polynomial_over_all_integrals(self):

        """Given which simple integrals are lost, we can now express
        the simple part as a function of the remaining simple
        integrals, and the complicated integral which has already been
        calculated.  The power corresponding to the non-simple
        integral will be zero here."""

        #bind globals
        _dof = self._alg.dof()
        _monomial = Polynomial.Monomial
        _simple = self._simple_part
        _lost = self._lost_integrals
        _kept = self._kept_integrals

        if self._non_simple_integral:
            n_integrals = len(_kept)+1
        else:
            n_integrals = len(_kept)
        self._integral_ring = PolynomialRing(n_integrals)
        _ring = self._integral_ring
        
        simple_part_in_integrals = _ring.zero()
        lost_integrals = Set()
        for m, c in _simple.powers_and_coefficients():
            powers = m.to_tuple()
            evens = powers[0::2]
            powers_without_lost = tuple([evens[src] for src in _kept])
            if self._non_simple_integral:
                #we add a trailing zero for the non-simple integral
                powers_without_lost = powers_without_lost+(0,)
            simple_part_in_integrals += _ring.monomial(powers_without_lost, c)
            if 1:
                assert evens == powers[1::2]
                for i, p in enumerate(powers):
                    if p != 0:
                        assert not (i in _lost)
                if self._non_simple_integral:
                    assert len(powers_without_lost) == len(_kept)+1
                else:
                    assert len(powers_without_lost) == len(_kept) == _dof

        self._simple_in_integrals = simple_part_in_integrals

    def express_both_parts_as_polynomial_over_all_integrals(self):

        """

        Express the whole polynomial (e.g., normalized Hamiltonian) as
        a polynomial over all the integrals.

        """

        simple_in_integrals = self._simple_in_integrals
        non_simple_integral = self._non_simple_integral

        total_in_integrals = simple_in_integrals.copy()
        if non_simple_integral:
            pows = (0,)*len(self._kept_integrals)+(1,)
            m = self._integral_ring.monomial(pows)
            total_in_integrals += m

        self._total_in_integrals = total_in_integrals

    def express_all_integrals_in_normal_form_coords(self):

        """

        We express all the integrals in terms of the normal form
        coordinates (complex, canonically-conjugate pairs), ensuring
        that the non-simple integral (if any) always appears last in
        the list.

        MANY OF THESE ROUTINES SHOULD PROBABLY BE GENERATORS, OR
        SHOULD BE ALLOWED TO USE GIVEN LIST-LIKE OBJECTS IN WHICH TO
        DO THEIR STORAGE.  OTHERWISE, MEMORY WILL BECOME A PROBLEM FOR
        LARGE SYSTEMS.

        """

        _alg = self._alg
        _monomial = _alg.monomial
        _kept = self._kept_integrals

        ints_in_coords = []
        for i in _kept:
            powers = [0,]*(_alg.n_vars())
            powers[2*i] = 1
            powers[2*i+1] = 1
            integral = _alg.zero()
            integral += _monomial(tuple(powers))
            ints_in_coords.append(integral)

        if self._non_simple_integral:
            ints_in_coords.append(self._non_simple_integral)
            
        self._integrals_in_coords = ints_in_coords
