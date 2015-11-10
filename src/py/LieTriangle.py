"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: LieTriangleGenerator

PURPOSE:

Implements the Lie Triangle for normalization of a complex diagonal
Hamiltonian with unknown generating function.

NOTES:

This version is implemented more like a [py-]generator, in that we can
ask for the next highest graded terms of both the normalised
Hamiltonian and the generator by providing the corresponding terms of
the original Hamiltonian.

The Deprit triangle is one of the methods for computing the
normalisation and coordinate changes.  There are others.  These
various methods are referred to by Murdock as Normal Form Formats.
Each format is either (a) iterative, or (b) recursive.  Each format is
either (i) generative (working via a generating function), or (ii)
direct.  The Deprit triangle, presented here, is of type (b)(i);
recursive and generative.

TODO/FUTURE:

(1) I suggest making the lie triangle accept an iterator-like object
    that provides a next() method for the input Hamiltonian (giving
    the terms of the correct grade --- note that we must still handle
    the factorial pre-factor within the triangle), and would itself
    provide a next() method so that the computation can be driven
    on-demand.

(2) Only the previous and current rows of the triangle (i.e., the
    intermediate quantities $H_{i}^{j}$ for (previous <= (i + j) <=
    current) are needed at each stage.  We can therefore safely delete
    (or release --- resource-wise) the earlier rows.

(3) The generating function list (i.e., the removed factorial
    pre-factor list of $W_i$'s) must be kept in order to perform the
    coordinate changes, etc.  However, for the Deprit triangle itself,
    we only require the previous generating function term [check
    this], so we could release earlier ones to a disk-based resource.

"""

import logging
logger = logging.getLogger()

from IsogradeInnerTaylorCoeffs import IsogradeInnerTaylorCoeffs
from Utility import binomial, factorial
from Powers import Powers

def _make_db_key(i, j):
    """

    Construct a key to store one of the intermediate Deprit triangle
    terms in a dictionary.  This function was originally intended for
    use with shelve files implementing persistent lists.

    """
    return '(%d,%d)' % (i, j)

class LieTriangle:
    """

    Computes and organizes normal form calculations Note: Generating
    functions w number from 0.

    """

    def __init__(self,
		 lie_alg,
		 quadratic_part,
                 state = None):
        """

        @param lie_alg: is a LieAlgebra (Classical or Semiclassical).

        @param quadratic_part: the quadratic part of the Hamiltonian
        to be normalised.  This serves to define the normal form
        style, i.e., the Lie operator, the partition method for the
        known terms in each triangle row, and the method for solution
        of the corresponding homological equation.

        A valid state consists of a 3-tuple of elements as follows:-

        @param h_ij_dict: a dictionary-like (shelves are acceptable
        also, as we don't use write-back semantics) object used to
        store the intermediate Deprit triangle quantities.

        @param w_i_list: a list-like object, used to store the inner
        (non-factorial denominator) terms of the generating function.

        @param current_row: the index, i, of the current row in the
        triangle.  This must be compatible with h_ij_dict and
        w_i_list.  For example, if we are beginning a new computation
        from scratch, then h_ij_dict is an empty dictionary-like
        object, w_i_list is an empty list-like object, and current_row
        is -1, indicating that we have not yet computed row 0 of the
        triangle.
        
        @attention: In order to get a polynomial for w from w_i_list,
        one has to sum the parts of w_list, each divided by the
        corresponding factorial, i.e., $w = \sum_i
        \mathrm{wlist}[i]/(i!)$.  We can do this via the helper class
        IsogradeInnerTaylorCoeffs, so that the logic is encapsulated
        in a single place.

	Correct input of current_row (the last one calculated) and h
	and w lists allows you to get the next row, thereby continuing
	an earlier calculation.

        TODO/FUTURE:

	We should probably encapsulate the dict, list, current_row
	into a single state object that can be persistent on disk, or
	loaded/saved.  If so, UUIDs would probably prove useful in
	file names?

        """
        #sanity check
        lie_alg.check_elt(quadratic_part)

	#store members
        self._lie_alg = lie_alg
        self._h_dc_2 = quadratic_part

        #handle default values for a new computation
        if state is None:
            state = ({}, [], -1)
        self._h_ij, self._w, self._current_row = state

	#check input list and dict against current_row value
        assert len(self._w) == (self._current_row + 1)
        for i in xrange(0, self._current_row + 1):
            for j in xrange(0, i + 1):
                assert self._h_ij.has_key(_make_db_key(j, i-j))

	#compute polynomial factorial conversions
        self._iso = IsogradeInnerTaylorCoeffs(self.lie_algebra(),
                                              offset=2)

        #process the quadratic part, preparing the Lie operator
        self.check_quadratic_part()
        self._freq = []
        self.extract_fundamental_frequencies()

    def get_state(self):
        return (self._h_ij, self._w, self._current_row)

    def lie_algebra(self):
        """

        Return the current Lie algebra.

        """
        return self._lie_alg

    def get_isograde_list_handler(self):
        """

        Return an iso-grade list handler for converting between
        polynomials and non-factorial-prefactor list representations.

        """
        return self._iso

    def filter(self, powers):
        """
    
        Filter must look only at the first 2*dof variables.  (This
        comment is a reminder that the Lie algebra might be a
        semiclassical one.)
    
        We filter out those monomials that may be factored into
        products of configuration with momentum variables.
    
        What to do with h-bar-only terms in the case of a
        semi-classical Lie algebra?  We include those too.
    
        """
        for i in xrange(self._lie_alg.dof()):
            if(powers[2*i+1] !=  powers[2*i]):
                #term in remainder
                return False
        #otherwise, term in normal form
        return True

    def check_quadratic_part(self):
        """

        We ensure that the complex input Hamiltonian has diagonal
        quadratic part.

        """
        alg = self.lie_algebra()
        
        #ensure basics
        alg.check_elt(self._h_dc_2)
        logger.info('Diagonal complex H has correct #variables.')
        
        #ensure no constant or linear terms
        assert alg.is_isograde(self._h_dc_2, 2)
        logger.info('Diagonal complex H has zero constant part.')
        logger.info('Diagonal complex H has zero linear part.')
        
        #ensure diagonal quadratic part
        assert alg.is_diagonal_polynomial(self._h_dc_2)
        logger.info('Diagonal complex H has diagonal quadratic part.')

    def extract_fundamental_frequencies(self):
        """

        We directly examine the quadratic part in order to determine
        the fundamental frequencies.

        It is much less error-prone to perform this examination than
        it is to infer the frequencies from the eigenvalues; if one
        were to change the form of the complexification, it would
        change these frequencies.  Extracting them here avoids any
        problems that such changes might cause and ensures that we
        deal with the correct Lie operator.

        """
        alg = self.lie_algebra()
        assert len(self._freq) == 0
        for i in xrange(0, alg.dof()):
            q_ind = 2*i
            p_ind = q_ind+1
            #classical and semi-classical
            powers = [0]*(alg.n_vars())
            powers[q_ind] = 1
            powers[p_ind] = 1
            self._freq.append(self._h_dc_2[Powers(tuple(powers))])
        assert len(self._freq) == alg.dof()
        logger.info('Diagonal complex H has fundamental frequencies:')
        for freq in self._freq:
            logger.info(freq)
    
    def compute_normal_form_and_generating_function(self, h_term):
        """

        Given the next iso-grade part of the original Hamiltonian,
        compute the next iso-grade parts of the normalised Hamiltonian
        and the generator.

        @param h_term: if the next row to be computed is i, then this
        shuld be the (i + 2)-iso-grade part of the original
        Hamiltonian, extracted directly from teh Hamiltonian, without
        any manipulation of factorial pre-factors.

        NOTES:

        Note that this method returns a pair of polynomials that may
        be added directly to the generator and the normalised
        Hamiltonian.

        These quantities are _NOT_ the $w_i$ and $k_i$; they are
        factorially-weighted versions of them.

        Internally, we store a list of the actual $w_i$.  If it is
        needed to make a polynomial from these, one must sum them
        being careful to include the proper factorial denominator,
        i.e., $K = \sum \frac{1}{n!}K_i$.

        """
        #prepare triangle row
        self._current_row += 1
        i = self._current_row
        assert (i >= 0)
        assert (len(self._w) == i)

        #prepare hamiltonian term, removing factorial division to give h_i
        assert self._lie_alg.is_isograde(h_term, i + 2)
        pre_factor = float(factorial(i))
        h_i = pre_factor * h_term
        self._h_ij[_make_db_key(i, 0)] = h_i
	
        #partition the known terms into normalised and remainder
        logger.info('Computing triangle row minus unknown term')
        known_terms = self.triangle_minus_unknown_term(0, i)
        logger.info('Partitioning the known terms')
        normalised, remainder = known_terms.partition(self.filter)
        k_i = normalised
        self._correct_row_using_partition(normalised, remainder)

        #solve homological equation
        w_i = self.solve_homological_eqn(remainder)
        self._w.append(w_i)
        self._check_homological_equation(remainder)
        logger.info('Step: %d, grade: %d', i, i+2)
        logger.info('H: %d, N: %d, R: %d',
                    len(known_terms), len(normalised), len(remainder))

        #return polynomial terms (must include factorial prefactor)
        return (1.0 / pre_factor) * k_i, (1.0 / pre_factor) * w_i

    def solve_homological_eqn(self, rem):
        """

        Find the solution of the Homological equation for given
        remainder and frequencies (_NOT_ raw eigen-values).

        This routine will always assume that the resulting inner
        product is okay; it is the responsibility of the partition
        routine to ensure that things operate correctly.  In other
        words, the partition filter, the inner-product, and the Lie
        operator, must together define a consistent normal form style.

        Therefore, at present, in the case of partial normalisation,
        there is some inefficiency here: the inner product might be
        calculated twice (once in the filter, and once again here).

        NOTES:

        We _must_ use the fundamental frequencies extracted from the
        quadratic part of the complex diagonal Hamiltonian here;
        recall that we are implementing the inverse of the Lie
        operator corresponding to this quadratic part (restricted to
        part of its image).  See also the note above regarding
        extraction of these frequencies.
        
        """
        result = self.lie_algebra().zero()
        for powers, coeff in rem.powers_and_coefficients():
            inner_prod = 0.0;
            for i in xrange(0, self.lie_algebra().dof()):
                inner_prod += (powers[2*i+1] - powers[2*i]) * self._freq[i]
            if abs(inner_prod) < 1.0e-5:
                logger.info('WARNING! small innner product %s'%inner_prod)
            result[powers] = -coeff/inner_prod
        return result

    def _check_homological_equation(self, rem):
        """

        For both debugging, and in order to cope with the natural
        reduction in numerical accuracy that occurs for large systems,
        we test that the term of the generating function, computed above,
        really does solve the corresponding homological equation for the
        given remainder terms to within acceptable numerical tolerance.

        """
        i = self._current_row
        k_2 = self._h_dc_2
        pbk = self.lie_algebra().bracket(k_2, self._w[i])
        err = (pbk + rem).l_infinity_norm()
        if 0:
            logger.info('Quadratic part')
            logger.info(self.lie_algebra().display(k_2))
        if 1:
            logger.info('Generating function term')
            logger.info(self.lie_algebra().display(self._w[i]))
        if 0:
            logger.info('Lie bracket')
            logger.info(pbk)
        logger.info('Error in homological equation: %s', err)
        assert err < (1.0 + rem.l_infinity_norm()) * 1.0e-14, err

    def _check_grade(self, i, j, h_term, w_term):
        """

        Check for consistency between the Deprit triangle indices and
        the supplied terms of the normalized function and the generating
        function.

        """
        gra0 = self.lie_algebra().grade(h_term)
        gra1 = self.lie_algebra().grade(w_term)
        assert (gra0+gra1-2) <= (i+j)+2
        if (gra0+gra1-2) != (i+j)+2:
            assert (gra0 == 0) or (gra1 == 0)

    def triangle_minus_unknown_term(self, i, j):
        """

        Compute one of the intermediate Deprit triangle quantities,
        $H_{i}^{j}$, ensuring that we do _not_ include the term from
        the Deprit recursion relation that would give rise to a Lie
        bracket between the quadratic part $H_{0}$ and the unknown
        generating function term $w_{i+j}$.

        The reason for omitting this term is so that we can handle it
        separately by solving the appropriate homological equation.

        Since this routine will be called recursively, we add an
        important assertion to ensure that the term we omit is,
        indeed, an unknown term (we do this by checking the length of
        the list representing the generating function).

        We always get the new generating function term at the
        last-but-one on the right of the lie triangle row for the
        current grade, enabling us to make the final step.

        NOTES:

        (1) that we must correct the row of the Deprit triangle later
            by removing the remainder terms from the final Lie
            bracket.

        (2) the following:

        H_i^j
        = H_{i+1}^{j-1}
        + \sum_{k=0}^i \bc{i}{k}\pb{H_{i-k}^{j-1}}{w_{k+1}}

        gra 2:
        w_0 - zero
        H_0^0 - given

        gra 3:
        w_1 - unknown
        H_1^0 - given
        H_0^1 = H_1^0 + (0,0){H_0^0, w_1} - partition left part, solve w_1

        gra 4:
        w_2 - unknown
        H_2^0 - given
        H_1^1 = H_2^0 + (1,0){H_1^0, w_1} + (1,1){H_0^0, w_2} - get w_2
        H_0^2 = H_1^1 + (0,0){H_0^1, w_1}

        """
        db_key = _make_db_key(i, j)
        alg = self.lie_algebra()
        if self._h_ij.has_key(db_key):
            pass
        else:
            if (j == 0):
                #this is more like assert (j != 0)
                assert self._h_ij.has_key(db_key)
            else:
                #we copy here for safety; ADB:
                temp = self.triangle_minus_unknown_term(i+1, j-1).copy()
                for k in xrange(0, i+1):
                    if (k+1) == (i+j):
                        #$\{H_0^0, w_{i+j}\}$
                        assert (i-k, j-1) == (0, 0)
                    else:
                        h_ik_j1 = self.triangle_minus_unknown_term(i-k, j-1)
                        w_k1 = self._w[k+1]
                        self._check_grade(i, j, h_ik_j1, w_k1)
                        temp += binomial(i, k) * alg.bracket(h_ik_j1, w_k1)
                #ONE MUST CORRECT THIS LATER BY REMOVING PART (BELOW)
                self._h_ij[db_key] = temp
        return self._h_ij[db_key]

    def _correct_row_using_partition(self, normal, remainder):
        """

        Make corrections to the row of the Deprit triangle, using the
        remainder partition of the known terms.

        NOTES:

        One the known terms have been partitioned into the normalised
        part and the remainder, we _must_ subtract away the remainder
        before proceeding; it does _not_ belong in the intermediate
        quantities of the Deprit triangle.

        This subtraction is absolutely vital!

        We could increase overall efficiency by partitioning each Hij
        contribution individually; we would not then need to subtract
        the remainder at the end.  However, in the big scheme of
        things, subtraction is pretty cheap.

        """
        i = self._current_row
        logger.info('Correcting the H_{i}^{j} using the remainder.')
        for j in xrange(1, i + 1):
            db_key = _make_db_key(i - j, j)
            #
            # we do not want to use writeback on dictionary-like shelves;
            # hence the following is _not_ an in-place operation:
            #
            self._h_ij[db_key] = self._h_ij[db_key] - remainder
        k_i = self._h_ij[_make_db_key(0, i)]
        err = (k_i - normal).l_infinity_norm()
        logger.info('Error in predicted normal form term: %s', err)
        if err > 1.0e-08:
            logger.info('WARNING!!! Error in predicted normal form term: %s',
                        err)
        assert err < 1.0e-05

    def triangle(self, i, j):
        """

        Here, we perform the Deprit triangle for a known generator.

        NOTES:

        This differs from the above method, in that the generating
        function is now known.  Thus, one does not need to perform
        partitioning, homological equation solving, or correcting of
        the row by removal of the remainder.

        """
        db_key = _make_db_key(i, j)
        alg = self.lie_algebra()
        if self._h_ij.has_key(db_key):
            pass
        else:
            if (j == 0):
                assert self._h_ij.has_key(db_key)
            else:
                #copy for safety; ADB:
                temp = self.triangle(i+1, j-1).copy()
                for k in xrange(0, i+1):
                    temp += binomial(i, k) * alg.bracket(self.triangle(i-k, j-1),
                                                         self._w[k+1])
                self._h_ij[db_key] = temp
        return self._h_ij[db_key]

