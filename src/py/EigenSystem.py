"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: EigenSystem

PURPOSE:

Eigensystems and Diagonalization of the quadratic part of the
Hamiltonian.

"""

import logging
logger = logging.getLogger() #'Diagonal')

from math import sqrt
from Accuracy import trunc, trunc_tuple
from Complex import real_part_of_vector
from LieAlgebra import LieAlgebra
from EigenValueVectorPair import EigenValueVectorPair
## Automatically adapted for numpy.oldnumeric Dec 16, 2008 by alter_code1.py
try:
    from numpy.oldnumeric.mlab import array, Complex, Float, zeros, eig
    from numpy.oldnumeric.linear_algebra import determinant, inverse
    from numpy.oldnumeric import transpose, matrixmultiply
except:
    from MLab import array, Complex, Float, zeros, eig
    from LinearAlgebra import determinant, inverse
    from Numeric import transpose, matrixmultiply
from Polynomial import Polynomial


def sort_by_ranking_function(seq, func=None):
    """

    @param seq: an iterable object.
    @param func: a function which returns a sort key for each item in seq.

    @return: a sorted list of the elements in seq, sorted with respect
    to the computed keys.

    """
    if func==None:
        def func(item):
            return item
    decorated = [(func(item), i, item) for i, item in enumerate(seq)]
    decorated.sort()
    return [dec[-1] for dec in decorated]

def rank_real_before_imag(x):
    """

    Sort complex numbers, real before imaginary, decreasing real
    magnitude, decreasing imaginary magnitude, pairing-up
    equal-magnitude but opposite-sign pairs with the positive partner
    first.

    """
    val, vec = x.val, x.vec
    c = complex(val)
    re, im = c.real, c.imag
    return (-abs(re), -abs(im), -re, -im) + tuple(real_part_of_vector(vec))

def rank_eigen_value_vector_pair(x):
    """

    Rank an EigenValueVectorPair by its eigenvalue, with reals before
    pure imaginaries.

    """
    return rank_real_before_imag(x)

def positive_member_first(left, right):
    """

    Given two EigenValueVectorPairs, left and right, return the pairs
    ordered so that the positive eigenvalue pair comes first.

    """
    #put +a before -a, +aJ before -aJ;
    lval, rval = left.val, right.val
    if not isinstance(lval, complex):
        lval = complex(lval, 0.0)
    if not isinstance(rval, complex):
        rval = complex(rval, 0.0)
    if lval.real == 0.0:
        assert rval.real == 0.0
        #pure-imaginary
        if lval.imag < 0.0:
            assert rval.imag > 0.0
            return (right, left)
        else:
            return (left, right)
    else:
        assert lval.imag == 0.0
        assert rval.imag == 0.0
        #pure real
        if lval.real < 0.0:
            assert rval.real > 0.0
            return (right, left)
        else:
            return (left, right)
    assert 0, 'you should never arrive here'

def is_pure_real(c):
    """

    Check that argument is a pure real (either a non-zero numerical
    float-like type, or a complex type with zero imaginary part and
    nonzero real part).

    """
    if isinstance(c, float):
        return c != 0.0
    assert isinstance(c, complex), repr(c)
    return ((c.imag == 0.0) and (c.real != 0.0))

def is_pure_imag(c):
    """

    Check that argument is a pure imaginary (a complex type with
    non-zero imaginary part and zero real part).

    """
    if isinstance(c, float):
        return False
    assert isinstance(c, complex), repr(c)
    return ((c.real == 0.0) and (c.imag != 0.0))

class EigenSystem:

    def __init__(self, eigen_value_vector_pairs, tolerance):
        self.val_vec_pairs = eigen_value_vector_pairs
        self.tolerance = tolerance

    def __repr__(self):
        return 'EigenSystem(%s, %s)'%(self.val_vec_pairs, self.tolerance)

    def analyse_eigen_system(self):
        self._truncate_eigen_values_and_vectors()
        self._sort_eigensystem_reals_first()
        self._collect_plus_minus_pairs()
        self._determine_equilibrium_type()

    def get_raw_eigen_value_vector_pairs(self):
        return self.val_vec_pairs

    def get_eigen_values_in_plus_minus_pairs(self):
        return [(p.val, m.val) for (p,m) in self.plus_minus_pairs]

    def get_plus_minus_pairs(self):
        return self.plus_minus_pairs

    def get_equilibrium_type(self):
        return self.equilibrium_type

    def _truncate_eigen_values_and_vectors(self):
        """

        Truncate the floating point representations of the eigenvalues
        and eigenvectors, so as to give eigenvalues in exact +/-
        pairs.

        """
        pairs = []
        for val_vec in self.val_vec_pairs:
            val = val_vec.val
            vec = val_vec.vec
            new_val = trunc(val, self.tolerance)
            new_vec = trunc_tuple(vec, self.tolerance)
            #logger.info('truncated eignenvalue', new_val)
            pairs.append(EigenValueVectorPair(new_val, new_vec))
        self.val_vec_pairs = pairs

    def _sort_eigensystem_reals_first(self):
        """

        Sort the eigensystem into a predictable ordering, with the
        pure real eigenvalues first.

        """
        pairs = self.val_vec_pairs
        sorted = sort_by_ranking_function(pairs, rank_eigen_value_vector_pair)
        self.val_vec_pairs = sorted

    def _collect_plus_minus_pairs(self):
        """

        Put the eigenvalues into pairs of (+a, -a) and (+bJ, -bJ).  We
        assume that two values (c,d) are really (+/- a, -/+a) if a
        differs from (-a) by less than the specified numerical
        tolerance.  This handles the siuation where LAPACK (or another
        library) does not produce exact pairs.

        """
        seq = self.val_vec_pairs
        tol = self.tolerance
        if len(seq)%2:
            raise IndexError, 'odd length sequence given to locate_pairs.'
        pairs = []
        for i_left in xrange(0, len(seq), 2):
            i_right = i_left+1
            left, right = seq[i_left], seq[i_right]
            left_val = trunc(complex(left.val), tol)
            right_val = trunc(complex(right.val), tol)
            err = abs(left_val-(-right_val))
            if err > tol:
                raise ValueError, 'eigenvalues not in +/- pairs, error %s'%err
            if (abs(left_val) == 0.0) or (abs(right_val) == 0.0):
                logger.info('zero eigenvalue found')
                if (abs(left_val) != 0.0) or (abs(right_val) != 0.0):
                    raise ValueError, 'zero not paired with zero'
            pairs.append(positive_member_first(left, right))
        self.plus_minus_pairs = pairs
        self._lie = LieAlgebra(len(pairs))

    def _determine_equilibrium_type(self):
        """

        Under the assumption the the eigensystem corresponds to an
        equilibrium point of saddles x centres type, determine the
        type of the equilibrium point.  The ordering here is that same
        as the ordering of the eigensystem into +/- pairs.

        """
        eq_type = ''
        for pos, neg in self.plus_minus_pairs:
            if is_pure_real(pos.val):
                assert is_pure_real(neg.val)
                eq_type += 's'
                #logger.info('saddle')
            else:
                if is_pure_imag(pos.val):
                    assert is_pure_imag(neg.val)
                    eq_type += 'c'
                    #logger.info('centre')
                else:
                    if (pos.val == 0.0):
                        assert (neg.val == 0.0)
                        eq_type += '0'
                        #logger.info('zero')
                    else:
                        raise ValueError, 'unknown eigenvalue combination.'
        assert len(eq_type) == len(self.plus_minus_pairs)
        logger.info('equilibrium type %s', eq_type)
        self.equilibrium_type = eq_type

# def sort_by_ranking_function(seq, func=None):
#     """

#     @param seq: an iterable object.
#     @param func: a function which returns a sort key for each item in seq.

#     @return: a sorted list of the elements in seq, sorted with respect
#     to the computed keys.

#     """
#     if func==None:
#         def func(item):
#             return item
#     decorated = [(func(item), i, item) for i, item in enumerate(seq)]
#     decorated.sort()
#     return [dec[-1] for dec in decorated]

# def rank_real_before_imag(x):
#     """

#     Sort complex numbers, real before imaginary, decreasing real
#     magnitude, decreasing imaginary magnitude, pairing-up
#     equal-magnitude but opposite-sign pairs with the positive partner
#     first.

#     """
#     c = complex(x)
#     re, im = c.real, c.imag
#     return (-abs(re), -abs(im), -re, -im)

# def rank_by_eigenvalue(x):
#     """

#     Rank an EigenValueVectorPair by its eigenvalue, with reals before
#     pure imaginaries.

#     """
#     return rank_real_before_imag(x.val)

# def positive_member_first(left, right):
#     """

#     Given two EigenValueVectorPairs, left and right, return the pairs
#     ordered so that the positive eigenvalue pair comes first.

#     """
#     #put +a before -a, +aJ before -aJ;
#     lval, rval = left.val, right.val
#     if not isinstance(lval, complex):
#         lval = complex(lval, 0.0)
#     if not isinstance(rval, complex):
#         rval = complex(rval, 0.0)
#     if lval.real == 0.0:
#         assert rval.real == 0.0
#         #pure-imaginary
#         if lval.imag < 0.0:
#             assert rval.imag > 0.0
#             return (right, left)
#         else:
#             return (left, right)
#     else:
#         assert lval.imag == 0.0
#         assert rval.imag == 0.0
#         #pure real
#         if lval.real < 0.0:
#             assert rval.real > 0.0
#             return (right, left)
#         else:
#             return (left, right)
#     assert 0, 'you should never arrive here'

# def is_pure_real(c):
#     """

#     Check that argument is a pure real (either a non-zero numerical
#     float-like type, or a complex type with zero imaginary part and
#     nonzero real part).

#     """
#     if isinstance(c, float):
#         return c != 0.0
#     assert isinstance(c, complex), repr(c)
#     return ((c.imag == 0.0) and (c.real != 0.0))

# def is_pure_imag(c):
#     """

#     Check that argument is a pure imaginary (a complex type with
#     non-zero imaginary part and zero real part).

#     """
#     if isinstance(c, float):
#         return False
#     assert isinstance(c, complex), repr(c)
#     return ((c.real == 0.0) and (c.imag != 0.0))

# def real_part_of_scalar(sca):
#     if not isinstance(sca, complex):
#         return sca
#     return sca.real

# def imag_part_of_scalar(sca):
#     if not isinstance(sca, complex):
#         return 0.0
#     return sca.imag

# def real_part_of_vector(vec):
#     return array(tuple([real_part_of_scalar(x) for x in vec]), Float)

# def imag_part_of_vector(vec):
#     return array(tuple([imag_part_of_scalar(x) for x in vec]), Float)

# def realify_vector(vec):
#     return real_part_of_vector(vec)+imag_part_of_vector(vec)

# class EigenSystem:

#     def __init__(self, eigen_value_vector_pairs, tolerance):
#         self.val_vec_pairs = eigen_value_vector_pairs
#         self.tolerance = tolerance

#     def __repr__(self):
#         return 'EigenSystem(%s, %s)'%(self.val_vec_pairs, self.tolerance)

#     def analyse_eigen_system(self):
#         self._truncate_eigen_values_and_vectors()
#         self._sort_eigensystem_reals_first()
#         self._collect_plus_minus_pairs()
#         self._determine_equilibrium_type()

#     def get_raw_eigen_value_vector_pairs(self):
#         return self.val_vec_pairs

#     def get_eigen_values_in_plus_minus_pairs(self):
#         return [(p.val, m.val) for (p,m) in self.plus_minus_pairs]

#     def get_plus_minus_pairs(self):
#         return self.plus_minus_pairs

#     def get_equilibrium_type(self):
#         return self.equilibrium_type

#     def _truncate_eigen_values_and_vectors(self):
#         """

#         Truncate the floating point representations of the eigenvalues
#         and eigenvectors, so as to give eigenvalues in exact +/-
#         pairs.

#         """
#         pairs = []
#         for val_vec in self.val_vec_pairs:
#             val = val_vec.val
#             vec = val_vec.vec
#             new_val = trunc(val, self.tolerance)
#             new_vec = trunc_tuple(vec, self.tolerance)
#             logger.info('truncated eignenvalue', new_val)
#             pairs.append(EigenValueVectorPair(new_val, new_vec))
#         self.val_vec_pairs = pairs

#     def _sort_eigensystem_reals_first(self):
#         """

#         Sort the eigensystem into a predictable ordering, with the
#         pure real eigenvalues first.

#         """
#         pairs = self.val_vec_pairs
#         sorted = sort_by_ranking_function(pairs, rank_by_eigenvalue)
#         self.val_vec_pairs = sorted

#     def _collect_plus_minus_pairs(self):
#         """

#         Put the eigenvalues into pairs of (+a, -a) and (+bJ, -bJ).  We
#         assume that two values (c,d) are really (+/- a, -/+a) if a
#         differs from (-a) by less than the specified numerical
#         tolerance.  This handles the siuation where LAPACK (or another
#         library) does not produce exact pairs.

#         """
#         seq = self.val_vec_pairs
#         tol = self.tolerance
#         if len(seq)%2:
#             raise IndexError, 'odd length sequence given to locate_pairs.'
#         pairs = []
#         for i_left in xrange(0, len(seq), 2):
#             i_right = i_left+1
#             left, right = seq[i_left], seq[i_right]
#             left_val = trunc(complex(left.val), tol)
#             right_val = trunc(complex(right.val), tol)
#             err = abs(left_val-(-right_val))
#             if err > tol:
#                 raise ValueError, 'eigenvalues not in +/- pairs, error %s'%err
#             if (abs(left_val) < tol) or (abs(right_val) < tol):
#                 raise ValueError, 'eigenvalues closer to zero than tolerance.'
#             pairs.append(positive_member_first(left, right))
#         self.plus_minus_pairs = pairs
#         self._lie = LieAlgebra(len(pairs))

#     def _determine_equilibrium_type(self):
#         """

#         Under the assumption the the eigensystem corresponds to an
#         equilibrium point of saddles x centres type, determine the
#         type of the equilibrium point.  The ordering here is that same
#         as the ordering of the eigensystem into +/- pairs.

#         """
#         eq_type = []
#         for pos, neg in self.plus_minus_pairs:
#             if is_pure_real(pos.val):
#                 assert is_pure_real(neg.val)
#                 eq_type.append(True)
#                 logger.info('saddle')
#             else:
#                 if is_pure_imag(pos.val):
#                     assert is_pure_imag(neg.val)
#                     eq_type.append(False)
#                     logger.info('centre')
#                 else:
#                     raise ValueError, 'eigenvalue neither real nor pure imag.'
#         assert len(eq_type) == len(self.plus_minus_pairs)
#         self.equilibrium_type = eq_type
