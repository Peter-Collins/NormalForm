"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: SparsePowers

PURPOSE:

Implement sparse powers (multi-indices), using the sparse-tuple
representation.

NOTES:

"""

from PowersBase import PowersBase
from SparseTuple import *

class SparsePowers(PowersBase):
    """

    A sparse representation of powers, using the SparseTuple class.
    We only store the non-zero powers.

    """

    def __init__(self, *args):
        self._powers = SparseTuple(*args)

    def copy(self):
        return SparsePowers(self._powers._len, self._powers._i_to_e)

    def __repr__(self):
        return 'SparsePowers(%d, %s)'%(self._powers._len,
                                       repr(self._powers._i_to_e))

    def __mul__(self, other):
        """

        Efficient multiplication for sparse power objects.

        """
        assert isinstance(other, SparsePowers)
        if self._powers._len == other._powers._len:
            temp = dict(self._powers._i_to_e)
            for i, p in other._powers._i_to_e.iteritems():
                temp[i] = temp.get(i, 0) + p
            return SparsePowers(self._powers._len, temp)
        else:
            raise IndexError, "Mismatched Powers lengths"

    def __pow__(self, other):
        """

        Efficient powers.

        """
        if other<0:
            raise ValueError
        temp = {}
        for i, p in self._powers.iteritems():
            assert p > 0
            temp[i] = p * other
        return SparsePowers(self._powers._len, temp)

    def __len__(self):
        return len(self._powers)

    def degree(self):
        """

        Gives the degree of the Powers.  Efficient for sparse powers.

        """
        return sum(self._powers._i_to_e.itervalues())

    def diff(self, var):
        """

        Returns a pair of multiplier, differentiated monomial.

        """
        if var<0 or var>=self._powers._len:
            raise IndexError(var)
        if self._powers.has_key(var):
            d = dict(self._powers._i_to_e)
            p = self._powers[var]
            assert p > 0
            d[var] = p-1
            return float(p), SparsePowers(self._powers._len, d)
        else:
            return 0.0, SparsePowers(self._powers._len, {}) #collapse for economy

    def diff_pow(self, var, pow):
        """

        Returns a pair of multiplier, differentiated monomial.

        """
        assert isinstance(var, int)
        assert isinstance(pow, int)
        if pow <= 0:
            return 1.0, self
        assert pow > 0
        if var<0 or var>=self._powers._len:
            raise IndexError(var)
        if self._powers.has_key(var):
            p = self._powers[var]
            if pow <= p:
                d = dict(self._powers._i_to_e)
                c = p
                for i in xrange(1, pow):
                    c *= (p-i)
                if pow == p:
                    del d[var]
                else:
                    d[var] -= pow
                return float(c), SparsePowers(self._powers._len, d)
        return 0.0, SparsePowers(self._powers._len, {}) #collapse for economy

    def to_tuple(self):
        return self._powers.to_tuple()

    def __hash__(self):
        """

        This is potentially very inefficient.  In C++ we will use
        Comparators instead of hashes.

        """
        return self._powers.__hash__()

    def __eq__(self, other):
        if isinstance(other, tuple):
            return self._powers == other
        return self._powers == other._powers

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self._powers < other._powers

    def __le__(self, other):
        return self._powers <= other._powers

    def __gt__(self, other):
        return self._powers > other._powers        

    def __ge__(self, other):
        return self._powers >= other._powers

    def __getitem__(self, i):
        return self._powers[i]

    def __getslice__(self, i, j):
        return self._powers.__getslice__(i, j)

    def __call__(self, args):
        """

        Call on numerical arguments.

        """
        if len(self._powers) == len(args):
            result = 1.0
            for i, p in self._powers._i_to_e.iteritems():
                result *= args[i] ** p
            return result
        else:
            raise IndexError, "Wrong number of arguments"
