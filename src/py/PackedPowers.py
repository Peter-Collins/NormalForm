"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: PackedPowers

PURPOSE:

Powers built using the packed tuple class.

NOTES:

This class is obsolete; we use MapPowers.

"""

from PackedTuple import *
from TuplePowers import TuplePowers

class PackedPowers(TuplePowers):
    """

    Pack powers into integers for a smaller representation using
    PackedTuple.

    Here, we override any operations to given them more efficient
    implementations on the underlying packed structure.

    """

    def __init__(self, powers=(1,), bits=4):
        assert isinstance(powers, tuple)
        self._powers = PackedTuple(powers, bits)

    def __repr__(self):
        return 'PackedPowers(%s, bits=%d)'%(repr(tuple(self._powers)),
                                            self._powers._bits)

    def __mul__(self, other):
        if len(self) == len(other):
            temp = []
            for i in range(len(self)):
                temp.append(self[i] + other[i])
            return PackedPowers(tuple(temp), bits=self._powers._bits) #ADB
        else:
            raise IndexError, "Mismatched Powers lengths"

    def __pow__(self, other):

        """Note: other must be an integer type > 0."""

        if other<0:
            raise ValueError
        temp = []
        for i in range(len(self)):
            temp.append(self[i] * other)
        return PackedPowers(tuple(temp), self._powers._bits)

    def degree(self):

        """Gives the degree of the Powers"""

        res = 0;
        #could use sum now that __iter__ is overloaded for _powers tuple?
        for i in xrange(0, len(self._powers)):
            res += self._powers[i]
        return res

    def diff(self, var):

        """Returns a pair of multiplier, differentiated monomial."""

        if var<0 or var>=len(self._powers):
            raise IndexError(var)
        t = self._powers.to_tuple()
        p = PackedPowers(t[0:var] + (max(0, t[var] - 1),) + t[var + 1:],
                         self._powers._bits)
        return t[var], p

    def to_tuple(self):
        return self._powers.to_tuple()
