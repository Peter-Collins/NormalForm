"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: TuplePowers

PURPOSE:

Implement powers (multi-indices) by a tuple of the exponents; this is
the so-called dense representation.

NOTES:

"""

from PowersBase import PowersBase

class TuplePowers(PowersBase):

    """

    Powers class.

      - Represents a product of pure powers.

      - There is no leading coefficient.

    """

    def __init__(self, *args):
        if len(args) <= 1:
            self.__init__from_elements(*args)
        else:
            if len(args) == 2:
                self.__init__from_length_and_dict(*args)
            else:
                raise IndexError, 'Bad Powers constructor call'

    def __init__from_elements(self, powers=(0,)):
        assert isinstance(powers, tuple)
        self._powers = powers

    def __init__from_length_and_dict(self, length, ind_to_pow):
        assert isinstance(length, int)
        assert isinstance(ind_to_pow, dict)
        powers = [0,]*length
        for i, p in ind_to_pow.iteritems():
            assert p >= 0
            if p != 0:
                powers[i] = p
        self._powers = tuple(powers)

    def __repr__(self):
        return 'TuplePowers(%s)'%repr(self._powers)

    def __cmp__(self, other):
        if len(self) == len(other):
            return self._powers == other._powers
        else:
            return False

    def __len__(self):
        return len(self._powers)

    def __hash__(self):
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

    def __call__(self, args):
        "args is a tuple that must be the same length as the monomial"
        if len(self._powers) == len(args):
            result = 1.0
            for i in range(len(args)):
                result *= args[i] ** self._powers[i]
            return result
        else:
            raise IndexError, "Wrong number of arguments"

    def __getslice__(self, i, j):
        return self._powers[i:j]

    def __getitem__(self, i):
        return self._powers[i]

    def __mul__(self, other):
        if not isinstance(other, TuplePowers):
            raise TypeError
        if len(self) == len(other):
            temp = []
            for i in range(len(self)):
                temp.append(self[i] + other[i])
            return TuplePowers(tuple(temp))
        else:
            raise IndexError, "Mismatched Powers lengths"

    def __pow__(self, other):
        "other must be an integer type > 0"
        if other<0:
            raise ValueError
        temp = []
        for i in range(len(self)):
            temp.append(self[i] * other)
        return TuplePowers(tuple(temp))

    def degree(self):
        "Gives the degree of the Powers"
        return sum(self._powers)

    def diff(self, var):
        if var<0 or var>=len(self._powers):
            raise IndexError
        return self._powers[var], TuplePowers(self._powers[0:var] + (max(0, self._powers[var] - 1),) + self._powers[var + 1:])

    def to_tuple(self):
        return tuple(self._powers)
