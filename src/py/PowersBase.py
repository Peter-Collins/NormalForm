"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: PowersBase

PURPOSE:

A product of coordinates raised to exponents, represented by a tuple
of the exponents.

NOTES:

These are really multi-indices with some extra methods to allow, for
example, evaluation of the pure powers object on a vector of
coefficients.

"""

class PowersBase:

    """

    Powers interface class.

      - Represents a product of pure powers.

      - There is no leading coefficient.

    """

    #def __init__(self, *args):
    #    raise NotImplementedError, 'Abstract base'

    def __repr__(self):
        raise NotImplementedError, 'Abstract base'

    def __len__(self):
        raise NotImplementedError, 'Abstract base'

    def __cmp__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __hash__(self):
        raise NotImplementedError, 'Abstract base'

    def __eq__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __le__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __gt__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __ge__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __call__(self, args):
        raise NotImplementedError, 'Abstract base'

    def __getslice__(self, i, j):
        raise NotImplementedError, 'Abstract base'

    def __getitem__(self, i):
        raise NotImplementedError, 'Abstract base'

    def __mul__(self, other):
        raise NotImplementedError, 'Abstract base'

    def __pow__(self, other):
        raise NotImplementedError, 'Abstract base'

    def degree(self):
        raise NotImplementedError, 'Abstract base'

    def diff(self, var):
        raise NotImplementedError, 'Abstract base'

    def to_tuple(self):
        raise NotImplementedError, 'Abstract base'

    def number_of_variables(self):
        return len(self)
