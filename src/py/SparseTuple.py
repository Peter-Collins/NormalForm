"""
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

A sparse representation of (non-negative) integer tuples; we only
store the non-zero elements.

"""

from math import *

class SparseTupleError(Exception):
    pass

class NegativeIntegerError(SparseTupleError):
    pass

class SparseTuple:

    """

    Sparse integer tuples.  Tuple of ints, only storing non-zero
    values.

      - stored as dictionary from index to non-zero value.

    """

    def __init__(self, *args):
        if len(args) <= 1:
            self.__init__from_elements(*args)
        else:
            if len(args) == 2:
                self.__init__from_length_and_dict(*args)
            else:
                raise IndexError, 'Bad SparsePowers constructor call'

    def __init__from_length_and_dict(self, length, ind_to_elt=None):
        """

        Usual constructor: construct SparseTuple from length and map
        from index to element.

        @param length: length of the tuple to be represented.
        @param ind_to_elt: map from index to (positive) element.

        """
        assert isinstance(length, int)
        assert length >= 0
        self._len = length
        self._i_to_e = {}
        if ind_to_elt:
            for i, e in ind_to_elt.iteritems():
                assert i >= 0 and i < length
                if e < 0:
                    raise NegativeIntegerError('negative value in sparse tuple')
                if e != 0:
                    self._i_to_e[i] = e

    def __init__from_elements(self, elements=(0,)):
        """

        Constructor: Create SparseTuple from iterable of elements.

        @param elements: non-negative integers.
        @type elements: iterable of non-negative integers.

        It would be much more efficient to be able to create a sparse
        tuple from a dictionary?  There is too much converting of
        lists to tuples to sparse tuple dicts going on.

        """
        self._i_to_e = {}
        self._len = 0
        for ind, elt in enumerate(elements):
            if elt < 0:
                raise NegativeIntegerError('negative value in sparse tuple')
            if elt != 0:
                self._i_to_e[ind] = elt
            self._len += 1

    def iteritems(self):
        return self._i_to_e.iteritems()

    def has_key(self, key):
        return self._i_to_e.has_key(key)

    def __getitem__(self, j):
        if isinstance(j, slice):
            return self.__getitem__slice(j)
        if j < 0:
            j += self._len
        if (j >= 0) and (j < self._len):
            return self._i_to_e.get(j, 0)
        raise IndexError('sparse tuple index out of range')

    def __getitem__slice(self, sl):
        start, stop, step = sl.start, sl.stop, sl.step
        res = {}
        ind = 0
        if start == None:
            start = 0
        if stop == None:
            stop = self._len
        if step == None:
            step = 1
        for k in xrange(start, stop, step):
            if self._i_to_e.has_key(k):
                res[ind] = self._i_to_e[k]
            ind += 1
        return SparseTuple(int(self._len/step), res)

    def __getslice__(self, start, stop):
        return self.__getitem__slice(slice(start, stop, 1))
        
    def __setitem__(self, j, item):
        if (j >= 0) and (j < len(self)):
            if item < 0:
                raise NegativeIntegerError(item)
            if self._i_to_e.has_key(j) and item == 0:
                del self._i_to_e[j] #remove zeros
            else:
                if item != 0:
                    self._i_to_e[j] = item
        else:
            raise IndexError('sparse tuple index out of range')

    def to_list(self):
        """

        Conversion to a list.

        """
        res = [0,]*self._len
        for i, p in self._i_to_e.iteritems():
            res[i] = p
        return res

    def to_tuple(self):
        """

        Conversion to a tuple.

        """
        return tuple(self.to_list())
    
    def __repr__(self):
        return 'SparseTuple(%d, %s)'%(self._len, self._i_to_e)

    def __cmp__(self, other):
        assert 0, 'use rich comparison only!'
        if len(self) == len(other):
            return self.to_tuple() == other.to_tuple()
        else:
            return False

    def __len__(self):
        return self._len

    def __hash__(self):
        """

        THIS IS TERRIBLY INEFFICIENT IF WE USE CONVERSION TO TUPLE?

        If we don't do this, though, we end up with objects that don't
        compare as keys to their TuplePowers counterparts.  For now,
        we'll take that risk, and hope that only one kind of power is
        in use at any time.

        Equal things must have same hash.  This is the only
        restriction; hash is used first, then equality when needed.

        """
        return hash(tuple(self._i_to_e.iterkeys()) + tuple(self._i_to_e.itervalues()))
        #return hash(self.to_tuple())

    def __iter__(self):
        for i in xrange(self._len):
            yield self[i]

    def __str__(self):
        if self._len == 1:
            res = str(self[0])+','
        else:
            res = ', '.join([str(x) for x in self])
        return '('+res+')'

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        if isinstance(other, tuple):
            return self.to_tuple() == other
        assert isinstance(other, SparseTuple)
        return self._i_to_e == other._i_to_e

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.to_tuple() < other.to_tuple()

    def __le__(self, other):
        return self.to_tuple() <= other.to_tuple()

    def __gt__(self, other):
        return self.to_tuple() > other.to_tuple()

    def __ge__(self, other):
        return self.to_tuple() >= other.to_tuple()

