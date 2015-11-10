"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: PackedTuple

PURPOSE:

Tuples of Ints packed together

NOTES:

This class is obsolete; we use hash-based maps.

16 = 0001 0000
15 = 1111
14 = 1110
13 = 1101
12 = 1100
11 = 1011
10 = 1010
9  = 1001
8  = 1000
7  = 0111
6  = 0110
5  = 0101
4  = 0100
3  = 0011
2  = 0010
1  = 0001
0  = 0000

17185 = 0100|0011|0010|0001
801   = 0000|0011|0010|0001
33    = 0000|0000|0010|0001
1     = 0000|0000|0000|0001

"""

from math import *

class PackedTupleError(Exception): pass
class NegativeIntegerError(PackedTupleError): pass
class InsufficientBitsError(PackedTupleError): pass
class NonPositiveBitsError(PackedTupleError): pass
class TooManyBitsError(PackedTupleError): pass

class PackedTuple:

    """

    Integer tuples packed to save space.

    """

    def __init__(self, terms, bits=4):

        """
        bits is the size to use for packing
        then maximum size of an entry in the tuple is 2^bits
        terms is an input integer tuple
        The packed numbers are stored in a Long integer
        n = len / npi
        """

        #assert isinstance(terms, tuple) #any sequence type okay
        if bits<=0:
            raise NonPositiveBitsError()
        if bits>30:
            raise TooManyBitsError()
        self._bits = bits
        max_size = self.max_size()
        self._len = len(terms)
        num_per_int = 31 / bits
        n = int(ceil(float(self._len) / float(num_per_int)))
        self._val = []
        for j in xrange(0, n):
            v = 0
            for i in xrange(0, num_per_int):
                m = j * num_per_int + i
                if (m < self._len):
                    if terms[m]<0:
                        msg = 'term: %s'%terms[m]
                        raise NegativeIntegerError(msg)
                    if terms[m]>max_size:
                        msg = 'term: %s bits: %s'%(terms[m], bits)
                        raise InsufficientBitsError(msg)
                    v += (terms[m] << (i * self._bits))
                else:
                    break
            self._val.append(v)

    def max_size(self):
        """ The maximum size for any term in the tuple """
        return (1 << self._bits) - 1

    def __getitem__(self, i):
        if (i >= 0) & (i < len(self)):
            num_per_int = 31 / self._bits
            n = i / num_per_int
            v = self._val[n]
            j = i - n * num_per_int
            return (v >> (j * self._bits)) & ((1 << self._bits) - 1)
        else:
            raise IndexError('packed tuple index out of range')

    def __setitem__(self, i, item):
        if (i >= 0) & (i < len(self)):
            if item<0:
                raise NegativeIntegerError(item)
            if item>self.max_size():
                raise InsufficientBitsError(item)
            num_per_int = 31 / self._bits
            n = i / num_per_int
            j = i - n * num_per_int
            cur_ival = (self[i] << (j * self._bits))
            ival = item << (j * self._bits)
            self._val[n] -= cur_ival
            self._val[n] += ival
        else:
            raise IndexError('packed tuple index out of range')

    def n_bits(self):
        return self._bits
        
    def to_list(self):
        res = []
        for i in xrange(0,self._len):
            res.append(self[i])
        return res

    def to_tuple(self):
        return tuple(self.to_list())
    
    def __getslice__(self, i, j):
        return PackedTuple(tuple((self.to_list())[i:j]), self._bits)
        
    def __repr__(self):
        if self._len == 1:
            res = '('+str(self[0])+',)'
        else:
            res = '('+(', '.join([str(x) for x in self]))+')' 
        return 'PackedTuple(%s, bits=%d)'%(res, self._bits)

    def __cmp__(self, other):
        if len(self) == len(other):
            #    return cmp(self._terms, other._terms)
            return self.to_tuple() == other.to_tuple()
        else:
            return False

    def __len__(self):
        return self._len

    def __hash__(self):
        return hash(self.to_tuple())

    def __iter__(self):
        for i in xrange(self._len):
            yield self[i]

    def __str__(self):
        """Convert to string.  ADB: take care with one-element tuples!"""
        if self._len == 1:
            res = str(self[0])+','
        else:
            res = ', '.join([str(x) for x in self])
        return '('+res+')'

    def __eq__(self, other):
        if isinstance(other, tuple):
            return self.to_tuple() == other
        return self.to_tuple() == other.to_tuple()

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
