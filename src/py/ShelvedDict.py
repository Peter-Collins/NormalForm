"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: ShelvedDict

PURPOSE:

Implement a dictionary-like object which resides in a file (shelf).

NOTES:

It is important to note that we do _NOT_ use the writeback flag on the
shelf.  Thus, in-place operations on dictionary access will not act as
expected.  One should therefore be careful with using these objects as
dict-like; they should be viewed as almost-dict-like!

For example, d[k] += v should be rewritten as d[k] = d[k] + v in the
case where the type of v is a mutable.  For immutable value types,
there is no difference between using ShelvedDict and an ordinary dict.

"""

import cPickle
import shelve
import os
import Uuid

def _key_to_str(key_obj):
    return cPickle.dumps(key_obj, -1)
def _str_to_key(str_obj):
    return cPickle.loads(str_obj)

class ShelvedDict:

    """

    We generate filenames for temporary files via UUID calls.
    We will only ask for UUIDs once per running python instance.
    We will increment a counter for the rest of the filename.

    If we choose to re-use filenames, to avoid creating too many temp
    file, then we *must* delete the corresponding files during
    __init__, otherwise their former contents will interfere.

    It is important to note that we do not use writeback in the shelf!
    This means that to change a value one should do shelf[key] = newval
    *not* shelf[key].mutate(), where mutate is a method of a mutable
    existing value.  (Of course, if we are storing immutables, then this
    does not arise.)

    """

    _count = 0
    _frees = []
    _db_name_pattern = '/tmp/shelved-dict--%s--%%d.db'%Uuid.uuidgen()

    def __init__(self, accesses_per_sync=10000):
        if ShelvedDict._frees:
            self._id = ShelvedDict._frees.pop()
        else:
            self._id = ShelvedDict._count
            ShelvedDict._count += 1
        db_name = (ShelvedDict._db_name_pattern)%self._id
        if os.path.exists(db_name):
            os.remove(db_name) #ensure a clean file!
        self._shelf = shelve.open(db_name, writeback=False, protocol=-1)
        assert not self._shelf.keys()
        self._accesses_per_sync = accesses_per_sync
        self._access = 0

    def __del__(self):
        ShelvedDict._frees.append(self._id)
        self._shelf.close()
        del self._shelf
        if self._access > self._accesses_per_sync:
            msgs = ['Shelf accesses: %d'%self._access,
                    'Per sync: %d'%self._accesses_per_sync,
                    'Syncs: %d'%self._access/self._accesses_per_sync]
            print ', '.join(msgs)

    def __getitem__(self, key_obj):
        #this returns a copy, so take care with mutable values
        return self._shelf[_key_to_str(key_obj)]

    def __setitem__(self, key_obj, value):
        self._shelf[_key_to_str(key_obj)] = value
        self._access += 1
        if self._access > self._accesses_per_sync:
            self._access = 0
            self._shelf.sync()

    def has_key(self, key):
        return self._shelf.has_key(_key_to_str(key))

    def iterkeys(self):
        for key in self._shelf.keys():
            yield _str_to_key(key)

    def keys(self):
        return list(self.iterkeys())
        
    def iteritems(self):
        for key, value in self._shelf.iteritems():
            yield (_str_to_key(key), value)

    def items(self):
        return list(self.iteritems())

    def itervalues(self):
        for key, value in self._shelf.iteritems():
            yield value

    def values(self):
        return list(self.itervalues())

    def __len__(self):
        return len(self._shelf)

    def get(self, key, default):
        return self._shelf.get(_key_to_str(key), default)

    def __delitem__(self, key):
        del self._shelf[_key_to_str(key)]

    def __str__(self):
        #needs to strip quotes
        return str(dict(self.iteritems()))

    def __repr__(self):
        #needs to strip quotes
        return repr(dict(self.iteritems()))

    def __eq__(self, other):
        if len(self.keys()) != len(other.keys()):
            return 0
        #we cannot use == here, because UserDict calls cmp
        #and complex numbers are not ordered!
        for k, v in self.iteritems():
            if not (other.has_key(k)):
                return 0
            if v != other[k]:
                return 0
        return 1

    def __ne__(self, other):
        return not self.__eq__(other)


