
"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: EigenValueVectorPair

PURPOSE:

Trivial class to hold a pair of eigenvalue and associated eigenvector.

"""

class EigenValueVectorPair:
    """

    A simple data structure holding a single eigenvalue and its
    corresponding eigenvector.

    """

    def __init__(self, eig_val, eig_vec):
        self.val = eig_val
        self.vec = eig_vec

    def __repr__(self):
        return 'EigenValueVectorPair(%s, %s)'%(repr(self.val), repr(self.vec))
