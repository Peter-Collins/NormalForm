"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: GradedInterface

PURPOSE:

Interface for graded algebras.

"""

class GradedInterface:
    def grade(self, elt):
        raise NotImplementedError
    def is_isograde(self, elt, grade=-1):
        raise NotImplementedError
    def isograde(self, elt, grade, up_to=None):
        raise NotImplementedError

