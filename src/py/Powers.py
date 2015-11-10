"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Powers

PURPOSE:

At present, this is a nasty mechanism for switching between powers
representations.

NOTES:

A program making use of powers imports this class in order to ensure a
single powers representation is in use in the program.  The need for
this arose when implementing hash for storage of disparate powers
representations in the same dictionary, or comparison based on hash
values.

There must be cleaner ways to accomplish this!

"""

from PowersBase import PowersBase

if 0:
    from TuplePowers import TuplePowers as Powers
else:
    from SparsePowers import SparsePowers as Powers
