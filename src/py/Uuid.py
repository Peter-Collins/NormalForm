"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Uuid

PURPOSE:

Generate unique identifiers for use as temporary file name prefixes.

NOTES:

This is important!  When using temporary files, we need to ensure
uniqueness among single threads in an application, and among multiple
copies of the same executable.

We use the linux "uuidgen" operating system call.

"""

import commands

def uuidgen():
    return commands.getoutput('uuidgen')

