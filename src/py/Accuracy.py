"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Accuracy

PURPOSE:

Various routines associated with numerical accuracy.

NOTES:

Probably, we should log all truncations and error sizes in order to
make our computations more meaningful.

"""

from MLab import array, Complex

def trunc(c, tolerance=1.0e-15):
    """

    Take a complex number, assumed to have either tiny real or tiny
    imaginary part, and produce a pure real or pure imaginary
    approximation to it.  If neither part seems small enough to
    neglect (i.e., if either part exceeds the optional tolerance
    parameter (default tolerance=1.0e-15) then return the original
    value.

    TODO: This would be better called purify.

    """
    if isinstance(c, float):
        re, im = c, 0.0
    else:
        re, im = c.real, c.imag
    if abs(re)<tolerance:
        re = 0.0
    if abs(im)<tolerance:
        im = 0.0
    return complex(re, im)

def trunc_tuple(t, tolerance):
    """

    Apply the purify (trunc) function to a tuple or array.

    """
    return array(tuple([trunc(c, tolerance) for c in t]), Complex)
