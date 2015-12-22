"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Complex

PURPOSE:

Implement a few helper methods for complex numbers and vectors.

NOTES:

We provide a pretty-printer for complex numbers that provides a
fixed-width format leaving spaces for zero real and imaginary parts.
This enables fast visual inspection of tables of coefficients.

"""

## Automatically adapted for numpy.oldnumeric Dec 16, 2008 by alter_code1.py
try:
    from numpy.oldnumeric.mlab import array, Float
except:
    from MLab import array, Float

def pretty_complex(com):
    """

    Pretty-printing for complex numbers.

    """
    co = complex(com)
    if co.imag == 0.0:
        res = ('(%+1.12e '%co.real+' '*(12+7)+' )')
    else:
        if co.real == 0.0:
            res = ('('+' '*(12+7)+' %+1.12ej)'%co.imag)
        else:
            res = ('(%+1.12e %+1.12ej)'%(co.real, co.imag))
    return res

def real_part_of_scalar(sca):
    """Extract real part of complex scalar."""
    if not isinstance(sca, complex):
        return sca
    return sca.real

def imag_part_of_scalar(sca):
    """Extract imaginary part of complex scalar."""
    if not isinstance(sca, complex):
        return 0.0
    return sca.imag

def real_part_of_vector(vec):
    """Extract real vector part of complex vector."""
    return array(tuple([real_part_of_scalar(elt) for elt in vec]), Float)

def imag_part_of_vector(vec):
    """Extract pure imaginary part of complex vector."""
    return array(tuple([imag_part_of_scalar(elt) for elt in vec]), Float)

def realify_vector(vec):
    """Is this used anywhere?"""
    return real_part_of_vector(vec)+imag_part_of_vector(vec)

