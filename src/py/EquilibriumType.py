"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: EquilibriumType

PURPOSE:

Define shorthand for equilibrium point types and provide
pretty-printing.

NOTES:

"""

eq_type_names = {'c': 'centre',
                 's': 'saddle',
                 '0': 'zero'}

def equilibrium_type_to_tex(eq_type):
    """

    @param eq_type: a string of s, c, 0.

    @return: a TeX representation of the type of the equilibrium
    point.

    """
    s_list = []
    names = eq_type_names
    count = 0
    prev = None
    are_saddles = list(eq_type)
    are_saddles.append(-1)
    for tp in are_saddles:
        if not (prev == None):
            if tp != prev:
                assert count > 0
                s = names[prev]
                if count > 1:
                    s = '%s$^{%d}$'%(s, count)
                s_list.append(s)
                count = 0
        count += 1
        prev = tp
    return '$\\times$'.join(s_list)

def equilibrium_type_to_str(eq_type):
    """

    @param eq_type:

    @return: a plaintext representation of the type of the equilibrium
    point.

    """
    s_list = []
    names = eq_type_names
    count = 0
    prev = None
    are_saddles = list(eq_type)
    are_saddles.append(-1)
    for tp in are_saddles:
        if not (prev == None):
            if tp != prev:
                assert count > 0
                s = names[prev]
                if count > 1:
                    s = '%s^(%d)'%(s, count)
                s_list.append(s)
                count = 0
        count += 1
        prev = tp
    return ' x '.join(s_list)

