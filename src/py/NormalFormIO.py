"""

AUTHOR: Peter Collins, 2006.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: NormalFormIO.py

PURPOSE:

Collect together software to read in and write Normal Forms
both from C++ and as ascii files from Mathematica.

NOTES:
    Moved from Polynomial.py and expanded to matrices and vectors
    We make xxpp an IO problem by always using xpxp internally
    Mathematica also uses different order for the planes
    Thus collect it here, maybe move low level bits in the future
"""

from Powers import Powers
from Polynomial import Polynomial

def xxpp_to_xpxp(xxpp):
    """

    Cope with different ordering conventions for coordinates and
    their canonically-conjugate momenta.  This more properly belongs
    in the LieAlgebra module.

    """
    assert not len(xxpp)%2
    dof = len(xxpp)/2
    xpxp = []
    for x, p in zip(xxpp[:dof], xxpp[dof:]):
        xpxp.append(x)
        xpxp.append(p)
    return xpxp

def xpxp_to_xxpp(xpxp):
    """

    Cope with different ordering conventions for coordinates and
    their canonically-conjugate momenta.  This more properly belongs
    in the LieAlgebra module.

    """
    assert not len(xpxp)%2
    return list(xpxp[0::2])+list(xpxp[1::2])

def reorder_dof(inp, order):
    """
    reorder the degrees of freedom on the fly
    """
    out = []
    n_ints = len(order)
    if len(inp) == n_ints:
        # only dof terms, e.g. integrals
        for i in xrange(n_ints):
            out.append(inp[order[i]])
    else:
        # xpxp list
        for i in xrange(n_ints):
            out.append(inp[order[i]*2])
            out.append(inp[order[i]*2 +1])
        # check for Semi-Classical
        if len(inp) == n_ints*2+1:
            out.append(inp[n_ints*2 +1])
        else:
            assert len(inp) == n_ints*2
    return out

def read_ascii_polynomial(istream, is_xxpp_format, order=None):
    """

    Read a polynomial from an input stream in ASCII format.

    @param istream: an iterable of input lines (e.g. file-like),

    @param is_xxpp_format: a boolean to specify whether the input file
    is in the old (Mathematica) format where all configuration
    coordinates come first, following by all the momenta, rather than
    the new even-odd format.

    @param order: reorder the degrees of freedom on the fly

    @return: polynomial.

    Note: this needs to be refactored; the xxpp logic belongs in a
    lie algebra of one sort or another.

    """
    line = istream.next()
    n_vars = int(line)
    assert n_vars >= 0
    line = istream.next()
    n_monomials = int(line)
    assert n_monomials >= 0
    p = Polynomial(n_vars)
    for i in xrange(n_monomials):
        line = istream.next()
        elts = line.split(' ')
        elts = elts[:n_vars]+[' '.join(elts[n_vars:])]
        assert len(elts) == n_vars+1
        powers_list = [int(e) for e in elts[:-1]]
        if is_xxpp_format:
            powers_list = xxpp_to_xpxp(powers_list)
        if order != None:
            powers_list = reorder_dof(powers_list, order)
        powers = tuple(powers_list)
        coeff = complex(elts[-1])
        m = coeff*Polynomial.Monomial(powers)
        assert not p.has_term(Powers(powers))
        p += m
    return p

def write_ascii_polynomial(ostream, poly, is_xxpp_format):
    """

    Real a polynomial from an input stream in ASCII format.

    @param ostream: an iterable of input lines (e.g. file-like),

    @param poly: a polynomial,

    @param is_xxpp_format: a boolean to specify whether the output
    file should be in the old (Mathematica) format where all
    configuration coordinates come first, following by all the
    momenta, rather than the new even-odd format.

    """
    ostream.write('%d\n'%poly.n_vars())
    ostream.write('%d\n'%len(poly))
    for pows, coeff in poly.powers_and_coefficients():
        powers_list = pows
        if is_xxpp_format:
            powers_list = xpxp_to_xxpp(powers_list)
        for p in powers_list:
            ostream.write('%d '%p)
        if isinstance(coeff, complex):
            ostream.write('%.17g%+.17gj\n'%(coeff.real, coeff.imag))
            #ostream.write('%s %+sj\n'%(repr(coeff.real), repr(coeff.imag)))
        else:
            ostream.write('%.17g\n'%coeff)
            #ostream.write('%s\n'%repr(coeff))

def write_file_polynomial(file_name, poly, is_xxpp_format):
    file_ostr = open(file_name, 'w')
    write_ascii_polynomial(file_ostr, poly.real(), is_xxpp_format)
    file_ostr.close()
    
def read_ascii_matrix(file_istr, is_xxpp_format):
    """
    Read a matrix from an input stream in ASCII format.
    """
    line = file_istr.next()
    n_rows = int(line)
    assert n_rows >= 0
    line = file_istr.next()
    n_cols = int(line)
    assert n_cols >= 0
    matrix = []
    for i in xrange(n_rows):
        line = file_istr.next()
        elts = line[:-1].split(' ')
        row = []
        for e in xrange(n_cols):
            row.append(float(elts[e]))
        if is_xxpp_format:
            row = xxpp_to_xpxp(row)
        assert n_cols == len(row)
        matrix.append(row)
    if is_xxpp_format:
        matrix = xxpp_to_xpxp(matrix)
    assert n_rows == len(matrix)
    return matrix

def write_file_mat(file_name, mat):
    from MLab import array, Complex, Float, zeros, eig
    # print array(mat, Float)
    res = []
    for r in mat:
        res.append(xpxp_to_xxpp(r))
    mat1 = xpxp_to_xxpp(res)
    #mat1 = xpxp_to_xxpp(self.reverse_order_dof(res, [1, 2, 0]))
    # print array(mat1, Float)
    write_ascii_matrix(file_name, mat1)

def write_ascii_matrix(file_name, mat):
    file_ostr = open(file_name, 'w')
    file_ostr.write('%d\n'%len(mat))
    file_ostr.write('%d\n'%len(mat[1]))
    for row in xrange(len(mat)):
        for col in mat[row]:
            file_ostr.write('%s '%repr(col))
        file_ostr.write('\n')

def read_ascii_vec_polynomials(istream, is_xxpp_format, order=None):
    """

    Read a vector of polynomials from an input stream in ASCII format.
    
    @param istream: an iterable of input lines (e.g. file-like),

    @param is_xxpp_format: a boolean to specify whether the input file
    is in the old (Mathematica) format where all configuration
    coordinates come first, following by all the momenta, rather than
    the new even-odd format.

    @return: vector of polynomials.
    """
    line = istream.next()
    n_vars = int(line)
    assert n_vars >= 0
    line = istream.next()
    n_polys = int(line)
    assert n_polys >= 0
    poly_vec = []
    for i in xrange(n_polys):
        poly_vec.append(read_ascii_polynomial(istream, is_xxpp_format,
                                              order=order))
        assert n_vars == poly_vec[i].n_vars()
    if is_xxpp_format:
        if n_polys == n_vars:
            poly_vec = xxpp_to_xpxp(poly_vec)
    if order != None:
        poly_vec = reorder_dof(poly_vec, order)
    return poly_vec

def write_file_vec_polynomials(file_name, poly_vec, is_xxpp_format):
    file_ostr = open(file_name, 'w')
    n_vars = poly_vec[1].n_vars()
    n_polys = len(poly_vec)
    file_ostr.write('%d\n'%n_vars)
    file_ostr.write('%d\n'%n_polys)
    print file_name, n_polys
    if is_xxpp_format:
        if n_polys == n_vars:
            poly_vec = xpxp_to_xxpp(poly_vec)
    for poly in poly_vec:
        write_ascii_polynomial(file_ostr, poly.real(), is_xxpp_format)
    file_ostr.close()




