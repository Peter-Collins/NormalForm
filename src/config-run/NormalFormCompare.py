"""

AUTHOR: Peter Collins, 2006.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: NormalFormCompare

PURPOSE:

Read in and compare Normal Forms calculated by this software with
acsii files from Mathematica when available.

NOTES:
The Mathematica axes ordering and the planes ordering need to be changed.
Use by default the dir
/home/maadb/PACKAGE-TO-HANDOVER/projects/normal-forms/hill/problem/l1/degree_18

"""


import sys
import os
# Get the executable path, not the run dir., and add ../py to $PYTHONPATH
exedir = sys.path[0]
sys.path.insert(1, os.path.join(exedir, "../py"))

from LieAlgebra import LieAlgebra
from NormalForm import NormalForm
from NormalFormIO import read_ascii_polynomial, xxpp_to_xpxp, reorder_dof
from NormalFormIO import read_ascii_matrix, read_ascii_vec_polynomials

from PolynomialRing import PolynomialRing
from PolynomialRingIO import PolynomialRingIO

from MLab import array, Complex, Float, zeros, eig
from LinearAlgebra import determinant, inverse
from Numeric import transpose, matrixmultiply

class NormalFormData:

    def __init__(self, dir_name, order_f=None,
                 is_xxpp_format=False, degree=None):
        self.dir_name = dir_name
        self.read_normal_form_data(self.dir_name,
                                   order_f=order_f,
                                   is_xxpp_format=is_xxpp_format,
                                   degree=degree)

    def __repr__(self):
        body = '%d' % (self._n_vars)
        if self.diag_to_equi:
            body += 'diag_to_equi\n'
            #body += 'diag_to_equi=%s' % repr(self.diag_to_equi)
        return 'NormalFormData:\n(%s)' % body

    def read_normal_form_data(self, dir_name, order_f=None,
                              is_xxpp_format=False, degree=None):
        """
        Read the NF data and reorder the variables xxpp to xpxp and
        the order of the NF planes as required for Mathematica data
        """
        def check_poly_degree(poly, name, degree):
            if degree:
                polyd = poly.degree() 
                if polyd != degree:
                    print "poly %s is degree %d not degree %d" % (name, polyd, degree)

        def check_vec_degree(poly_vec, name, degree):
            for i in xrange(len(poly_vec)):
                check_poly_degree(poly_vec[i], name+"[%d]"%i, degree)

        def read_file_vec_polynomials(
            dir_name, file_name, is_xxpp_format, order=None, degree=degree):
            file_istr = open(os.path.join(dir_name, file_name), 'r')
            vec_poly = read_ascii_vec_polynomials(
                file_istr, is_xxpp_format, order=order)
            file_istr.close()
            check_vec_degree(vec_poly, file_name, degree)
            return vec_poly
        # see what we have
        file_istr = open(os.path.join(dir_name, "norm_to_ints.vec"), 'r')
        line = file_istr.next()
        n_vars = int(line)
        assert n_vars >= 0
        line = file_istr.next()
        n_ints = int(line)
        assert n_ints >= 0
        file_istr.close()
        self._n_vars = n_vars
        assert n_ints == n_vars / 2
        
        dof = n_ints
        #dof = p.n_vars() / 2
        lie = LieAlgebra(dof)

        #nf = NormalForm(lie, h_er)
        # TODO If there's nothing here then do a run to generate data

        #These are dim dof, is_xxpp_format doesn't apply
        file_istr = open(os.path.join(dir_name, "ints_to_freq.vec"), 'r')
        self.ints_to_freq = read_ascii_vec_polynomials(file_istr, False)
        file_istr.close()
        print 'Primary frequencies:', [poly((0.0,)*n_ints)
                                       for poly in self.ints_to_freq]
        my_order_f = [poly((0.0,)*n_ints) for poly in self.ints_to_freq]
        # Don't do reordering for first NF only second with order_f set
        new_order = None
        # Find the order in which the supplied frequencies are found
        if order_f:
            degree = check_poly_degree(self.ints_to_freq[0], "the second normal form", degree)
            degree = 18
            new_order = []
            for i in xrange(n_ints):
                for j in xrange(n_ints):
                    if (abs(my_order_f[j] - order_f[i]) < 1.0e-4):
                        new_order.append(j)
            print "new_order", new_order
            assert len(new_order) == n_ints
            #RE-Read in new order, These are dim dof- no is_xxpp_format
            file_istr = open(os.path.join(dir_name, "ints_to_freq.vec"), 'r')
            self.ints_to_freq = read_ascii_vec_polynomials(file_istr, False,
                                                           order=new_order)
            file_istr.close()
            print 'Primary frequencies:', [poly((0.0,)*n_ints)
                                           for poly in self.ints_to_freq]
        # equi_to_tham - equi coords are not reordered until diagonalised
        file_istr = open(os.path.join(dir_name, "equi_to_tham.pol"), 'r')
        self.equi_to_tham = read_ascii_polynomial(
            file_istr, is_xxpp_format=is_xxpp_format)
        file_istr.close()
        check_poly_degree(self.equi_to_tham, "equi_to_tham.pol", degree)
        self.equi_to_tvec = read_file_vec_polynomials(
            dir_name, "equi_to_tvec.vec",
            is_xxpp_format=is_xxpp_format, degree=degree-1)
        # now reorder to the Mathematica order of the planes
        file_istr = open(os.path.join(dir_name, "ints_to_tham.pol"), 'r')
        self.ints_to_tham = read_ascii_polynomial(file_istr, False,
                                                       order=new_order)
        file_istr.close()
        check_poly_degree(self.ints_to_tham, "ints_to_tham.pol", degree/2)
        # The matrices are different they do the reordering on diag side
        file_istr = open(os.path.join(dir_name, "diag_to_equi.mat"), 'r')
        mat = read_ascii_matrix(file_istr, is_xxpp_format)
        file_istr.close()
        # reorder input variables
        if order_f:
            for i in xrange(len(mat)):
                mat[i] = reorder_dof(mat[i], new_order)
        assert n_vars == len(mat)
        assert n_vars == len(mat[1])
        self.diag_to_equi = mat
        file_istr = open(os.path.join(dir_name, "equi_to_diag.mat"), 'r')
        mat = read_ascii_matrix(file_istr, is_xxpp_format)
        file_istr.close()
        # reorder output variables
        if order_f:
            mat = reorder_dof(mat, new_order)
        assert n_vars == len(mat)
        assert n_vars == len(mat[1])
        self.equi_to_diag = mat
        # reorder all for norm_to_diag, diag_to_norm
        self.norm_to_diag = read_file_vec_polynomials(
            dir_name, "norm_to_diag.vec",
            is_xxpp_format=is_xxpp_format, order=new_order, degree=degree-1)
        self.diag_to_norm = read_file_vec_polynomials(
            dir_name, "diag_to_norm.vec",
            is_xxpp_format=is_xxpp_format, order=new_order, degree=degree-1)
        # norm_to_ints has is_xxpp_format on input vars but not o/p in r_a_v_p
        self.norm_to_ints = read_file_vec_polynomials(
            dir_name, "norm_to_ints.vec",
            is_xxpp_format=is_xxpp_format, order=new_order, degree=2)
    
def poly_vec_substitute(A_in_terms_of_B, B_in_terms_of_C):
    return [poly.substitute(B_in_terms_of_C) for poly in A_in_terms_of_B]

def poly_vec_isograde(poly_vec, grade):
    n_vars = poly_vec[0].n_vars()
    ring = PolynomialRing(n_vars)
    return [ring.isograde(poly, 0, up_to=grade+1) for poly in poly_vec]

def compare_poly(poly1, poly2, comment="", grade=None):
    if grade:
        ring = PolynomialRing(poly1.n_vars())
        poly = ring.isograde(poly1 - poly2, 0, up_to=grade+1)
    else:
        poly = poly1 - poly2
    if poly.n_terms() != 0:
        print poly.l_infinity_norm(), '    \t', \
        poly.l1_norm() / (poly.n_terms() * 2), \
              '    \t', comment
    else:
        print "0 \t\t\t"*2, comment

def compare_poly_vec(poly_v1, poly_v2, comment="", grade=4):
    n_polys = len(poly_v1)
    assert n_polys == len(poly_v2)
    n_vars = poly_v1[0].n_vars()
    assert n_vars == poly_v2[0].n_vars()

    ring = PolynomialRing(n_vars)
    for i in xrange(n_polys):
        poly_vec1 = ring.isograde(poly_v1[i], 0, up_to=grade+1)
        poly_vec2 = ring.isograde(poly_v2[i], 0, up_to=grade+1)
        compare_poly(poly_vec1, poly_vec2, comment+"[%d]" %(i))

def compare_normal_form(dir_name1, dir_name2, grade=4):
    def compare_matrix(dia, mat1, mat2, comment):
        print dia.matrix_norm(array(mat1, Float)-array(mat2, Float)), comment

    print "Comparing data in Normal form files upto grade %d" %(grade)
    print
    print "Reading python data"

    first = NormalFormData(dir_name1, is_xxpp_format=True, degree=grade)

    ring = PolynomialRing(first.equi_to_tham.n_vars())
    
    print
    print "Comparing python data with cpp files upto grade %d" %(grade)
    print "l_infinity_norm \t l1_norm \t\t polynomials"
    ringIO = PolynomialRingIO(ring)
    file_istr = open("hill_l1_18--norm_to_diag.vpol", 'r')
    pv_norm_to_diag = ringIO.read_sexp_vector_of_polynomials(file_istr)
    compare_poly_vec(first.norm_to_diag, 
                     pv_norm_to_diag, "norm_to_diag", grade=grade-1)
    file_istr = open("hill_l1_18--diag_to_norm.vpol", 'r')
    pv_diag_to_norm = ringIO.read_sexp_vector_of_polynomials(file_istr)
    compare_poly_vec(first.diag_to_norm,
                     pv_diag_to_norm, "diag_to_norm", grade=grade-1)
    print
    print "Reading mathematica data"
    n_ints = len(first.ints_to_freq)
    # get the frequencies to find the order of the planes
    order_f=[poly((0.0,)*n_ints) for poly in first.ints_to_freq]
    second = NormalFormData(dir_name2, order_f=order_f,
                            is_xxpp_format=True, degree=grade)
    
    from Diagonal import Diagonalizer
    lie = LieAlgebra(n_ints)
    dia = Diagonalizer(lie)
    grade_ints = grade / 2
    dia.matrix_is_symplectic(array(first.diag_to_equi, Float))
    dia.matrix_is_symplectic(array(second.diag_to_equi, Float))
    dia.matrix_is_symplectic(array(first.equi_to_diag, Float))
    dia.matrix_is_symplectic(array(second.equi_to_diag, Float))

    # For the case develloped, Hill, there is a 45deg rotation between the
    # diagonalised coordinates in each of the centre planes. Thus:
    # These matrices are different
    #compare_matrix(dia, first.diag_to_equi, second.diag_to_equi, "diag_to_equi")
    #compare_matrix(dia, first.equi_to_diag, second.equi_to_diag, "equi_to_diag")
    # We neeed to convert between the diagonal planes and back to
    # compare the nonlinear normalisation plolynomials
    # second.diag_to_first.diag = first.diag_in_terms_of_second.diag =
    fd_in_sd = dia.matrix_as_vector_of_row_polynomials(matrixmultiply(
        array(first.equi_to_diag, Float),array(second.diag_to_equi, Float)))
    sd_in_fd = dia.matrix_as_vector_of_row_polynomials(matrixmultiply(
        array(second.equi_to_diag, Float),array(first.diag_to_equi, Float)))

    print
    print "Comparing mathematica data with cpp files upto grade %d" %(grade-1)
    compare_poly_vec(pv_norm_to_diag,
                     poly_vec_substitute(fd_in_sd, poly_vec_substitute(
        poly_vec_isograde(second.norm_to_diag, grade-1), sd_in_fd)),
                     "norm_to_diag", grade=grade-1)
    print "Comparing mathematica data with cpp files upto grade %d" %(grade-1)
    compare_poly_vec(pv_diag_to_norm,
                     poly_vec_substitute(fd_in_sd, poly_vec_substitute(
        poly_vec_isograde(second.diag_to_norm, grade-1), sd_in_fd)),
                     "diag_to_norm", grade=grade-1)

    print
    print "Comparing mathematica data with python upto grade %d" %(grade)

    compare_poly(first.equi_to_tham, second.equi_to_tham,
                 "equi_to_tham", grade=grade)
    compare_poly_vec(first.ints_to_freq , second.ints_to_freq , "ints_to_freq",
                     grade=grade_ints-1)
    ring_ints = PolynomialRing(second.ints_to_tham.n_vars())
    poly_2 = ring_ints.isograde(second.ints_to_tham, 0, up_to=grade_ints+1)
    compare_poly(first.ints_to_tham , poly_2 , "ints_to_tham")
    compare_poly_vec(first.norm_to_ints , second.norm_to_ints ,
                     "norm_to_ints", grade=grade)

    second.diag_to_norm = poly_vec_isograde(second.diag_to_norm, grade)
    second.norm_to_diag = poly_vec_isograde(second.norm_to_diag, grade)
    compare_poly_vec(first.norm_to_diag, 
                     poly_vec_substitute(fd_in_sd, poly_vec_substitute(
        second.norm_to_diag, sd_in_fd)), "norm_to_diag", grade=grade-1)
    compare_poly_vec(first.diag_to_norm,
                     poly_vec_substitute(fd_in_sd, poly_vec_substitute(
        second.diag_to_norm, sd_in_fd)), "diag_to_norm", grade=grade-1)
    compare_poly_vec(first.diag_to_norm,
                     second.diag_to_norm,"diag_to_norm", grade=grade-1)
    compare_poly_vec(first.diag_to_norm,
                     poly_vec_substitute(fd_in_sd, poly_vec_substitute(
        second.diag_to_norm, sd_in_fd)), "diag_to_norm", grade=grade)
    compare_poly_vec(first.equi_to_tvec,
                     second.equi_to_tvec,"equi_to_tvec", grade=grade-1)

def main():
    if len(sys.argv)>1:
        import os
        degree = 0
        if os.path.isdir(sys.argv[1]):
            firstdir = sys.argv[1]
            if len(sys.argv)==3:
                degree = int(sys.argv[2])
            if len(sys.argv)==4:
                seconddir = sys.argv[2]
                degree = int(sys.argv[3])
        else:
            firstdir = "/home/maprcc/projects/nf/nf--unified--0.5/src/config-run/Hill_runs/nftest2"
            degree = int(sys.argv[1])
    if degree < 1 | degree > 20:
        degree = 6

    # Do a python comparison of normal forms in two directories
    # compare_normal_form(firstdir, "/localhome/maprcc/projects/normal-forms/hill/problem/l1/degree_18/" , grade=degree)
    if len(sys.argv)==4:
        compare_normal_form(firstdir, seconddir , grade=degree)
    else:
        compare_normal_form(firstdir, "../../test/hill_problem_l1_degree_18/" , grade=degree)

if __name__ == '__main__':
    main()
   
