"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: NormalForm

PURPOSE:

The various steps in the overall normalization algorithm, collected
together into one class.

NOTES:

The steps are as follows:-

 1. confirmation that we have a valid equilibrium point,

 2. diagonalization of the equilibrium Hamiltonian into real normal
    form,

 3. complexification of the real diagonal system into the complex one,

 4. normalization to complex normal form,

 5. extraction of the complex integrals,

 6. realification of the normal form,

 7. realification of the integrals,

 8. confirmation that the generating function is invariant under the
    coordinate changes,

 9. confirmation that the coordinate changes really do transform the
    diagonal Hamiltonians into their normalized counterparts and vice
    versa.

During this process, the class performs various tests on the validity
of the computations and keeps a log of all results.

"""

import unittest

import logging
logger = logging.getLogger() #('NormalForm')

import LieTriangle
from PolynomialRing import PolynomialRing
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from Powers import Powers
from EigenSystem import EigenSystem
from Diagonal import Diagonalizer, Complexifier
from IntegralExtractor import IntegralExtractor
from CoordinateChange import CoordinateChange
from NormalFormIO import *

from sys import stdout

def truncated_poly(lie_algebra, p, max_grade=None):
    logger.info('Pretty-printing polynomial...')
    if max_grade == None:
        max_grade = lie_algebra.grade(p)
    s = []
    for grade in xrange(0, max_grade+1):
        g = lie_algebra.isograde(p, grade)
        if g:
            s.append('[grade %d]%s'%(grade, lie_algebra.display(g)))
    return '\n'.join(s)

tol_equi_real_coeffs = 1.0e-17
tol_diag_real_coeffs = 1.0e-12
tol_diag_real_imag_part = 1.0e-12
tol_diag_comp_off_diag_part = 5.0e-8
tol_comp_generator_through_lie_triangle = 5.0e-8
tol_norm_real_imag_part = 5.0e-8

class NormalForm:

    def __init__(self, lie_algebra, h_er):
        logger.info('Initializing...')
        self.alg = lie_algebra
        self.h_er = h_er

        logger.info('...done initializing.')

    def get_lie_algebra(self):
        return self.alg

    def set_tolerance(self, tolerance):
        self.tolerance = tolerance

    def confirm_equilibrium_point(self):
        logger.info('Confirming equilibrium point...')

        alg = self.alg

        #logger.info('Real Equilibrium Hamiltonian (terms up to grade')
        #logger.info(truncated_poly(self.alg, self.h_er))
        logger.info('Equilibrium Hamiltonian has %d terms'%self.h_er.n_terms())
        gra = alg.grade(self.h_er)
        if alg.isograde(self.h_er, 0):
            logger.error('WARNING! constant term on input Hamiltonian removed')
            logger.error('size %s', alg.isograde(self.h_er, 0).l_infinity_norm())
            self.h_er = alg.isograde(self.h_er, 1, alg.grade(self.h_er)+1)
        if alg.isograde(self.h_er, 1):
            logger.error('WARNING! linear term on input Hamiltonian removed')
            logger.error('size %s', alg.isograde(self.h_er, 1).l_infinity_norm())
            self.h_er = alg.isograde(self.h_er, 2, alg.grade(self.h_er)+1)
        assert not alg.isograde(self.h_er, 0, 2)
        assert gra == alg.grade(self.h_er)

        logger.info('WARNING! Removing small coefficients of size <=%s'%tol_equi_real_coeffs)
        self.h_er = self.h_er.with_small_coeffs_removed(tol_equi_real_coeffs)
        logger.info('Real Equilibrium Hamiltonian has %d terms'%self.h_er.n_terms())

        logger.info('...done confirming equilibrium point.')

    def diagonalize(self):
        logger.info('Diagonalizing...')
        self.dia = Diagonalizer(self.alg)

        logger.info('Equilibrium Hamiltonian has %d terms'%self.h_er.n_terms())

        logger.info('Computing the eigensystem...')
        self.eig = self.dia.compute_eigen_system(self.h_er, self.tolerance)
        logger.info('...done computing the eigensystem.')

        logger.info('Computing the diagonal change...')
        self.dia.compute_diagonal_change()
        logger.info('...done computing the diagonal change.')

        logger.info('Measuring error in symplectic identities...')
        mat = self.dia.get_matrix_diag_to_equi()
        assert self.dia.matrix_is_symplectic(mat)
        logger.info('...done measuring error in symplectic identities.')

        logger.info('Converting the matrix into vec<poly> representation...')
        er_in_terms_of_dr = self.dia.matrix_as_vector_of_row_polynomials(mat)

        self.er_in_terms_of_dr = er_in_terms_of_dr
        logger.info('... done converting the matrix into vec<poly> repn.')

        logger.info('Applying diagonalizing change to equilibrium H...')
        self.h_dr = self.h_er.substitute(er_in_terms_of_dr)
        logger.info('...done applying diagonalizing change to equilibrium H.')
        del er_in_terms_of_dr

        logger.info('Real diagonal Hamiltonian has %d terms'%self.h_dr.n_terms())

        err = self.h_dr.imag().l_infinity_norm()
        logger.info('Size of imaginary part: %s', err)
        assert err < tol_diag_real_imag_part

        logger.info('WARNING! Removing imaginary terms')
        self.h_dr = self.h_dr.real()
        logger.info('Real diagonal Hamiltonian has %d terms'%self.h_dr.n_terms())

        logger.info('WARNING! Removing small coefficients of size <=%s'%tol_diag_real_coeffs)
        self.h_dr = self.h_dr.with_small_coeffs_removed(tol_diag_real_coeffs)
        logger.info('Real diagonal Hamiltonian has %d terms'%self.h_dr.n_terms())

        #logger.info('Real Diagonal Hamiltonian (terms up to gra 3)')
        #logger.info(truncated_poly(self.alg, self.h_dr))
        logger.info('Quadratic part')
        logger.info(truncated_poly(self.alg, self.alg.isograde(self.h_dr, 2)))
        logger.info('Terms in quadratic part %d'%(self.alg.isograde(self.h_dr, 2).n_terms()))
        
        logger.info('...done diagonalizing.')

    def complexify(self):
        logger.info('Complexifying...')

        self.eq_type = self.eig.get_equilibrium_type()
        logger.info('Equilibrium type %s'%self.eq_type)

        logger.info('Terms in real diagonal %d'%self.h_dr.n_terms())

        self.com = Complexifier(self.alg, self.eq_type)
        self.r_in_terms_of_c = self.com.calc_sub_complex_into_real()
        self.c_in_terms_of_r = self.com.calc_sub_real_into_complex()

        logger.info('Performing complex substitution...')
        self.h_dc = self.h_dr.substitute(self.r_in_terms_of_c)
        logger.info('Terms in complex diagonal %d'%self.h_dc.n_terms())
        
        logger.info('Removing small coefficients...')
        self.h_dc = self.h_dc.with_small_coeffs_removed(self.tolerance)
        logger.info('Terms in complex diagonal %d'%self.h_dc.n_terms())

        h_dc_2 = self.alg.isograde(self.h_dc, 2)
        off_diag = h_dc_2 - self.alg.diagonal_part_of_polynomial(h_dc_2)
        err = off_diag.l_infinity_norm()
        logger.info('Size of off-diagonal part of quadratic part of diagonal h: %s', err)
        assert err < tol_diag_comp_off_diag_part

        logger.info('Removing off-diagonal part')
        self.h_dc -= off_diag
        del off_diag
        logger.info('Terms in complex diagonal %d'%self.h_dc.n_terms())
        
        h_dc_2 = self.alg.isograde(self.h_dc, 2)
        logger.info('Quadratic part')
        logger.info(self.alg.display(h_dc_2))
        
        assert (self.alg.is_diagonal_polynomial(h_dc_2))
        del h_dc_2
        self.eig_vals = self.eig.get_eigen_values_in_plus_minus_pairs()

        logger.info('Complex diagonal Hamiltonian has %d terms'%self.h_dc.n_terms())
        logger.info('...done complexifying.')

    def cvec2r_to_rvec2r(self, p):
        """

        @param p: a polynomial map from a vector of complex
        coordinates to a real scalar.

        @return: a polynomial map from the real version of the
        coordinates to the real scalar.

        """
        return p.substitute(self.c_in_terms_of_r)
        
    def cvec2cvec_to_rvec2rvec(self, vp):
        """

        @param vp: vector of polynomial maps from vector of complex
        coordinates to a complex coordinate.

        @return: vector of polynomial maps from vector of reals to
        real.

        """
        cvec_in_rvec = []
        for i in xrange(len(vp)):
            cvec_in_rvec.append(vp[i].substitute(self.c_in_terms_of_r))
        cvec_in_cvec = []
        for r_in_c_i in self.r_in_terms_of_c:
            cvec_in_cvec.append(r_in_c_i.substitute(cvec_in_rvec))
        return cvec_in_cvec

    def normalize(self):
        logger.info('Normalizing...')
        
        steps = self.desired_grade-1 #1+h_dc.grade()-2

        self.w_c_list = []
        self.h_nc_list = []
        self.h_nc_dict = {}

        self.h_nc = self.alg.zero()
        self.w_c = self.alg.zero()
        h_dc_2 = self.alg.isograde(self.h_dc, 2)
        state = (self.h_nc_dict, self.w_c_list, -1)
        tri = LieTriangle.LieTriangle(self.alg, h_dc_2, state)
            
        #
        # note that, in the following, h_term is _NOT_ h_i for the
        # Deprit triangle; it is (h_i)/factorial(i) and this is
        # taken care of _WITHIN_ the triangle procedure.
        #

        for step in xrange(0, steps + 1):
            h_term = self.alg.isograde(self.h_dc, step + 2)
            k_term, w_term = tri.compute_normal_form_and_generating_function(h_term)
            logger.info('Resulting normal form terms:')
            logger.info(truncated_poly(self.alg, k_term))
            self.h_nc += k_term
            self.w_c += w_term

        logger.info('Complex normal form:')
        logger.info(truncated_poly(self.alg, self.h_nc))

        self.iso_nf = tri.get_isograde_list_handler()

        #h_nc_new = tri.compute_lie_transform_given_generating_function([], steps)
        #logger.info('Normalized using computed generating function:')
        #logger.info(self.alg.display(h_nc_new))

    def check_generating_function_invariant(self):
        logger.info('confirming generating function invariant...')

        grade = self.alg.grade
        #1+self.h_dc.grade()-3 #check
        steps = self.desired_grade-2

        logger.info('surely we can use independent no. steps in coord change,')
        logger.info('than we did in the nf computation? CARE!!!')
        logger.info('the generating function is only invariant with correct steps.')
        
        cc = CoordinateChange(self.alg, self.w_c_list)
        res1 = []
        cc.express_norm_in_diag(self.w_c, res1, {}, steps)
        res2 = []
        cc.express_diag_in_norm(self.w_c, res2, {}, steps)
        iso_cc = cc.get_isograde_list_handler()

        for res in [res1, res2]:
            w_c_new = iso_cc.list_to_poly(res)
            logger.info('Grade %d to %d', grade(self.w_c), grade(w_c_new))
            gra = min((grade(self.w_c), grade(w_c_new)))
            err = self.alg.isograde((w_c_new-self.w_c), 0, gra+1).l_infinity_norm()
            logger.info('Change in like-grade part of generating function: %s', err)
            assert err < tol_comp_generator_through_lie_triangle

        del cc
        del res1
        del res2
        del iso_cc
        del w_c_new

    def check_diag_transforms_to_norm(self):
        logger.info('confirming that diagonal Hamiltonian transforms to normal...')
        grade = self.alg.grade
        
        steps = self.desired_grade-2 #1+self.h_dc.grade()-3
        cc = CoordinateChange(self.alg, self.w_c_list)
        iso_cc = cc.get_isograde_list_handler()
        
        res1 = []
        cc.express_norm_in_diag(self.h_nc, res1, {}, steps)
        h_dc_new = iso_cc.list_to_poly(res1)
        gra = min((grade(self.h_dc), grade(h_dc_new)))
        err = self.alg.isograde(self.h_dc-h_dc_new, 0, gra+1).l_infinity_norm()
        logger.info('Error expressing H_{nc} in x_{dc}')
        logger.info('Grade %d to %d',
                    grade(self.h_dc),
                    grade(h_dc_new))
        logger.info(err)

        logger.info('confirming that normal Hamiltonian transforms to diagonal...')
        
        res2 = []
        cc.express_diag_in_norm(self.h_dc, res2, {}, steps)
        h_nc_new = iso_cc.list_to_poly(res2)
        gra = min((grade(self.h_nc), grade(h_nc_new)))
        err = self.alg.isograde(self.h_nc-h_nc_new, 0, gra+1).l_infinity_norm()
        logger.info('Error expressing H_{dc} in x_{nc}')
        logger.info('Grade %d to %d', grade(self.h_nc), grade(h_nc_new))
        logger.info(err)

        del cc
        del iso_cc
        del h_dc_new
        del h_nc_new
        del res1
        del res2

    def normalize_real(self):
        logger.info('computing real normalization...')
        
        self.h_nr = self.h_nc.substitute(self.c_in_terms_of_r)
        err = self.h_nr.imag().l_infinity_norm()
        tol_norm_real_imag_part = (1.0 + self.h_nr.l_infinity_norm()) * 1.0e-10
        logger.info('l infinity norm of real normal form imaginary part: %s', err)
        assert err < tol_norm_real_imag_part
        self.h_nr = self.h_nr.real()
        logger.info('Real normal form:')
        logger.info(truncated_poly(self.alg, self.h_nr))

        self.w_r_list = [p.substitute(self.c_in_terms_of_r) for p in self.w_c_list]
        #self.w_r = self.iso_nf.list_to_poly(self.w_r_list)

    def extract_integrals(self):
        logger.info('extracting the integrals...')
        
        int_ext = IntegralExtractor(self.alg)
        int_ext.set_complex_normal_hamiltonian(self.h_nc)

        int_ext.find_lost_simple_integrals_and_non_simple_integral()
        logger.info('Lost integrals:')
        logger.info(int_ext._lost_integrals)
        assert ((not int_ext._lost_integrals), int_ext._lost_integrals)
        logger.info('Real simple part:')
        logger.info((self.cvec2r_to_rvec2r(int_ext._simple_part)))
        logger.info('Real non-simple part:')
        logger.info((self.cvec2r_to_rvec2r(int_ext._non_simple_integral)))

        int_ext.list_the_simple_integrals()
        logger.info('Kept integrals:')
        logger.info(int_ext._kept_integrals)

        int_ext.express_simple_part_as_polynomial_over_all_integrals()
        logger.info('Simple part in terms of integrals:')
        logger.info((int_ext._simple_in_integrals))

        int_ext.express_both_parts_as_polynomial_over_all_integrals()
        logger.info('Total in terms of integrals')
        logger.info((int_ext._total_in_integrals))

        int_ext.express_all_integrals_in_normal_form_coords()
        logger.info('Complex Integrals in complex normal form coords:')
        logger.info((int_ext._integrals_in_coords))
        logger.info('Complex Integrals in real normal form coords:')
        ic_in_nr = [self.cvec2r_to_rvec2r(int_i) for int_i in int_ext._integrals_in_coords]
        logger.info(ic_in_nr)

        #convert simple part in terms of complex integrals,
        #to simple part in terms of real integrals.
        def simple_int_is_saddle(real_integral):
            return len(real_integral) == 1
        def simple_int_is_centre(real_integral):
            return len(real_integral) == 2
        comp_ints_in_real_ints = []
        real_integrals = []
        n_integrals = len(int_ext._integrals_in_coords)
        ring_ints = PolynomialRing(n_integrals)
        for i, integral in zip(int_ext._kept_integrals,
                               int_ext._integrals_in_coords):
            # This is the complex integral in real coordinates
            real_integral = integral.substitute(self.c_in_terms_of_r)
            mon = ring_ints.coordinate_monomial(i)
            if simple_int_is_centre(real_integral):
                # Note that a centre real integral = i * complex integral
                comp_ints_in_real_ints.append((-1.0J)*mon)
                real_integrals.append((1.0J)*real_integral)
                logger.info('centre')
            elif simple_int_is_saddle(real_integral):
                comp_ints_in_real_ints.append(mon)
                real_integrals.append(real_integral)
                logger.info('saddle')
            else:
                assert 0, 'Simple integral, neither saddle nor centre!'
        h_ir = int_ext._simple_in_integrals.substitute(comp_ints_in_real_ints)
        logger.info('Simple part in terms of the real simple integrals:')
        logger.info(h_ir)

        logger.info('Real integrals in real nf coords:')
        logger.info(real_integrals)
        
        #we should check this by converting the simple part to real
        h_ir_check = self.cvec2r_to_rvec2r(int_ext._simple_part)
        err = (h_ir_check-h_ir.substitute(real_integrals)).l_infinity_norm()
        logger.info('Error in composite of real integrals and real simple part')
        logger.info(err)
        self.h_ir = h_ir
        self.real_integrals = real_integrals

    def compute_diag_to_norm(self):
        logger.info('computing coordinate changes...')

        steps = self.desired_grade-2 #1+self.h_dc.grade()-3 #check!
        res = []
        for s in xrange(self.alg.n_vars()):
            #do we request a coordinate change in h-bar???
            logger.info('complex coord change, component %d of %d', s, self.alg.n_vars())
            xs = self.alg.coordinate_monomial(s)
            sp = [0,]*(self.alg.n_vars())
            sp[s] = 1

            #use complex nf
            logger.info('COMPLEX component of x_{nc} expressed IN x_{dc} (diag TO norm):')
            x_i_list = []
            x_ij_dict = {}
            cc = CoordinateChange(self.alg, self.w_c_list)
            cc.express_norm_in_diag(xs, x_i_list, x_ij_dict, steps)

            #near-identity?
            assert (len(x_i_list[0]) == 1)
            assert (x_i_list[0][Powers(tuple(sp))] == 1.0)
            for gra_minus_1, pol in enumerate(x_i_list):
                gra = self.alg.grade(pol)
                assert ((gra == 0) or (gra == gra_minus_1+1))

            iso_cc = cc.get_isograde_list_handler()
            #for i, polc in enumerate(x_i_list):
                #logger.info(truncated_poly(self.alg, iso_cc.inner_taylor_to_poly(polc, i)))

            poly = iso_cc.list_to_poly(x_i_list)
            res.append(poly)
        self.nc_in_dc = res
     

    def compute_diag_to_norm_real_using_w_real(self):
        logger.info('computing real coordinate changes using real generator...')

        steps = self.desired_grade-2 #self.h_dc.grade()-2 #check!
        res = []
        for s in xrange(self.alg.n_vars()):
            logger.info('complex coord change, component %d of %d', s, self.alg.n_vars())
            xs = self.alg.coordinate_monomial(s)
            sp = [0,]*(self.alg.n_vars())
            sp[s] = 1

            #use real nf
            logger.info('REAL component of x_{nr} expressed IN x_{dc} (diag TO norm):')
            x_i_list = []
            x_ij_dict = {}
            cc = CoordinateChange(self.alg, self.w_r_list)
            cc.express_norm_in_diag(xs, x_i_list, x_ij_dict, steps)

            #near-identity?
            assert (len(x_i_list[0]) == 1)
            assert (x_i_list[0][Powers(tuple(sp))] == 1.0)
            for gra_minus_1, pol in enumerate(x_i_list):
                gra = self.alg.grade(pol)
                assert ((gra == 0) or (gra == gra_minus_1+1))

            iso_cc = cc.get_isograde_list_handler()
            for i, polc in enumerate(x_i_list):
                err = polc.imag().l_infinity_norm()
                logger.info('l_infinity_norm of imaginary part: %s', err)
                tol_norm_real_imag_part = (1.0 + self.w_r_list[i].l_infinity_norm()) * 1.0e-10
                assert err < tol_norm_real_imag_part
                #logger.info(truncated_poly(self.alg, iso_cc.inner_taylor_to_poly(polc.real(), i)))
            poly = iso_cc.list_to_poly(x_i_list)
            res.append(poly.real())
        self.nr_in_dr_via_wr = res

    def compute_diag_to_norm_real_using_conversion(self):
        logger.info('computing real coordinate changes using real conversion...')

        res = self.cvec2cvec_to_rvec2rvec(self.nc_in_dc)
        self.nr_in_dr_via_nc_in_dc = res
        logger.info('ERRORS in real generator vs. conversion:')
        for m0, m1 in zip(self.nr_in_dr_via_nc_in_dc, self.nr_in_dr_via_wr):
            err = (m1-m0).l_infinity_norm()
            logger.info(err)

    def compute_norm_to_diag(self):
        logger.info('computing inverse coordinate changes...')

        steps = self.desired_grade-2 #1+self.h_dc.grade()-3 #check!
        res = []
        for s in xrange(self.alg.n_vars()):
            logger.info('complex coord change, component %d of %d', s, self.alg.n_vars())
            xs = self.alg.coordinate_monomial(s)
            sp = [0,]*(self.alg.n_vars())
            sp[s] = 1

            #use complex nf
            logger.info('COMPLEX component of x_{dc} expressed IN x_{nc} (norm TO diag):')
            x_i_list = []
            x_ij_dict = {}
            cc = CoordinateChange(self.alg, self.w_c_list)
            cc.express_diag_in_norm(xs, x_i_list, x_ij_dict, steps)

            #near-identity?
            assert (len(x_i_list[0]) == 1)
            assert (x_i_list[0][Powers(tuple(sp))] == 1.0)
            for gra_minus_1, pol in enumerate(x_i_list):
                gra = self.alg.grade(pol)
                assert ((gra == 0) or (gra == gra_minus_1+1))

            iso_cc = cc.get_isograde_list_handler()
            #for i, polc in enumerate(x_i_list):
                #logger.info(truncated_poly(self.alg, iso_cc.inner_taylor_to_poly(polc, i)))

            poly = iso_cc.list_to_poly(x_i_list)
            res.append(poly)
        self.dc_in_nc = res
     
    def compute_norm_to_diag_real_using_w_real(self):
        logger.info('computing real inverse coordinate change using real generator...')
        
        steps = self.desired_grade-2 #self.h_dc.grade()-2 #check!
        res = []
        for s in xrange(self.alg.n_vars()):
            logger.info('complex coord change, component %d of %d', s, self.alg.n_vars())
            xs = self.alg.coordinate_monomial(s)
            sp = [0,]*(self.alg.n_vars())
            sp[s] = 1

            #use real nf
            logger.info('REAL component of x_{dc} expressed IN x_{nr} (norm TO diag):')
            x_i_list = []
            x_ij_dict = {}
            cc = CoordinateChange(self.alg, self.w_r_list)
            cc.express_diag_in_norm(xs, x_i_list, x_ij_dict, steps)

            #near-identity?
            assert (len(x_i_list[0]) == 1)
            assert (x_i_list[0][Powers(tuple(sp))] == 1.0)
            for gra_minus_1, pol in enumerate(x_i_list):
                gra = self.alg.grade(pol)
                assert ((gra == 0) or (gra == gra_minus_1+1))

            iso_cc = cc.get_isograde_list_handler()
            for i, polc in enumerate(x_i_list):
                err = polc.imag().l_infinity_norm()
                logger.info('l_infinity_norm of imaginary part: %s', err)
                tol_norm_real_imag_part = (1.0 + self.w_r_list[i].l_infinity_norm()) * 1.0e-10
                assert err < tol_norm_real_imag_part
            poly = iso_cc.list_to_poly(x_i_list)
            res.append(poly.real())
        self.dr_in_nr_via_wr = res

    def compute_norm_to_diag_real_using_conversion(self):
        logger.info('computing real inverse coordinate changes using real conversion...')

        res = self.cvec2cvec_to_rvec2rvec(self.dc_in_nc)
        self.nr_in_dr_via_nc_in_dc = res
        logger.info('ERRORS in real generator vs. conversion:')
        for m0, m1 in zip(self.nr_in_dr_via_nc_in_dc, self.dr_in_nr_via_wr):
            err = (m1-m0).l_infinity_norm()
            logger.info(err)

    def perform_all_computations(self, desired_grade):
        
        logger.info('main computation started...')

        #flags for extra computations:
        do_check_generating_function_invariant = False
        do_check_diag_transforms_to_norm = False
        do_compute_complex_diag_vs_norm = False
        do_check_diag_vs_norm_using_conversion = False
        
        #setup:
        self.desired_grade = desired_grade
        self.h_er = self.alg.isograde(self.h_er, 0, desired_grade+1)
        self.confirm_equilibrium_point()

        #main steps:
        self.diagonalize()
        self.complexify()
        self.normalize()
        self.normalize_real()
        self.extract_integrals()

        #optional computations:
        if do_check_generating_function_invariant:
            self.check_generating_function_invariant()
        if do_check_diag_transforms_to_norm:
            self.check_diag_transforms_to_norm()
        if do_compute_complex_diag_vs_norm:
            self.compute_diag_to_norm()
            self.compute_norm_to_diag()

        #coordinate change maps:
        self.compute_diag_to_norm_real_using_w_real()
        self.compute_norm_to_diag_real_using_w_real()

        #optional computations:
        if do_check_diag_vs_norm_using_conversion:
            self.compute_diag_to_norm_real_using_conversion()
            self.compute_norm_to_diag_real_using_conversion()

        logger.info('...main computation done.')

        #h_nr_check = self.h_dr.substitute(self.dr_in_nr_via_wr)
        #print "(self.h_dr.substitute(self.dr_in_nr_via_wr)-self.h_nr).l_infinity_norm()\n", h_nr_check.l_infinity_norm(), self.h_nr.l_infinity_norm(), (h_nr_check-self.h_nr).l_infinity_norm(), "\n", h_nr_check.l1_norm(), self.h_nr.l1_norm(), (h_nr_check-self.h_nr).l1_norm()
        #h_dr_check = self.h_nr.substitute(self.nr_in_dr_via_wr)
        #print "(self.h_nr.substitute(self.nr_in_dr_via_wr)-self.h_dr).l_infinity_norm()\n", h_dr_check.l_infinity_norm(), self.h_dr.l_infinity_norm(), (h_dr_check-self.h_dr).l_infinity_norm(), "\n", h_dr_check.l1_norm(), self.h_dr.l1_norm(), (h_dr_check-self.h_dr).l1_norm()
        #print "truncated to degree 5"
        #print self.alg.isograde(h_dr_check-self.h_dr, 1, up_to=5).l1_norm()
        #print self.alg.isograde(h_nr_check-self.h_nr, 1, up_to=5).l1_norm()
        self.write_out_integrals()
        self.write_out_transforms()

    def write_out_transforms(self):
        #write out in Mathematica xxpp format
        is_xxpp_format = True
        write_file_mat("diag_to_equi.mat", self.dia.matrix_diag_to_equi)
        write_file_mat("equi_to_diag.mat",
                           self.dia.get_matrix_equi_to_diag())
        write_file_vec_polynomials("norm_to_diag.vec",
                                       self.dr_in_nr_via_wr, is_xxpp_format)
        write_file_vec_polynomials("diag_to_norm.vec",
                                       self.nr_in_dr_via_wr, is_xxpp_format)
        write_file_vec_polynomials("norm_to_ints.vec",
                                       self.real_integrals, is_xxpp_format)
        write_file_polynomial("equi_to_tham.pol",
                                   self.h_er, is_xxpp_format)

    def write_out_integrals(self):
        #write out in Mathematica xxpp format
        is_xxpp_format = True
        # equi_to_tvec is J.Hessian
        hess = []
        for i in xrange(self.h_er.n_vars()):
            hess.append(self.h_er.diff(i))
        e2tv = [ poly.substitute(hess) for poly in  self.dia.matrix_as_vector_of_row_polynomials(self.dia.skew_symmetric_matrix())]
        write_file_vec_polynomials("equi_to_tvec.vec",e2tv,is_xxpp_format)


        #These are dim dof, is_xxpp_format doesn't apply
        i2f = []
        for i in xrange(self.h_ir.n_vars()):
            i2f.append(self.h_ir.diff(i))
        write_file_vec_polynomials("ints_to_freq.vec", i2f, False)
                       
        file_ostr = open("ints_to_tham.pol", 'w')
        write_ascii_polynomial(file_ostr, self.h_ir.real(), False)
        file_ostr.close()

        logger.info('...done')
