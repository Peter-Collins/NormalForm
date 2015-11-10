"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: SemiclassicalNormalForm

PURPOSE:

The semi-classical implementation of the normal form algorithm, which
includes quantization of the Hamiltonian.

NOTES:

We inherit from NormalForm and override the perform_all_computations
method to include the quantize/de-quantize steps.

When I say "quantize", here, what I really mean is to take all
polynomial expressions and augment them by moving from the classical
Lie algebra to the corresponding semi-classical Lie algebra, by
addition of the h-bar variable.

"""

from NormalForm import NormalForm
from Diagonal import Complexifier
from LieAlgebra import SemiclassicalLieAlgebra, ClassicalToSemiclassical
import Polynomial

import logging
logger = logging.getLogger() #('SemiclassicalNormalForm')

class SemiclassicalNormalForm(NormalForm):

    def __init__(self, classical_lie_algebra, h_er):
        NormalForm.__init__(self, classical_lie_algebra, h_er)

    def quantize(self):
        cla_to_sem = ClassicalToSemiclassical(self.get_lie_algebra())

        #keep classical versions
        self.alg_classical = self.alg
        self.com_classical = self.com
        self.c_in_terms_of_r_classical = self.c_in_terms_of_r
        self.r_in_terms_of_c_classical = self.r_in_terms_of_c

        # if we have h_dc convert it
        if hasattr(self,"h_dc"):
            self.h_dc_classical = self.h_dc
            self.h_dc = cla_to_sem(self.h_dc)

        #make semiclassical versions
        self.alg = cla_to_sem.semi_classical_algebra()
        n_vars = self.alg.n_vars()
        self.com = Complexifier(self.alg, self.eq_type)
        self.r_in_terms_of_c = self.com.calc_sub_complex_into_real()
        self.c_in_terms_of_r = self.com.calc_sub_real_into_complex()

        #make polys for use by cplusplus
        self.h_er_classical = self.h_er
        self.h_er = cla_to_sem(self.h_er)
        self.er_in_terms_of_dr_classical = self.er_in_terms_of_dr
        er_in_terms_of_dr = []
        for i in xrange(len(self.er_in_terms_of_dr)):
            er_in_terms_of_dr.append(cla_to_sem(self.er_in_terms_of_dr[i]))
        er_in_terms_of_dr.append(Polynomial.Polynomial.CoordinateMonomial(n_vars,n_vars-1))
        self.er_in_terms_of_dr = er_in_terms_of_dr

    def perform_all_computations(self, desired_grade):
        """

        Overrides classical normal form.

        """
        
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
        self.quantize()
        self.normalize()

        #the following need semiclassical real-complex relations...
        
        self.normalize_real()
        #self.extract_integrals()

        #optional computations:
        #if do_check_generating_function_invariant:
        #    self.check_generating_function_invariant()
        #if do_check_diag_transforms_to_norm:
        #    self.check_diag_transforms_to_norm()
        #if do_compute_complex_diag_vs_norm:
        #    self.compute_diag_to_norm()
        #    self.compute_norm_to_diag()

        #coordinate change maps:
        #self.compute_diag_to_norm_real_using_w_real()
        #self.compute_norm_to_diag_real_using_w_real()

        #optional computations:
        #if do_check_diag_vs_norm_using_conversion:
        #    self.compute_diag_to_norm_real_using_conversion()
        #    self.compute_norm_to_diag_real_using_conversion()

        logger.info('...main computation done.')
