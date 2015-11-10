"""

AUTHOR: Peter Collins, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: RunConfig

PURPOSE:

Provide configuration and run control, with various options, for the
normalization algorithms.

NOTES:

These are based on and adapted from a combination of:
projects/nf--unified--0.5/src/py/NormalFormTest.py.
projects/nf--unified--0.5/examples/system-bath/src/SystemBathExample.py

"""

# This config gives the same results as SystemBathExample.py
config = { "tolerance" : 5.0e-12 ,
           "degree" : 4 ,
           "max_degree" : 20 ,
           "n_bath_modes" : 8 ,
           "run_normal_form_cpp" : True ,
           "runprofile" : False ,
           "system" : "SystemBath" }

# We need the software location in the PYTHONPATH maybe from bin/setup.sh
import sys
import os
# Look for a path from the environment and add to $PYTHONPATH
if os.environ.has_key("NFBASEDIR"):
    exedir = os.environ["NFBASEDIR"]
    exedir = os.path.normpath(os.path.join(exedir, "nf--unified--0.5/src/"))
# See if it is in a local user's filespace
elif os.path.exists(os.path.expanduser("~/.NFBASEDIR")):
    f_istr = open(os.path.expanduser("~/.NFBASEDIR"), 'r')
    # read and drop the newline
    exedir = f_istr.readline()[:-1]
    f_istr.close()
    exedir = os.path.normpath(os.path.join(exedir, "nf--unified--0.5/src/"))
# We may be running an example from the svn checkout
else:
    # Get the executable path, not the run dir., and add ../py to $PYTHONPATH
    exedir = sys.path[0]
    exedir = os.path.normpath(os.path.join(exedir, ".."))
print "setting path to", exedir
sys.path.insert(1, os.path.join(exedir, "py"))
exedir = os.path.normpath(os.path.join(exedir, "config-run"))

class NfConfig:
    """
    The config above is used if you run this file:
    python RunConfig.py

    OR in another file do:

    import RunConfig
    config = { "tolerance" : 5.0e-14 , "degree" : 10 , "system" : "Crtbp" }
    RunConfig.NfConfig(config).run()

    Default values are setup below - you can override any of them as above
    
    Make the setup for any one of:
            "Hill" : 			Hill(),			degree = 6
            "RydbergCircular" : 	RydbergCircular(),
            "Hcn" : 			Hcn(),			degree = 10
            "RydbergCrossed" : 		RydbergCrossed(),
            "SingleSaddle" : 		SingleSaddle(),
            "SingleCentre" : 		SingleCentre(),
            "SingleSaddleBad" : 	SingleSaddleBad(),
            "AnotherSaddle" : 		AnotherSaddle(),
            "Rftbp" : 			Rftbp(),
            "Crtbp" : 			Crtbp(),
            "SystemBath" : 		SystemBath(degree)

    """
    def __init__(self, aconfig, norm_form=None ):
        self._aconfig = { "tolerance" : 5.0e-14 ,
                          "degree" : 10 ,
                          "max_degree" : 20 ,
                          "logfile" : 'NormalFormTest.log' ,
                          "do_stream" : True ,
                          "runprofile" : False ,
                          "compute_diagonalisation" : True ,
                          "run_normal_form_python" : False ,
                          "run_normal_form_cpp" : False ,
                          "system" : "SystemBath" }
        self._aconfig.update(aconfig)
        # use a norm_form if given one, this may be SemiClassical
        self._nf = norm_form
        #start logging
        #logfile_name = file_prefix + '.log'
        setup_logging(self._aconfig["logfile"],
                      do_stream=self._aconfig["do_stream"])

    def run_examp(self):
        name = self._aconfig["system"]
        import NfExample
        if "SystemBath" == name:
            examp = NfExample.SystemBath(self._aconfig["n_bath_modes"])
        else:
            examp = eval("NfExample." + self._aconfig["system"] +"()")
        self.run_syst(examp)

    def run_syst(self,ex):
        global example
        example=ex
        if self._aconfig["runprofile"] :
            self.run_profile()
        else:
            self.run_nf()

    def run_nf(self):

        #input normal form computation parameters
        degree = self._aconfig["degree"]
        max_degree = self._aconfig["max_degree"]
        #tolerance is not really unified in the code.
        #tolerance = 5.0e-14
        tolerance = self._aconfig["tolerance"]

        lie = example.get_lie_algebra()
        h_er = example.get_h_er()

        #truncate at requied grade
        h_er_trunc = lie.isograde(h_er, 0, max_degree + 1)
        # Use the, poss SemiClassical, nf given or make one
        if self._nf == None:
            nf = NormalForm(lie, h_er_trunc)
        else:
            nf = self._nf

        #perform diagonalisation, normal form, coordinate changes, etc.
        #compute_diagonalisation
        if self._aconfig["compute_diagonalisation"]:
            diagonalisation_degree = max_degree  # 2 & in line below + 1
            compute_diagonalisation(nf, lie, h_er_trunc,
                                    diagonalisation_degree, tolerance,
                                    example.prefix)
        #run_normal_form_ C plus plus
        if self._aconfig["run_normal_form_cpp"] :
            acommand= os.path.join(exedir, "MakeNormalForm.exe") + \
                      " %s %i %i" % (example.prefix, nf.alg.n_vars(), degree)
            if self._aconfig["do_stream"] == False:
                acommand=acommand+" > Cpp%s 2>&1" % (self._aconfig["logfile"])
            print acommand
            os.system(acommand)

            if isinstance(nf, SemiclassicalNormalForm):
                return

            from IsogradeInnerTaylorCoeffs import IsogradeInnerTaylorCoeffs

            print "using C++ output for integral files"
            #compute polynomial factorial conversions
            _iso = IsogradeInnerTaylorCoeffs(nf.alg, offset=2)

            file_prefix = '%s--k-from-nc' % (example.prefix)
            ring_io = PolynomialRingIO(nf.alg)
            p_list = []
            for grade in xrange(2, degree + 1):
                file_name = '%s--grade-%d.pol' % (file_prefix, grade)
                f_istr = open(file_name, 'r')
                p_grade = ring_io.read_sexp_polynomial(f_istr)
                f_istr.close()
                p_list.append(p_grade)
            nf.h_nc = _iso.list_to_poly(p_list, poly=None)
            nf.extract_integrals()

            file_name = '%s--norm_to_diag.vpol' % (example.prefix)
            f_istr = open(file_name, 'r')
            nf.dr_in_nr_via_wr=ring_io.read_sexp_vector_of_polynomials(f_istr)
            f_istr.close()
            file_name = '%s--diag_to_norm.vpol' % (example.prefix)
            f_istr = open(file_name, 'r')
            nf.nr_in_dr_via_wr=ring_io.read_sexp_vector_of_polynomials(f_istr)
            f_istr.close()
            nf.write_out_integrals()
            nf.write_out_transforms()

        #run_normal_form_python
        if self._aconfig["run_normal_form_python"] :
            compute_normal_form(nf, lie, h_er_trunc, degree, tolerance)

    def run_profile(self):
        using_hotshot = False

        if using_hotshot:
            import hotshot
            import hotshot.stats
        else:
            import profile
            import pstats
        profile_name = 'NormalFormTest.prof'
        proportion_worst_results = 0.125

        create_empty_file(profile_name)
    
        if using_hotshot:
            profiler = hotshot.Profile(profile_name)
            benchtime = profiler.runcall(run_nf) #main call
            profiler.close()
            stats = hotshot.stats.load(profile_name)
        else:
            profile.run('run_nf()', profile_name) #main call
            stats = pstats.Stats(profile_name)

        #output profile
        stats.strip_dirs()
        stats.sort_stats('time', 'cum', 'calls')
        print '[1] Statistics:'
        stats.print_stats(proportion_worst_results)
        print '[2] Callers for the above:'
        stats.print_callers(proportion_worst_results)


def setup_logging(logfile_name, do_stream=False):
    global logger
    import logging
    from sys import stderr
    #clear any existing log file:
    create_empty_file(logfile_name)
    file_handler = logging.FileHandler(logfile_name)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(stderr)
    stream_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    if do_stream:
        logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

#main imports
from Powers import Powers
from LieAlgebra import LieAlgebra
from Polynomial import Polynomial
from NormalFormIO import read_ascii_polynomial
use_generator = 1
if use_generator:
    from NormalForm import NormalForm
else:
    from NormalFormOld import NormalForm
from Diagonal import Complexifier
from SystemBath import SystemBath, new_random_system_bath
from PolynomialRingIO import PolynomialRingIO
from SemiclassicalNormalForm import SemiclassicalNormalForm

def compute_diagonalisation(nf, lie_algebra, h_er_trunc, max_degree, tolerance, system_prefix):
    """

    Here, just go as far as we can with double precision, and dump the
    output files.

    """

    nf.set_tolerance(tolerance)
    desired_grade = max_degree
    nf.desired_grade = max_degree
    logger.info('main computation started...')
    #setup:
    nf.h_er = nf.alg.isograde(nf.h_er, 0, desired_grade+1)
    nf.confirm_equilibrium_point()
    #output max grade
    print "Set up to degree", desired_grade,
    #output h_er
    print "with hamiltonian degree", nf.h_er.degree()
    #diagonalize
    nf.diagonalize()
    #compute the complexification maps
    eq_type = nf.eig.get_equilibrium_type()
    nf.eq_type = eq_type
    com = Complexifier(lie_algebra, eq_type)
    nf.com = com
    r_from_c = com.calc_sub_complex_into_real()
    c_from_r = com.calc_sub_real_into_complex()
    nf.r_in_terms_of_c = r_from_c
    nf.c_in_terms_of_r = c_from_r

    # convert to semi classical and then output files
    if isinstance(nf, SemiclassicalNormalForm):
        nf.quantize()
        r_from_c = nf.r_in_terms_of_c
        c_from_r = nf.c_in_terms_of_r
    #make a polynomial ring io, poss SemiClassical
    lie_algebra = nf.alg
    ring_io = PolynomialRingIO(lie_algebra)
    #dump the Taylor series
    dump_isogrades(nf.alg, nf.h_er, max_degree, example.prefix, 'h-from-e')
    #output vector of polynomials for er in terms of dr
    er_from_dr = nf.er_in_terms_of_dr
    file_name = '%s--e-from-d.vpol' % system_prefix
    file_ostr = open(file_name, 'w')
    ring_io.write_sexp_vector_of_polynomials(file_ostr, er_from_dr)
    file_ostr.close()
    #output equilibrium type
    print nf.eig.get_equilibrium_type()
    #output the complexification maps
    file_name = '%s--r-from-c.vpol' % system_prefix
    file_ostr = open(file_name, 'w')
    ring_io.write_sexp_vector_of_polynomials(file_ostr, r_from_c)
    file_ostr.close()
    file_name = '%s--c-from-r.vpol' % system_prefix
    file_ostr = open(file_name, 'w')
    ring_io.write_sexp_vector_of_polynomials(file_ostr, c_from_r)
    file_ostr.close()

def dump_isogrades(graded_algebra, poly, max_grade, system_prefix,
                   poly_name, min_grade=0, min_grade_number=0):
    logger.info('dumping %s series...' % (poly_name))
    file_prefix = '%s--%s' % (system_prefix, poly_name)
    ring_io = PolynomialRingIO(graded_algebra)
    for grade in xrange(min_grade, max_grade + 1):
        file_name = '%s--grade-%d.pol' % (file_prefix,
                                          grade-min_grade+min_grade_number)
        file_ostr = open(file_name, 'w')
        p_grade = graded_algebra.isograde(poly, grade)
        ring_io.write_sexp_polynomial(file_ostr, p_grade)
        file_ostr.close()
    logger.info('...done')

def compute_normal_form(nf, lie_algebra, h_er_trunc, max_degree, tolerance):
    #perform computations
    nf.set_tolerance(tolerance)
    nf.perform_all_computations(max_degree)
    #dump the Normal Form series
    dump_isogrades(nf.alg, nf.h_nc, max_degree, example.prefix,
                   'k-from-nc-py', min_grade=2, min_grade_number=2)
    #dump the Generating Function series
    dump_isogrades(nf.alg, nf.w_c, max_degree, example.prefix,
                   'w--py', min_grade=2)

def create_empty_file(name):
    file = open(name, 'w')
    file.close()

if __name__ == '__main__':
    thisconfig = NfConfig(config)
    thisconfig.run_examp()

