"""

AUTHOR: Peter Collins, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: RunEckartasMM

PURPOSE:

Example run of the semi-classical normal form procedure, applied to
the EckartasMM system from Holger.

NOTES:

See run_config for configuration options
This is also a test file for a semi-classical test run of both the
python and C++ code, it should run in about 2m27s(c++) + 111m0s(py)
"""
import RunConfig
import NfExample

config = { "tolerance" : 5.0e-14 ,
           "do_stream" : False ,
           "degree" : 4 ,
           "max_degree" : 4 ,
           "compute_diagonalisation" : True ,
           "run_normal_form_python" : False ,
           "run_normal_form_cpp" : True ,
           "system" : "EckartasMM" }

from Powers import Powers
from LieAlgebra import LieAlgebra
from SemiclassicalNormalForm import SemiclassicalNormalForm
from Polynomial import Polynomial
from NormalFormIO import read_ascii_polynomial
import os
import sys
import logging
logger = logging.getLogger()

class EckartasMM(NfExample.NfExample):
    """

    EckartasMM semi-classical system equations for 3-DoF.

    """
    def __init__(self):
        name = "eckartasmm_equi_to_tham.pol"
        logger.info('reading')
        logger.info(name)
        # Get the executable path, not the run dir., and add ../../test
        exedir = sys.path[0]
        in_file = open(os.path.join(exedir, "../../test", name), 'r')
        p = read_ascii_polynomial(in_file,
                                  is_xxpp_format=False)
        in_file.close()
        logger.info('done')
        dof = 3
        self.lie = LieAlgebra(dof)
        self.h_er = p
        grade = 10
        self.prefix='EckartasMM--semi-classical'
        logfile_name = 'EckartasMM--semi-classical--grade-%d.log'%grade
        h_er_trunc = self.lie.isograde(self.h_er, 0, grade+1)
        self.nf = SemiclassicalNormalForm(self.lie, h_er_trunc)
        self.nf.set_tolerance(config["tolerance"])  

def main():
    # Do a python diagonalisation and a C++ run
    syst = EckartasMM()
    RunConfig.NfConfig(config,norm_form=syst.nf).run_syst(syst)

    del syst

    degree = 2
    if len(sys.argv)>1:
        degree = int(sys.argv[1])

    # Now do a python run
    config["compute_diagonalisation"] = False
    config["run_normal_form_python"] = True
    config["run_normal_form_cpp"] = False
    config["degree"] = degree
    syst = EckartasMM()
    RunConfig.NfConfig(config,norm_form=syst.nf).run_syst(syst)


if __name__ == '__main__':
    main()
