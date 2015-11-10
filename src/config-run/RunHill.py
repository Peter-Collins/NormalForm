"""

AUTHOR: Peter Collins, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: RunHill

PURPOSE:

A sample setup and configuration for the normalization algorithms.

NOTES:

See RunConfig.py for configuration options

"""
import sys
import RunConfig
degree = 6
if len(sys.argv)>1:
    degree = int(sys.argv[1])
# pull things into the global context for profile
# from RunConfig import run_nf
# degree 6 runs in about 2m, 8 in 20m, 10 in 2h
config = { "tolerance" : 5.0e-14 , "degree" : degree , "system" : "Hill" ,
           "do_stream" : False ,
           "compute_diagonalisation" : True ,
           "run_normal_form_python" : False ,
           "run_normal_form_cpp" : True }

RunConfig.NfConfig(config).run_examp()

# Now do a python run if degree is < 7
config["compute_diagonalisation"] = False
config["run_normal_form_python"] = True
config["run_normal_form_cpp"] = False
if degree < 7:
    RunConfig.NfConfig(config).run_examp()
