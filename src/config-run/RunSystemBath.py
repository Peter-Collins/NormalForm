"""

AUTHOR: Peter Collins, 2006.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: RunSystemBath

PURPOSE:

A setup and configuration for the SystemBath normalization algorithms.
Read in the degree and the number of bath modes as parameters.

NOTES:

See RunConfig.py for configuration options

"""
import sys
import RunConfig

n_bath_modes = 8
degree = 4
if len(sys.argv)>1:
    degree = int(sys.argv[1])
    if len(sys.argv)>2:
        n_bath_modes = int(sys.argv[2])
                
if len(sys.argv)>3 and sys.argv[3]=="python":
    print "running python only"
    config = { "tolerance" : 5.0e-12 ,
               "degree" : degree ,
               "max_degree" : 20 ,
               "n_bath_modes" : n_bath_modes ,
               "do_stream" : False ,
               "compute_diagonalisation" : False ,
               "run_normal_form_python" : True ,
               "run_normal_form_cpp" : False ,
               "runprofile" : False ,
               "system" : "SystemBath" }
else:
    config = { "tolerance" : 5.0e-12 ,
               "degree" : degree ,
               "max_degree" : 20 ,
               "n_bath_modes" : n_bath_modes ,
               "do_stream" : False ,
               "compute_diagonalisation" : True ,
               "run_normal_form_python" : False ,
               "run_normal_form_cpp" : True ,
               "runprofile" : False ,
               "system" : "SystemBath" }

RunConfig.NfConfig(config).run_examp()
