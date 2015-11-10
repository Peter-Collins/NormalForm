"""

AUTHOR: Peter Collins, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: RunCrtbp

PURPOSE:

A sample setup and configuration for the normalization algorithms.

NOTES:

See run_config for configuration options

"""
import RunConfig
# pull things into the global context for profile
from RunConfig import run_nf
config = { "tolerance" : 5.0e-14 , "degree" : 10 , "system" : "Crtbp" }

RunConfig.NfConfig(config).run_examp()
