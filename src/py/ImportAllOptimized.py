"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: ImportAllOptimized

PURPOSE:

Precompile all Python modules, optionally in optimised compile mode.

NOTES:

"""

import os
import sys
import glob

def importAllModulesInPath_(file_path, do_optimize=0):
    cmd = os.system
    file_names = glob.glob(os.path.join(file_path, '*.py'))
    for file_name in file_names:
        parts = file_name.split('/')
        module_name = parts[-1][:-3]
        print file_name, module_name
        if do_optimize:
            cmd("python2.3 -O -c 'import %s\n'"% module_name)
        cmd("python2.3 -c 'import %s\n'"% module_name)

if __name__ == '__main__':
    do_optimize = 0
    if len(sys.argv)>1:
        if sys.argv[1] == '-O':
            paths = sys.argv[2:]
            do_optimize = 1
            print 'optimized'
        else:
            print 'standard'
            paths = sys.argv[1:]
        for path in paths:
            importAllModulesInPath_(path, do_optimize=do_optimize)

