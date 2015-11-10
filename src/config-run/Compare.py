#! python
"""

AUTHOR: Peter Collins, 2006.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Compare.py

PURPOSE:

Read in and compare results calculated by the python and C++ software 

NOTES:

This is a set of utilities for TestNF which ought to go somewhere else.
It's not pretty, don't look.

# you can run this by eg:
# export PYTHONPATH=/home/maprcc/projects/subversion/mainco/NormalForm/nf--unified--0.5/src/py: /home/maprcc/projects/subversion/mainco/NormalForm/nf--unified--0.5/src/config-run
# echo "import compare ; compare.compare_py_cpp()" | python

"""

import glob
import string

#//file_name = '/home/maprcc/projects/nf/nf--unified--0.5/src/config-run/EckartasMM_runs/cpp7/EckartasMM--semi-classical--w--inner-4.pol'

from PolynomialRingIO import PolynomialRingIO
from LieAlgebra import SemiclassicalLieAlgebra

def read_sexp_file(file_name):
    sem = SemiclassicalLieAlgebra(3)
    ring_io = PolynomialRingIO(sem)
    file_istr = open(file_name, 'r')
    return ring_io.read_sexp_polynomial(file_istr)

def compare_poly(poly1, poly2, comment=""):
    poly = poly1 - poly2
    if poly.n_terms() != 0:
        print poly.l_infinity_norm(), '\t', poly.l1_norm() / (poly.n_terms() * 2), \
              '\t', comment

def compare_files(n_vars, file_name1, file_name2, comment="", factor = 1.0):

    from LieAlgebra import LieAlgebra
    from SemiclassicalNormalForm import SemiclassicalNormalForm
    dof = n_vars / 2
    if n_vars % 2:
        lie = SemiclassicalLieAlgebra(dof)
    else:
        lie = LieAlgebra(dof)
    ring_io = PolynomialRingIO(lie)

    file_istr = open(file_name1, 'r')
    poly1=ring_io.read_sexp_polynomial(file_istr)

    file_istr = open(file_name2, 'r')
    poly2=ring_io.read_sexp_polynomial(file_istr)
    compare_poly(factor * poly1, poly2, comment=comment)

def read_n_vars():
    file_name = glob.glob('*r-from-c.vpol')[0]
    file_prefix = file_name.split("r-f")[0]
    file_istr = open(file_name, 'r')
    for line in file_istr.readlines():
        words = string.split(line[1:-2])
        if words[0] == "num-polynomials":
            n_vars = words[1]
        if words[0] == "num-variables":
            assert ( n_vars == words[1] )
            file_istr.close
            return int(n_vars), file_prefix

def compare_py_cpp():
    # works in current dir
    # not files are not the same data, one is inner
    n_vars, file_prefix = read_n_vars() #'hill_l1_18--k-from-nc'
    #num_degrees=len(glob.glob('*w--inner*'))
    num_degrees=len(glob.glob('*k-from-nc-py--*'))
    print "maximum diff            average diff"
    for grade in xrange(2, 2+num_degrees):
        file_name1 = '%sk-from-nc--grade-%d.pol' % (file_prefix, grade)
        file_name2 = '%sk-from-nc-py--grade-%d.pol' % (file_prefix, grade)
        compare_files(n_vars, file_name1, file_name2, comment=file_name1)
    factorial = 1.0
    #file_prefix = 'hill_l1_18--'
    num_degrees=len(glob.glob('*w--py--*'))
    for grade in xrange(0, num_degrees):
        file_name1 = '%sw--inner-%d.pol' % (file_prefix, grade)
        file_name2 = '%sw--py--grade-%d.pol' % (file_prefix, grade)
        # multiply first poly (cpp) by 1/n!
        compare_files(n_vars, file_name1,
                      file_name2, comment=file_name1, factor = 1.0/factorial)
        factorial = factorial * (grade +1)
