"""

AUTHOR: Dr. Andrew David Burbanks, 2005
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: CoordinateChange

PURPOSE:

Compute the coordinate changes relating complex diagonal to complex
normal form coordinates.

NOTES:

For parallel computation, one can specify for which variable index we
want to compute the coordinate change.

"""

# parallel; we specify the index that we wish to compute.

from Polynomial import Polynomial
from IsogradeInnerTaylorCoeffs import IsogradeInnerTaylorCoeffs
from Utility import binomial, factorial
from LieTriangle import _make_db_key

class CoordinateChange:
   
    """

    This class is responsible for taking scalar functions written in
    one set of variables (either the diagonal complex or normal form
    complex) and expressing them in terms of the other set of
    variables, making use of the generating function that was computed
    during the normalization procedure.  In particular, it can be used
    to express the coordinates themselves (via the scalar function
    projecting onto each coordinate value) in the other coordinate
    system, thereby computing the coordinate-change maps as vectors of
    polynomials.

    """

    def __init__(self, lie_alg, w_i_list):
        self._alg = lie_alg
        self._w_i = w_i_list
        self._iso = IsogradeInnerTaylorCoeffs(self._alg, offset=1) #1

    def get_isograde_list_handler(self):
        return self._iso

    def triangle_diag_in_norm(self, i, j):
        """

        Express scalar function of diagonal coordinates in terms of
        normal form coordinates.

        """
        db_key = _make_db_key(i, j)
        bracket = self._alg.bracket
        if self._x_ij.has_key(db_key):
            pass
        else:
            if (j == 0):
                self._x_ij[db_key] = self._iso.poly_to_inner_taylor(self._f, i)
            else:
                temp = self.triangle_diag_in_norm(i+1, j-1).copy()
                for k in xrange(0, i+1):
                    if (i+1-k) >= len(self._w_i):
                        assert 0, 'out of range, (k+1)=%d'%(k+1)
                    else:
                        w_term = self._w_i[i+1-k]
                    temp += binomial(i, k) * bracket(self.triangle_diag_in_norm(k, j-1), w_term)
                self._x_ij[db_key] = temp
        return self._x_ij[db_key]

    def express_diag_in_norm(self, f, x_i_list, x_ij_dict, n):
        self._f = f
        self._x_i = x_i_list
        self._x_ij = x_ij_dict
        for i in xrange(0, n + 1):
            self._x_i.append(self.triangle_diag_in_norm(0, i))

    def triangle_norm_in_diag(self, i, j):
        """

        Express normal form coordinate expression in terms of diagonal
        coordinates.

        """
        db_key = _make_db_key(i, j)
        bracket = self._alg.bracket
        if self._x_ij.has_key(db_key):
            pass
        else:
            if (i == 0):
                self._x_ij[db_key] = self._iso.poly_to_inner_taylor(self._f, j)
            else:
                temp = self.triangle_norm_in_diag(i-1, j+1).copy()
                for k in xrange(0, i):
                    if (k+1) >= len(self._w_i):
                        assert 0, 'out of range, (k+1)=%d'%(k+1)
                    else:
                        w_term = self._w_i[k+1]
                    temp -= binomial(i-1, k) * bracket(self.triangle_norm_in_diag(i-k-1, j), w_term) #NOTE: MUST BE MINUS!!!
                self._x_ij[db_key] = temp
        return self._x_ij[db_key]

    def express_norm_in_diag(self, f, x_i_list, x_ij_dict, n):
        self._f = f
        self._x_i = x_i_list
        self._x_ij = x_ij_dict
        for i in xrange(0, n + 1):
            self._x_i.append(self.triangle_norm_in_diag(i, 0))
