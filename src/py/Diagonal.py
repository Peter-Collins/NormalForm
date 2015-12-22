"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Diagonal

PURPOSE:

Eigensystems and Diagonalization of the quadratic part of the
Hamiltonian.

NOTES:

Here, we need to look at what happens with zero eigenvalues.

"""

import logging
logger = logging.getLogger() #'Diagonal')

from math import sqrt
from Accuracy import trunc, trunc_tuple
from Complex import real_part_of_vector, imag_part_of_vector
from LieAlgebra import LieAlgebra
from EigenValueVectorPair import EigenValueVectorPair
from EigenSystem import EigenSystem
## Automatically adapted for numpy.oldnumeric Dec 16, 2008 by alter_code1.py
try:
    from numpy.oldnumeric.mlab import array, Complex, Float, zeros
    from numpy.oldnumeric.linear_algebra import determinant, inverse
    from numpy.oldnumeric import transpose, matrixmultiply
    # numpy.oldnumeric.eig seems to return eigvecs as the transpose of the MLab eig
    from numpy.oldnumeric.mlab import eig as eigT
    def eig(linear_matrix):
        eigvals, eigvecs = eigT(linear_matrix)
        eigvecs = transpose(eigvecs)
        return eigvals, eigvecs

except:
    from MLab import array, Complex, Float, zeros, eig
    from LinearAlgebra import determinant, inverse
    from Numeric import transpose, matrixmultiply
from Polynomial import Polynomial

class Diagonalizer:

    def __init__(self, lie_algebra):
        self._lie = lie_algebra
        self._eig = None
        
    def dof(self):
        return self._lie.dof()

    def get_lie_algebra(self):
        return self._lie

    def n_vars(self):
        return self._lie.n_vars()

    def compute_eigen_system(self, h_er, tolerance):
        raw_eig_pairs = self.eigenvalue_eigenvector_pairs(h_er)
        self._eig = EigenSystem(raw_eig_pairs, tolerance=tolerance)
        return self._eig

    def compute_diagonal_change(self):
        """

        For the given eigenvalue/eigenvector pairs of the matrix
        corresponding to the quadratic part of a Hamiltonian (i.e.,
        the linear part of Hamilton's equations at about an
        equilibrium point), compute the corresponding real linear
        symplectic change that brings the quadratic part of the
        Hamiltonian into real diagonal form.

        """
        self._eig.analyse_eigen_system()
        self._compute_symplectic_matrix_diag_to_equi()

    def get_matrix_diag_to_equi(self):
        """

        The matrix which maps a vetor in the real diagonal coordinate
        system into the corresponding vector in the real equilibrium
        coordinate system.  One can take a polynomial expression in
        terms of the real equilibrium coordinates and convert it into
        an expression in the real diagonal coordinates by:-

         1. express the matrix as a vector of linear row-polynomials,
            denoted equi_in_diag, i.e., equilibrium coordinates in
            terms of diagonal ones.

         2. poly_in_diag = poly_in_equi.substitute(equi_in_diag).

        """
        return self.matrix_diag_to_equi

    def get_matrix_equi_to_diag(self):
        """

        The matrix which maps a vetor in the real equilibrium
        coordinate system into the corresponding vector in the real
        diagonal coordinate system.  One can take a polynomial
        expression in terms of the real diagonal coordinates and
        convert it into an expression in the real equilibrium
        coordinates by:-

         1. express the matrix as a vector of linear row-polynomials,
            denoted diag_in_equi, i.e., diagonal coordinates in terms
            of equilibrium ones.

         2. poly_in_equi = poly_in_diag.substitute(diag_in_equi).

        """
        return inverse(self.matrix_diag_to_equi)

    def zero_matrix(self, size=None):
        """

        The $2d\\times 2d$ zero matrix.

        @param size: (optional) size of the zero matrix.

        """
        if size == None:
            d = self._lie.n_vars()
        else:
            d = size
        result = zeros((d, d), Float)
        return result

    def identity_matrix(self, size=None):
        """

        The $2d\\times 2d$ identity matrix.

        @param size: (optional) size of the identity matrix.

        """
        if size == None:
            d = self._lie.n_vars()
        else:
            d = size
        result = zeros((d, d), Float)
        for i in xrange(0, d):
            result[i, i] = +1.0
        return result

    def skew_symmetric_matrix(self):
        """

        The $2d\\times 2d$ skew-symmetric matrix $J$.
        
        For our coordinate ordering (namely, evens are coordinates and
        odds are the conjugate momenta), this is zero except on
        disjoint $2\\times 2$ blocks along the diagonal which carry
        the usual $2\\times 2$ skew-symmetric matrix $J_2$.

        """
        #matrices are accessed by row first, then column:
        d = self._lie.n_vars()
        result = zeros((d, d), Float)
        for qi in xrange(0, d, 2):
            pi = qi+1
            result[qi, pi] = +1.0
            result[pi, qi] = -1.0
        return result

    def vector_l2_norm(self, vec):
        """

        The Euclidean norm of a vector.

        @param vec: an iterable of elements which know abs and pow.

        """
        res = 0.0
        for elt in vec:
            res += abs(elt)**2.0
        return sqrt(res)

    def matrix_norm(self, mat):
        """

        A matrix norm for any shape and size of matrix.

        @param mat: an iter of iters of elements which know abs and
        pow.

        """
        res = 0.0
        for row in mat:
            for elt in row:
                res += abs(elt)**2.0
        return res

    def hamiltons_equations_rhs(self, h):
        """

        The right-hand side of Hamilton's equations for the given
        polynomial Hamiltonian (using terms of all degrees).

        """
        grad_h = self._lie.grad(h)
        result = []
        for qi in xrange(0, self.n_vars(), 2):
            pi = qi+1
            result.append(grad_h[pi])
            result.append(-grad_h[qi])
        return tuple(result)

    def hessian_matrix_of_quadratic_part(self, h):
        """

        The Hessian matrix corresponding to the quadratic part of a
        general polynomial Hamiltonian.

        """
        assert isinstance(h, Polynomial)
        assert h.n_vars()==self._lie.n_vars()
        h_2 = h.homogeneous(2)
        assert h_2.degree()==2 or h_2==self._lie.zero()
        d = self._lie.n_vars()
        result = zeros((d, d), Complex)
        for row in xrange(0, d):
            diff_row = h_2.diff(row)
            for col in xrange(0, d):
                diff_row_col = diff_row.diff(col)
                assert diff_row_col.is_constant()
                if len(diff_row_col) > 0:
                    result[row, col] = diff_row_col.coeffs_as_list()[0]
                else:
                    result[row, col] = 0.0
        return result

    def linear_matrix(self, h):
        """

        Returns the $d\\times d$ matrix that defines the linear part
        of Hamilton's equations for the given polynomial Hamiltonian.

        """
        J = self.skew_symmetric_matrix()
        hess_h_2 = self.hessian_matrix_of_quadratic_part(h)
        return matrixmultiply(J, hess_h_2)

    def eigenvalue_eigenvector_pairs(self, h):
        """

        Compute the eigensystem corresponding to the quadratic part of
        a given Hamiltonian, i.e., to the linear part of Hamilton's
        equations.

        @param h: a polynomial Hamiltonian.
        
        @return: the eigen-(values, vectors) tuple from the quadratic
        part of the Hamiltonian.

        """
        linear_matrix = self.linear_matrix(h)
        eigvals, eigvecs = eig(linear_matrix)
        assert len(eigvals)==self._lie.n_vars()
        assert len(eigvecs)==self._lie.n_vars()
        pairs = []
        for val, vec in zip(eigvals, eigvecs):
            pairs.append(EigenValueVectorPair(val, vec))
        self.check_eigen_value_vector_pairs(linear_matrix, pairs)
        return pairs

    def check_eigen_value_vector_pairs(self, mat, pairs):
        """

        Check the accuracy of eigenvalue/eigenvector pairs for a given
        matrix.

        @param mat: a square matrix.

        @param pairs: an iterable of EigenValueVectorPair.

        """
        for pair in pairs:
            logger.info('Eigenvalue %s', pair.val)
            #logger.info('Eigenvector %s', pair.vec)
            m_v = matrixmultiply(mat, pair.vec)
            e_v = pair.val * pair.vec
            err = self.vector_l2_norm(m_v - e_v)
            logger.info('Error in eigenvalue/eigenvector pair %s', err)
            assert err < 1.0e-12
            if isinstance(pair.val, complex):
                err = min(abs(pair.val.real), abs(pair.val.imag))
            else:
                err = 0.0
            logger.info('Error from pure real/imag %s', err)
            assert err < 1.0e-12
            pair.val = self.purify_eigenvalue(pair.val)

    def purify_eigenvalue(self, val):
        """

        Force a complex eigenvalue to be pure real or pure imaginary.

        """
        if isinstance(val, complex):
            if abs(val.real) > abs(val.imag):
                return complex(val.real, 0.0)
            else:
                return complex(0.0, val.imag)
        else:
            return val

    def _compute_symplectic_matrix_diag_to_equi(self):
        """

        Construct the change of coordinates between the equilibrium
        and real-diagonal coordinates.

        Using the (generally complex) eigenvectors, corresponding to
        the quadratic part of the Hamiltonian, we build a real
        symplectic linear coordinate transformation on the phase space
        (expressed as a matrix), under which the quadratic part is
        transformed into real normal form.

        It is important to note that any non-zero multiple of an
        eigenvector is also an eigenvector.  Therefore, although
        certain eigenvalues (those corresponding to centres) will come
        in pure imaginary complex conjugate pairs, the assocaited
        eigenvectors need not be complex conjugate pairs --- although,
        we could express them as such by multiplication by suitable
        complex scale factors --- instead, we may have any complex
        (including pure-imaginary) multiple of the complex conjugate
        eigenvector.

        Care must therefore be taken in the choice of the eigenvectors
        for the symplectic matrix, especially in the case where the
        original system already contains subsystems that are in real
        normal form.

        To cope with the above, we use the following scheme:

         1. In the case of a saddle (pure real +/- eigenvalues, $(+e,
         -e)$) we simply use whatever eigenvectors, $(v_0, v_1)$, we
         are given and scale (via a real scale factor) these suitably
         to give a real symplectic subsystem.

         2. In the case of a centre (pure imaginary complex conjugate
            pair of eigenvalues), we denote the given eigenvectors by
            $v_0,v_1$.  We then use $v_0':=\re(v_0)$ and
            $v_1':=\im(v_0)$, i.e., we only use one of the complex
            eigenvectors, taking its real and imaginary parts.  Again,
            this choice of real vectors will be suitably scaled so as
            to give a symplectic subsystem.

        It can be confirmed that the action of the linearised version
        of Hamilton's equations on the above choice of real vectors
        spans the same two-dimensional subspace as on the original
        complex ones.

        The chosen vectors form the columns of a matrix $M$.  Given a
        vector in the real diagonal coordinate system, $x_{dr}$, one
        can then compute the equilibrium coordinate system vector from
        $Mx_{dr}=x_{er}$.  So the matrix is the transformation from
        DIAG to EQUI.

        todo: In the case of pairs of zero eigenvalues, perhaps we should
        emit an error, insisting that the user first projects (via a
        momentum-independent projector) to a subspace of the non-null
        space?

        """
        logger.info('Constructing symplectic change:')
        pm_val_vec_pairs = self._eig.plus_minus_pairs
        rearranged_pairs = []
        columns = []
        for i in xrange(0, self.dof()):
            i0 = 2*i
            i1 = i0+1
            plus, minus = pm_val_vec_pairs[i]
            if self._eig.equilibrium_type[i] == 's':
                logger.info('Saddle-type eigenvalues')
                a = real_part_of_vector(plus.vec) #real anyway
                b = real_part_of_vector(minus.vec) #real anyway
            else:
                if self._eig.equilibrium_type[i] == 'c':
                    logger.info('Centre-type eigenvalues')
                    a = real_part_of_vector(plus.vec)
                    b = imag_part_of_vector(plus.vec) #N.B. *same* vector!
                else:
                    logger.info('Zero-type eigenvalues')
                    assert (self._eig.equilibrium_type[i] == '0')
                    #todo: assert imag part is close to zero
                    #todo: assert plus close to minus
                    assert 0, 'do something clever to original hamiltonian, please'
            #logger.info('Original eigenvectors:')
            #logger.info(plus.vec)
            #logger.info(minus.vec)
            #logger.info('Choice of real vectors spanning real subspace:')
            #logger.info(a)
            #logger.info(b)
            c = 0.0
            for k in xrange(0, self.dof()):
                k0 = 2*k
                k1 = k0+1
                c += (a[k0]*b[k1] - a[k1]*b[k0])
            if (abs(c) < 1.0e-12):
                assert (self._eig.equilibrium_type[i] == '0')
                assert 0, 'do something clever to original hamiltonian, please'
            else:
                if c < 0.0:
                    if self._eig.equilibrium_type[i] == 's':
                        a_final = a
                        b_final = -b
                        c_final = -c #ADB!                        
                        logger.info('Negative Multiplier in saddle: %s', c_final)
                        logger.info('Negating second vector')
                    else:
                        #Swapping pairs to make symplectic matrix:
                        rearranged_pairs.append((minus, plus))
                        a_final = b
                        b_final = a
                        c_final = -c #ADB!
                        logger.info('Negative Multiplier in centre: %s', c_final)
                        logger.info('NEGATIVE EIGENVALUE, Swap two vectors')
                else:
                    rearranged_pairs.append((plus, minus))
                    a_final = a
                    b_final = b
                    c_final = +c #ADB!
            logger.info('Multiplier present in eigenvector matrix: %s', c_final)
            scale = +pow(1.0/c_final, 0.5) #+/- #ADB!
            logger.info('Scale factor to produce symplectic matrix: %s', scale)
            columns.append(scale*a_final)
            columns.append(scale*b_final)
        self.final_pairs = rearranged_pairs
        #we must transpose the result, since arrays are created row-wise:
        mat = transpose(array(columns, Float))
        #logger.info('Candidate symplectic matrix')
        #logger.info(mat)
        assert self.matrix_is_symplectic(mat)
        self.matrix_diag_to_equi = mat

    def matrix_is_symplectic(self, m, tolerance=1.0e-12):
        """

        Confirm that a given matrix $M$ is symplectic to within
        numerical tolerance.

        This is done by taking the 4 submatrices:

         1. $A = M[0::2, 0::2]$, i.e., configuration coordinates only,

         2. $B = M[0::2, 1::2]$, i.e., configuration rows, momenta cols,

         3. $C = M[1::2, 0::2]$, i.e., momenta rows, configuration cols,

         4. $D = M[1::2, 1::2]$, i.e., momenta only,

        and verifying the following symplectic identities:-
        
         1. $MJM^{T} = J$, the $2n\\times 2n$ sympletic matrix,
        
         2. $AD^{T}-BC^{T} = I$, the $n\\times n$ identity,
            
         3. $AB^{T}-BA^{T} = Z$, the $n\\times n$ zero,
            
         4. $CD^{T}-DC^{T} = Z$.

        Finally, we confirm that $\\det{M} = 1$.

        """
        det = determinant(m)
        j = self.skew_symmetric_matrix()
        approx_j = matrixmultiply(m, matrixmultiply(j, transpose(m)))
        a = m[0::2, 0::2] #even, even
        b = m[0::2, 1::2] #even, odd
        c = m[1::2, 0::2] #odd, even
        d = m[1::2, 1::2] #odd, odd
        i = self.identity_matrix(self.dof())
        approx_i = matrixmultiply(a, transpose(d)) - matrixmultiply(b, transpose(c))
        approx_z0 = matrixmultiply(a, transpose(b)) - matrixmultiply(b, transpose(a))
        approx_z1 = matrixmultiply(c, transpose(d)) - matrixmultiply(d, transpose(c))
        norm = self.matrix_norm
        logger.info('Matrix from diagonal to equilibrium coordinates:')
        logger.info('[output supressed]') #logger.info(m)
        logger.info('error in determinant:            %s', abs(det-1.0))
        logger.info('error in symplectic identity #1: %s', norm(approx_j - j))
        logger.info('error in symplectic identity #2: %s', norm(approx_i - i))
        logger.info('error in symplectic identity #3: %s', norm(approx_z0))
        logger.info('error in symplectic identity #4: %s', norm(approx_z1))
        okay = True
        if not (abs(det-1.0) < tolerance):
            okay = False
        if not (norm(approx_j - j) < tolerance):
            okay = False
        if not (norm(approx_i - i) < tolerance):
            okay = False
        if not (norm(approx_z0) < tolerance):
            okay = False
        if not (norm(approx_z1) < tolerance):
            okay = False
        return okay

    def matrix_as_vector_of_row_polynomials(self, mat):
        """

        Express a matrix as a tuple of linear polynomials, each
        representing one row.

        In order to express a polynomial in terms of a new coordinate
        system, related to the existing one via a linear change of
        variables (expressed as a matrix), we convert the linear
        transformation to a vector of (linear) polynomials that may be
        substituted into the expression of the original polynomial.

        """
        n_vars = self.n_vars()
        res = []
        for irow in xrange(n_vars):
            lin = Polynomial(n_vars)
            row = mat[irow]
            for icol in xrange(n_vars):
                mon = self._lie.coordinate_monomial(icol)
                lin += row[icol]*mon
            res.append(lin)
        return tuple(res)

class Complexifier:
    def __init__(self, lie_alg, eq_type):
        """

        @param lie_alg: a lie algebra (classical or semiclassical).
        
        @param eq_type: a length $d$ iterable of booleans, each of
        which specifies whether the corresponding degree of freedom
        represents a saddle (True) or a centre (False).

        """
        self.alg = lie_alg
        self.eq_type = eq_type
    def get_lie_algebra(self):
        return self.alg
    def get_eq_type(self):
        return self.eq_type
    def calc_sub_real_into_complex(self):
        """

        Get an expression for complex coordinates in terms of their
        real counterparts.

        This expression should be substituted into an expression in
        terms of the complex coordinates in order to give the
        equivalent expression in terms of real coordinates.

        @return: a tuple of polynomials, each expressing one complex
        coordinate in terms of all of the real coordinates.

        """
        subs = []
        s2 = 1.0/sqrt(2.0)
        for i in xrange(self.alg.dof()):
            if self.eq_type[i] == 'c':
                q = s2 * self.alg.q(i)
                p = s2 * self.alg.p(i)
                # These routines are now both symplectic (det=1)
                # from nf--classical--0.9/py/LieAlgebra.py  
                #the following is correct; see, e.g. Roman's notes.
                #subs.append(q + (1.0J)*p)
                #subs.append((1.0J)*q + p)
                #this is the original form from nf--unified--0.5
                #the following is correct; see, e.g. Roman's notes.
                #subs.append(p - (1.0J)*q)
                #subs.append(q - (1.0J)*p)
                # This agrees with notes by Andy and Holger
                subs.append(q - (1.0J)*p)
                subs.append(p - (1.0J)*q)
            else:
                assert (self.eq_type[i] == 's' or self.eq_type[i] == '0')
                subs.append(self.alg.q(i))
                subs.append(self.alg.p(i))
        if len(subs) == self.alg.n_vars() - 1:
            subs.append(self.alg.coordinate_monomial(-1))
        assert len(subs)==self.alg.n_vars()
        return tuple(subs)
    def calc_sub_complex_into_real(self):
        """

        Get an expression for real coordinates in terms of their
        complex counterparts.

        This expression should be substituted into an expression in
        terms of the real coordinates in order to give the equivalent
        expression in terms of complex coordinates.

        @return: a tuple of polynomials, each expressing one real
        coordinate in terms of all of the complex coordinates.

        """
        subs = []
        s2 = 1.0/sqrt(2.0)
        for i in xrange(self.alg.dof()):
            if self.eq_type[i] == 'c':
                q = s2 * self.alg.q(i)
                p = s2 * self.alg.p(i)
                # This must be the inverse of the previous routine
                # from nf--classical--0.9/py/LieAlgebra.py  
                #subs.append(q + (-1.0J)*p)
                #subs.append((-1.0J)*q + p)
                #this is the original form from nf--unified--0.5
                #subs.append(p + (1.0J)*q)
                #subs.append(q + (1.0J)*p)
                # This agrees with notes by Andy and Holger
                subs.append(q + (1.0J)*p)
                subs.append(p + (1.0J)*q)
            else:
                assert (self.eq_type[i] == 's' or self.eq_type[i] == '0')
                subs.append(self.alg.q(i))
                subs.append(self.alg.p(i))
        if len(subs) == self.alg.n_vars() - 1:
            subs.append(self.alg.coordinate_monomial(-1))
        assert len(subs)==self.alg.n_vars()
        return tuple(subs)

