"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: LieAlgebra

PURPOSE:

Classes to implement graded Lie algebras, based on polynomial rings.

NOTES:

Here, we supply both an abstract interface, and classical and
semi-classical implementations.  The interface insists that derived
classes provide a Lie bracket; for the classical case this is the
Poisson bracket, where as for the semi-classical case it is the Moyal
bracket.

Note that we assume an even/odd mapping of configuration space and
momentum space into a single index space.  In the C++ rewrite of this
code, the information on the form of such an embedding is gathered in
a single place, namely via inlined functions i->iQ(i) and i->iP(i)
and, for the semiclassical Lie algebra, and additional inlined
function iHBar(void).  In the Python code, however, we don't implement
this much flexibility, and enforce the even/odd mapping.

The implementation of the Moyal bracket would be better done
explicitly, rather than in terms of the Moyal product.  For now, we
stick with the implementation in terms of the product, as it has been
well tested.

CREDITS:

The Moyal bracket implementation used here is based on one from Holger
Waalkens, 2005, in terms of the Moyal product.

"""

import logging
logger = logging.getLogger() #('LieAlgebra')

from Numeric import zeros, Complex, matrixmultiply, Float, transpose
from Polynomial import Polynomial
from Powers import Powers
from math import sqrt
from Utility import factorial, binomial
from GradedInterface import GradedInterface
from PolynomialRing import PolynomialRingInterface

class PolynomialLieAlgebraBase(PolynomialRingInterface):
    """

    A symplectic lie algebra structure built on polynomials.

    From PolynomialRingInterface, we inherit:

    .one()

    .zero()

    .coordinate_monomial(i)

    .grad(poly)

    All of which call the .n_vars() method that we must override in
    derived classes; it will depend on .dof(): a classical Lie algebra
    has .vars()->2*.dof() and the semiclassical will have
    .vars()->1+2*.dof().

    """
    def __init__(self, dof):
        self._dof = dof

    #lie algebra:
    def dof(self):
        return self._dof
    def q(self, i):
        return self.coordinate_monomial(2*i)
    def p(self, i):
        return self.coordinate_monomial(2*i+1)
    def bracket(self, a, b):
        raise NotImplementedError
    def is_diagonal_polynomial(self, p):
        """

        In this context, being diagonal means having matching powers on the q
        and p components, i.e. on the configuration and momenta.  Hence
        diagonality belongs in the PolynomialLieAlgebra class.

        """
        self.check_elt(p)
        n_vars = 2*self.dof() #okay for both derived classes
        for m, c in p.powers_and_coefficients():
            found_nonzero = False
            for qi in xrange(0, n_vars, 2):
                pi = qi+1
                if m[qi] != m[pi]:
                    return False
                if m[qi] > 0:
                    if found_nonzero:
                        return False
                    else:
                        found_nonzero = True
        return True
    def diagonal_part_of_polynomial(self, p):
        self.check_elt(p)
        res = self.zero()
        for m, c in p.powers_and_coefficients():
            store = True
            found_nonzero = False
            for qi in xrange(0, 2*self.dof(), 2):
                pi = qi+1
                if m[qi] != m[pi]:
                    store = False
                if m[qi] > 0:
                    if found_nonzero:
                        store = False
                    else:
                        found_nonzero = True
            if store:
                res[m] = c
        return res

class LieAlgebra(PolynomialLieAlgebraBase):
    """

    The classical Lie Algebra structure is implemented as a polynomial
    ring of dimension twice the degrees of freedom, with configuration
    and momentum coordinates alternating.

    The bracket forming the algebra operation is the usual Poisson
    bracket.

    """
    def __init__(self, dof):
        PolynomialLieAlgebraBase.__init__(self, dof)

    #polynomial ring interface:
    def n_vars(self):
        return 2*self.dof()

    #lie algebra interface:
    def bracket(self, a, b):
        return self.poisson_bracket(a, b)
    def poisson_bracket(self, a, b):
        """

        Note that the definition of the Poisson bracket here assumes
        that we are using the even/odd mapping of configuration and
        momentum variables into the polynomial ring.

        """
        self.check_elt(a)
        self.check_elt(b)
        res = self.zero()
        for i in range(0, self.n_vars(), 2):
            p = (a.diff(i)*b.diff(i + 1))
            q = (a.diff(i + 1)*b.diff(i))
            res += (p-q)
        return res

    def display(self, pol_p, prefixes=None):
        if not prefixes:
            prefixes = ['q', 'p']
        assert len(prefixes) == 2
        names = ['%s_{%d}'%(prefixes[i%2], i//2) for i in xrange(self.n_vars())]
        assert len(names) == self.n_vars(), names
        return PolynomialRingInterface.display(self, pol_p, names=names)

class SemiclassicalLieAlgebra(PolynomialLieAlgebraBase):
    """

    The semi-classical Lie Algebra structure is implemented just as
    the classical one, as a polynomial ring, but augmented by an extra
    h-bar variable and thus of one more dimension than its classical
    counterpart.  The mapping between power indices and the
    configuration and momentum components is the same as for the
    classical case (configurations are even, momenta are odd).

    In addition, we have an additional monomial that we can request,
    namely h-bar.

    The bracket here is the Moyal bracket for polynomials.

    """
    def __init__(self, dof):
        PolynomialLieAlgebraBase.__init__(self, dof)

    #polynomial ring interface:
    def n_vars(self):
        return 1+2*self.dof()

    #extra method:
    def grade_of_powers(self, powers):
        """

        The grade of a single powers object in the semiclassical Lie
        Algebra is the degree of the configuration and momentum parts
        plus twice the degree in h-bar.  Here, we implement this by
        taking the degree of the entire powers tuple and then count
        the degree in h-bar one more time.

        """
        assert len(powers) == self.n_vars()
        return powers.degree() + powers[-1]

    def bracket(self, a, b):
        """

        For the semiclassical algebra, the Lie bracket is the Moyal
        bracket.  Here, we implement this in terms of the Moyal
        product.  It would be more efficient to do it directly, since
        the implementation in terms of the product computes many terms
        which in fact cancel each other out.

        """
        return self.moyal_bracket(a, b)

    def moyal_bracket(self, pol_a, pol_b, tolerance=5.0e-14):
        return self.moyal_bracket_via_product(pol_a, pol_b, tolerance)

    def moyal_product(self, pol_a, pol_b, tolerance=1.0e-15):
        """

        This implementation of the Moyal product is adapted from one
        by Dr. Holger Waalkens, 2005.

        """
        self.check_elt(pol_a)
        self.check_elt(pol_b)
        n_max = min(self.grade(pol_a), self.grade(pol_b))
        d = self.dof()
        result = self.zero()
        for n in xrange(0, n_max+1):
            for po_a, co_a in pol_a.powers_and_coefficients():
                for po_b, co_b in pol_b.powers_and_coefficients():
                    m = [0,]*(2*d+1)
                    if (n==0):
                        #usual product:
                        for k in xrange(0, 2*d+1):
                            m[k] = po_a[k] + po_b[k]
                        c = co_a * co_b
                        result += self.monomial(Powers(m), c)
                        #end case n==0:
                    else:
                        for N in xrange(0, (n+1)**(2*d)):
                            Nnew = N
                            for k in xrange(2*d, d, -1):
                                ii = (k - d)-1 #integer div
                                assert 0<=ii<d, ii
                                ip = 2*ii+1
                                assert 0<=ip<2*d
                                m[ip] = int(Nnew) // int((n+1)**(k-1))
                                Nnew -= m[ip]*((n+1)**(k-1))
                            for k in xrange(d, 0, -1):
                                ii = k - 1 #integer div
                                assert 0<=ii<d, ii
                                iq = 2*ii
                                assert 0<=iq<2*d
                                m[iq] = int(Nnew) // int((n+1)**(k-1))
                                Nnew -= m[iq]*((n+1)**(k-1))
                            diff_orderis = 0
                            for k in xrange(0, 2*d):
                                diff_orderis += m[k]
                            if (diff_orderis == n):
                                cc = (0.0+0.5J)**(n)
                                n_mono = [0,]*(2*d+1)
                                n_mono[-1] = n + po_a[-1] + po_b[-1]
                                cc *= (co_a*co_b)
                                for k in xrange(0, 2*d):
                                    cc /= float(factorial(m[k]))
                                sign_exp = 0
                                for k in xrange(1, 2*d, 2): #ip
                                    sign_exp += m[k]
                                cc *= float((-1)**sign_exp)
                                for k in xrange(0, 2*d-1, 2): #iq
                                    iq = k
                                    ip = k+1
                                    pak, qbk = po_a[ip], po_b[iq]
                                    if (m[iq] <= pak) and (m[iq] <= qbk):
                                        n_mono[ip] += pak - m[iq]
                                        n_mono[iq] += qbk - m[iq]
                                        cc *= float(factorial(pak))/float(factorial(pak-m[iq]))
                                        cc *= float(factorial(qbk))/float(factorial(qbk-m[iq]))
                                    else:
                                        cc = 0.0+0.0J
                                    qak, pbk = po_a[iq], po_b[ip]
                                    if (m[ip] <= qak) and (m[ip] <= pbk):
                                        n_mono[iq] += qak - m[ip]
                                        n_mono[ip] += pbk - m[ip]
                                        cc *= float(factorial(qak))/float(factorial(qak-m[ip]))
                                        cc *= float(factorial(pbk))/float(factorial(pbk-m[ip]))
                                    else:
                                        cc = 0.0+0.0J
                                if (abs(cc) > tolerance):
                                    result[Powers(n_mono)] += cc
                            #end diff_order==n
                        #end for N
                    #end case n!=0
                #end pol_b
            #end pol_a
        #end for n
        return result

    def moyal_bracket_via_product(self, pol_a, pol_b, tolerance=5.0e-14):
        """

        This implementation of the Moyal bracket in terms of the product
        involves a division by h-bar.  The implementation of this
        kind of division would be better placed in a method somewhere;
        for now, we leave it here with a clear label.

        """
        ab = self.moyal_product(pol_a, pol_b, tolerance=tolerance)
        ba = self.moyal_product(pol_b, pol_a, tolerance=tolerance)
        hbar_result = (0.0+1.0J)*(ab - ba)
        #now divide by hbar
        result = self.zero()
        for po, co in hbar_result.powers_and_coefficients():
            if abs(co) > tolerance:
                power_hbar = po[-1]
                assert power_hbar > 0, power_hbar
                po_div_hbar = (po.to_tuple())[:-1]+(power_hbar-1,)
                result[Powers(po_div_hbar)] += co
        return result
    
    def grade(self, elt):
        """Overrides GradedInterface."""
        self.check_elt(elt)
        gra_pow = self.grade_of_powers
        res = 0
        for m, c in elt.powers_and_coefficients():
            assert c != 0.0 #security
            gra = gra_pow(m)
            if gra > res:
                res = gra
        return res
    def is_isograde(self, elt, grade=-1):
        """Overrides GradedInterface."""
        self.check_elt(elt)
        gra_pow = self.grade_of_powers
        gra = grade
        for m, c in elt.powers_and_coefficients():
            assert c != 0.0 #security
            if gra == -1:
                gra = gra_pow(m)
            else:
                if gra != gra_pow(m):
                    return False
        return True
    def isograde(self, elt, grade, up_to=None):
        """Overrides GradedInterface."""
        self.check_elt(elt)
        gra_pow = self.grade_of_powers
        if up_to==None:
            up_to = grade+1
        assert grade < up_to
        res = self.zero()
        for m, c in elt.powers_and_coefficients():
            assert c != 0.0 #security
            if grade <= gra_pow(m) < up_to:
                res[m] = c
        return res

    def h_bar(self, pow=1):
        """Semiclassical h-bar variable."""
        return self.coordinate_monomial(-1, pow)

    def display(self, pol_p, prefixes=None):
        """Pretty-print a polynomial."""
        self.check_elt(pol_p)
        if not prefixes:
            prefixes = ['x', 'p']
        assert len(prefixes) == 2
        names = ['%s_{%d}'%(prefixes[i%2], i//2) for i in xrange(2*self.dof())]
        names.append('\\bar{h}')
        assert len(names) == self.n_vars()
        return PolynomialRingInterface.display(self, pol_p, names=names)

    #further generality would be accomplished by, for example:
    #
    #def qi(self, i):
    #    raise NotImplementedError #override in derived classes
    #def pi(self, i):
    #    raise NotImplementedError #override in derived classes
    #
    #def q(self, i):
    #    return self.coordinate_monomial(self.qi(i))
    #def p(self, i):
    #    return self.coordinate_monomial(self.pi(i))
    #
    #def hi(self):
    #    return NotImplementedError
    #def h_bar(self):
    #    return self.coordinate_monomial(self.hi())
    #def grade_of_powers(self, powers):
    #    assert len(powers) == self.n_vars()
    #    return powers.degree() + powers[self.hi()]
    #
    #...etc.

class ClassicalToSemiclassical:
    """

    A convenient class for converting polynomials into semi-classical
    form by augmenting them with an extra h-bar variable.

    """
    def __init__(self, classical_algebra):
        self.cla = classical_algebra
        self.sem = SemiclassicalLieAlgebra(self.cla.dof())
    def classical_algebra(self):
        return self.cla
    def semi_classical_algebra(self):
        return self.sem
    def conversion_substitution(self):
        dst = self.sem
        src = self.cla
        for i in xrange(src.n_vars()):
            yield dst.coordinate_monomial(i)
    def __call__(self, p):
        self.cla.check_elt(p)
        p_sem = p.substitute(list(self.conversion_substitution()))
        self.sem.check_elt(p_sem)
        return p_sem

