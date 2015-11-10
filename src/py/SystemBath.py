"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: SystemBath

PURPOSE:

Used to generate System Bath Hamiltonian, given number of bath modes.

NOTES:

The system defined by the Hamiltonian is (2n+2)-dimensional, and is
ordered (s, p_s, x, p_x, y, p_y, ...).

There is no need for Taylor expansion, as we have the Hamiltonian in
explicit polynomial form.

"""

from math import *
from random import *

from Polynomial import *
from LieAlgebra import LieAlgebra

class SystemBath:
    """

    The original (_not_ mass-weighted) system bath.

    The system-bath model represents a 'system' part: a symmetric
    quartic double-well potential, coupled to a number of 'bath
    modes': harmonic oscillators.  The coupling is achieved via a
    bilinear coupling between the configuration space coordinate of
    the system and the conjugate momenta of each of the bath modes.
    The resulting Hamiltonian is a polynomial of degree 4 in the phase
    space coordinates.

    With this version, the client must specify all of the following:-

    @param n_bath_modes: (non-negative int).
    @param system_mass: (positive real).
    @param imag_harmonic_frequency_at_barrier: (real; imag part of pure imag).
    @param reciprocal_barrier_height_above_well_bottom: (positive real).
    @param bath_masses: (seq of n_bath_modes positive reals).
    @param bath_frequencies: (seq of n_bath_modes reals).
    @param bath_coupling_constants: (seq of n_bath_modes reals).

    """

    def __init__(self,
                 n_bath_modes,
                 system_mass,
                 imag_harmonic_frequency_at_barrier,
                 reciprocal_barrier_height_above_well_bottom,
                 bath_masses,
                 bath_frequencies,
                 bath_coupling_constants):

        assert n_bath_modes>=0
        assert system_mass >= 0.0
        assert abs(imag_harmonic_frequency_at_barrier) > 0.0
        assert reciprocal_barrier_height_above_well_bottom >= 0.0
        assert len(bath_masses) == n_bath_modes
        assert len(bath_frequencies) == n_bath_modes
        assert len(bath_coupling_constants) == n_bath_modes

        for f, g in zip(bath_frequencies[:-1], bath_frequencies[1:]):
            assert f < g

        self._m_s = system_mass #system mass
        self._omega_b = imag_harmonic_frequency_at_barrier
        self._v_0_sh = reciprocal_barrier_height_above_well_bottom
        self._n = n_bath_modes
        self._c = bath_coupling_constants #to the system s coordinate.
        self._w = bath_frequencies
        self._m = bath_masses
        self._lie = LieAlgebra(n_bath_modes+1)

        #$a = (-1/2)m_s\omega_b^2.$
        self._a = -0.5*self._m_s*(self._omega_b**2)
        #$b = \frac{m_s^2\omega_b^4}{16V_0}.$
        self._b = ((self._m_s**2) * (self._omega_b**4))/(16.0*self._v_0_sh)

    def lie_algebra(self):
        """

        Return the Lie algebra on which the polynomials will be
        constructed.  For N bath modes, this has (N+1)-dof.

        """
        return self._lie

    def hamiltonian_real(self):
        """

        Calculate the real Hamiltonian for the system-bath model.

        """
        #Establish some convenient notation:
        n = self._n
        a = self._a
        b = self._b
        m_s = self._m_s
        c = self._c
        w = self._w
        m = self._m
        q_s = self._lie.q(0)
        p_s = self._lie.p(0)

        #Compute some constants:
        coeff_q_s = a
        for i in xrange(0, len(c)):
            coeff_q_s += (c[i]**2.0)/(2.0 * m[i] * (w[i]**2.0))
        coeff_p_s = 1.0/(2.0 * m_s)

        coeff_q_bath = []
        coeff_p_bath = []
        for i in xrange(0, len(c)):
            coeff_q_bath.append(0.5 * m[i] * (w[i]**2.0))
            coeff_p_bath.append(1.0/(2.0 * m[i]))

        #Sanity checks:
        assert n >= 0, 'Need zero or more bath modes.'
        assert len(c) == n, 'Need a coupling constant for each bath mode.'
        assert len(coeff_q_bath) == n, 'Need constant for each bath config.'
        assert len(coeff_p_bath) == n, 'Need constant for each bath momentum.'

        #System part:
        h_system = coeff_p_s * (p_s**2)
        h_system += coeff_q_s * (q_s**2)
        h_system += b * (q_s**4)

        #Bath part:
        h_bath = self._lie.zero()
        for i in xrange(len(c)):
            bath_dof = i+1
            h_bath += coeff_q_bath[i] * (self._lie.q(bath_dof)**2)
            h_bath += coeff_p_bath[i] * (self._lie.p(bath_dof)**2)

        #Coupling part:
        h_coupling = self._lie.zero()
        for i, c_i in enumerate(c):
            bath_dof = i+1
            h_coupling += -c_i * (self._lie.q(bath_dof) * q_s)

        #Complete Hamiltonian:
        h = h_system + h_bath + h_coupling

        #Sanity checks:
        assert h.degree() == 4
        assert h.n_vars() == 2*n+2
        assert len(h) == (3) + (2*n) + (n) #system+bath+coupling
        
        return h

class MassWeightedSystemBath:
    """

    The system-bath model represents a 'system' part (a symmetric
    quartic double-well potential) coupled to a number of 'bath
    modes' (harmonic oscillators).  The coupling is achieved via a
    bilinear coupling between the configuration space coordinate of
    the system and the conjugate momenta of each of the bath modes.
    The resulting Hamiltonian is a polynomial of degree 4 in the phase
    space coordinates.

    """

    def __init__(self,
                 n_bath_modes,
                 imag_harmonic_frequency_at_barrier,
                 reciprocal_barrier_height_above_well_bottom,
                 damping_strength,
                 bath_cutoff_frequency,
                 bath_masses,
                 bath_frequencies,
                 bath_compound_coupling_constants):
        """

        Construct a mass-weighted system bath given the values of the
        parameters and the compound coupling constants.

        @param n_bath_modes: (non-negative int).
        @param imag_harmonic_frequency_at_barrier: (real; im part of pure im).
        @param reciprocal_barrier_height_above_well_bottom: (positive real).
        @param damping_strength: (real).
        @param bath_cutoff_frequency: (real, <<bath_frequencies[-1]).
        @param bath_masses: (seq of n_bath_modes positive reals).
        @param bath_frequencies: (increasing seq of n_bath_modes reals).
        @param bath_compound_coupling_constants: (seq of n_bath_modes reals).

        """

        #check inputs
        assert n_bath_modes>=0
        assert abs(imag_harmonic_frequency_at_barrier) > 0.0
        assert reciprocal_barrier_height_above_well_bottom >= 0.0
        assert len(bath_masses) == n_bath_modes
        assert len(bath_frequencies) == n_bath_modes
        assert len(bath_compound_coupling_constants) == n_bath_modes

        #ensure that the bath frequencies are increasing
        for f, g in zip(bath_frequencies[:-1], bath_frequencies[1:]):
            assert f < g

        #store member variables
        self._n = n_bath_modes
        self._omega_b = imag_harmonic_frequency_at_barrier
        self._v_0_sh = reciprocal_barrier_height_above_well_bottom
        self._eta = damping_strength
        self._omega_c = bath_cutoff_frequency
        self._c_star = bath_compound_coupling_constants #to system coord
        self._w = bath_omegas
        self._lie = LieAlgebra(n_bath_modes+1)

    def compute_compound_constants(n_bath_modes,
                                   damping_strength,
                                   bath_cutoff_frequency,
                                   bath_frequencies):
        """

        Compute the compound coupling constants.

        @param n_bath_modes: (non-negative int).
        @param damping_strength: (real).
        @param bath_cutoff_frequency: (real, <<bath_frequencies[-1]).
        @param bath_frequencies: (seq of n_bath_modes reals increasing).

        @return: bath_compound_coupling_constants (seq of n_bath_modes reals).

        """

        #check inputs
        assert n_bath_modes>=0
        assert len(bath_frequencies) == n_bath_modes
        for f, g in zip(bath_frequencies[:-1], bath_frequencies[1:]):
            assert f < g
        assert bath_frequencies[-1] > bath_cutoff_frequency

        #accumulate compound frequencies
        c_star = []
        omega_c = bath_cutoff_frequency
        eta = damping_strength
        for jm1, omega_j in enumerate(bath_frequencies):
            c = (-2.0/(pi*(jm1+1.0)))*eta*omega_c
            d = ((omega_j+omega_c)*exp(-omega_j/omega_c) - omega_c)
            c_star.append(c*d)
        return c_star
    compute_compound_constants = staticmethod(compute_compound_constants)

    def bath_spectral_density_function(self, omega):
        """

        The bath is defined in terms of a continuous spectral density
        function, which has the Ohmic form with an exponential cutoff.

        For infinite bath cutoff frequency, $\omega_c$, the bath is
        strictly Ohmic, i.e., the friction kernel becomes a delta
        function in the time domain, and the classical dynamics of the
        system coordinate are described by the ordinary Langevin
        equation.  In that case, $eta$ (the damping strength) is the
        classically measurable friction coefficient.

        However, for finite values of the bath cutoff frequency, the
        friction kernel is nonlocal, which introduces memory effects
        into the Generalized Langevin Equation (GLE).

        """
        return self.eta * omega * exp(-omega/self._omega_c)

    def hamiltonian_real(self):
        """

        Calculate the real Hamiltonian for the system-bath model.

        """
        #establish some convenient notation:
        n = self._n
        w = self._w
        c_star = self._c_star

        #sanity checks:
        assert n >= 0, 'Need zero or more bath modes.'
        assert len(c_star) == n, 'Need a coupling constant for each bath mode.'

        #system coefficients:
        a = -0.5*(self._omega_b**2)
        b = (self._omega_b**4)/(16.0*self._v_0_sh)
        coeff_q_s = a
        for i in xrange(0, len(c_star)):
            coeff_q_s += c_star[i]/(2.0 * (w[i]))
        coeff_p_s = 1.0/2.0

        #system part:
        q_s = self._lie.q(0)
        p_s = self._lie.p(0)
        h_system = coeff_p_s * (p_s**2)
        h_system += coeff_q_s * (q_s**2)
        h_system += b * (q_s**4)

        #bath coefficients:
        coeff_q_bath = []
        coeff_p_bath = []
        for i in xrange(0, len(c_star)):
            coeff_q_bath.append(0.5 * (w[i]**2.0))
            coeff_p_bath.append(1.0/2.0)

        #sanity checks:
        assert len(coeff_q_bath) == n, 'Need constant for each bath config.'
        assert len(coeff_p_bath) == n, 'Need constant for each bath momentum.'

        #bath part:
        h_bath = self._lie.zero()
        for i in xrange(len(c_star)):
            bath_dof = i+1
            h_bath += coeff_q_bath[i] * (self._lie.q(bath_dof)**2)
            h_bath += coeff_p_bath[i] * (self._lie.p(bath_dof)**2)

        #coupling part:
        h_coupling = self._lie.zero()
        for i, c_i in enumerate(c_star):
            bath_dof = i+1
            h_coupling += -sqrt(c_i*w[i]) * (self._lie.q(bath_dof)*q_s)

        #complete Hamiltonian:
        h = h_system + h_bath + h_coupling

        #sanity checks:
        assert h.degree() == 4
        assert h.n_vars() == 2*n+2
        assert len(h) == (3) + (2*n) + (n) #system+bath+coupling
        
        return h

def new_random_system_bath(n_bath_modes,
                           system_mass,
                           imag_harmonic_frequency_at_barrier,
                           reciprocal_barrier_height_above_well_bottom,
                           random_seed):
    #check inputs
    assert n_bath_modes >= 0
    assert system_mass >= 0.0

    #unitialize random number generator
    seed(random_seed)

    #generate parameters
    bath_masses = []
    bath_omegas = []
    bath_coupling_constants = []
    for i in xrange(0, n_bath_modes):
        bath_coupling_constants.append(uniform(0.001, 0.5))
        bath_masses.append(uniform(0.5, 3.6))
        bath_omegas.append(gauss(0.0, 2.0))

    #sort frequencies into increasing order
    bath_omegas.sort()

    #instantiate the system bath
    sb = SystemBath(n_bath_modes,
                    system_mass,
                    imag_harmonic_frequency_at_barrier,
                    reciprocal_barrier_height_above_well_bottom,
                    bath_masses,
                    bath_omegas,
                    bath_coupling_constants)
    return sb

