from thermodynamic_properties.chem_constants import R_si_units
from thermodynamic_properties.critical_constants import CriticalConstants
from thermodynamic_properties.util import percent_difference
from thermodynamic_properties.cp_ig import CpIdealGas
import matplotlib.pyplot as plt
import numpy as np


class Cubic(CriticalConstants, CpIdealGas):
    """Generic Cubic Equation of State
    Defined as in :cite:`Perry`

    .. math::
        P = \\frac{RT}{V - b} - \\frac{a(T)}{(V + \\epsilon b)(V + \\sigma b)}

    where :math:`\\epsilon` and :math:`\\sigma` are pure numbers--the same for all substances.
    :math:`a(T)` and :math:`b` are substance-dependent.

    :param R: gas constant, set to SI units
    :type R: float, hard-coded
    :param sigma: Parameter defined by specific equation of state, :math:`\\sigma`
    :type sigma: float
    :param epsilon: Parameter defined by specific equation of state, :math:`\\epsilon`
    :type epsilon: float
    :param Omega: Parameter defined by specific equation of state, :math:`\\Omega`
    :type Omega: float
    :param Psi: Parameter defined by specific equation of state, :math:`\\Psi`
    :type Psi: float
    :param tol: tolerance for percent difference between compressilibity factor calculated iteratively, set to 0.01
    :type tol: float

    """
    def __init__(self, sigma: float, epsilon: float, Omega: float, Psi: float,
                 dippr_no: str = None, compound_name: str = None, cas_number: str = None, **kwargs):
        """

        :param kwargs: used in :ref:`CpIdealGas`
        """
        CriticalConstants.__init__(self, dippr_no, compound_name, cas_number)
        CpIdealGas.__init__(self, dippr_no, compound_name, cas_number, **kwargs)
        self.R = R_si_units
        self.sigma = sigma
        self.epsilon = epsilon
        self.Omega = Omega
        self.Psi = Psi
        self.b = self.Omega*self.R*self.T_c/self.P_c
        self.tol = 0.01

    def get_default_natural_log(self):
        return np.log

    """Expressions"""
    def alpha_expr(self, T_r):
        """An empirical expression, specific to a particular form of the equation of state

        :param T_r: reduced temperature (T/Tc), dimensionless
        :return: :math:`\\alpha (T_r)`
        """
        raise NotImplementedError

    def beta_expr(self, T, P):
        """

        .. math::
            \\beta = \\Omega\\frac{P_\\mathrm{r}}{T_\\mathrm{r}}

        :param T:  temperature in K
        :param P: pressure in Pa
        :return: :math:`\\beta`
        """
        P_r = P/self.P_c
        T_r = T/self.T_c
        return self.Omega*P_r/T_r

    def q_expr(self, T):
        """

        .. math::
            q = \\frac{\\Psi \\alpha(T_\\mathrm{r})}{\\Omega T_\\mathrm{r}}

        :param T: temperautre in K
        """
        T_r = T / self.T_c
        return self.Psi*self.alpha_expr(T_r)/self.Omega/T_r

    def a_expr(self, T):
        """

        .. math::
            a(T) = \\Psi \\frac{\\alpha(T_r)R^2T_{\\mathrm{c}}^2}{P_{\\mathrm{c}}}

        :param T: temperautre in K
        """
        return self.Psi*self.alpha_expr(T/self.T_c)*self.R*self.R*self.T_c*self.T_c/self.P_c

    def Z_expr(self, P, V, T):
        return P*V/self.R/T

    def I_expr(self, P, V, T, log=None):
        """

        .. math::
            I = \\frac{1}{\\sigma - \\epsilon}\\ln{\\left(
                \\frac{Z + \\sigma \\beta}{Z + \\epsilon \\beta}
            \\right)}

        :param log: function to be used for natural logarithm, defaults to :ref:`np.log`
        :type log: callable, optional
        """
        if log is None:
            log = self.get_default_natural_log()
        Z = self.Z_expr(P, V, T)
        B = self.beta_expr(T, P)
        return log(
            (Z + self.sigma*B)/(Z + self.epsilon*B)
        ) / (self.sigma - self.epsilon)

    def d_ln_alpha_d_ln_Tr(self, T_r):
        """

        :param T_r: reduced temperature [dimensionless]
        :return: Expression for :math:`\\frac{\\mathrm{d} \\ln \\alpha(T_\\mathrm{r})}{\\mathrm{d} \\ln T_\\mathrm{r}}`
        """
        raise NotImplementedError

    def H_R_RT_expr(self, P, V, T, log=None):
        """Dimensionless residual enthalpy

        .. math::
            \\frac{H^\\mathrm{R}}{RT} = Z - 1 + \\left[\\frac{\\mathrm{d} \\ln \\alpha(T_\\mathrm{r})}{\\mathrm{d} \\ln T_\\mathrm{r}} - 1\\right]qI

        :return: Expression for residual enthalpy (divided by RT) -- dimensionless
        """
        return self.Z_expr(P, V, T) - 1 + (self.d_ln_alpha_d_ln_Tr(T/self.T_c) - 1)*self.q_expr(T)*self.I_expr(P, V, T, log=log)

    def G_R_RT_expr(self, P, V, T, log=None):
        """Dimensionless residual gibbs

        .. math::
            \\frac{G^\\mathrm{R}}{RT} = Z - 1 + \\ln(Z-\\beta) - qI

        :return: Expression for residual gibbs free (divided by RT) -- dimensionless
        """
        if log is None:
            log = self.get_default_natural_log()
        Z = self.Z_expr(P, V, T)
        return Z - 1 - log(Z - self.beta_expr(T, P)) - self.q_expr(T)*self.I_expr(P, V, T, log=log)

    def S_R_R_expr(self, P, V, T, log=None):
        """Dimensionless residual entropy

        .. math::
            \\frac{S^\\mathrm{R}}{R} = \\ln(Z-\\beta) + \\frac{\\mathrm{d} \\ln \\alpha(T_\\mathrm{r})}{\\mathrm{d} \\ln T_\\mathrm{r}} qI

        :return: Expression for residual entropy (divided by R) -- dimensionless
        """
        if log is None:
            log = self.get_default_natural_log()
        Z = self.Z_expr(P, V, T)
        return log(Z - self.beta_expr(T, P)) + self.d_ln_alpha_d_ln_Tr(T/self.T_c)*self.q_expr(T)*self.I_expr(P, V, T, log=log)

    def H_expr(self, P, V, T, T_ref, val_ref=0., log=None):
        """Expression for fluid enthalpy

        :param P: pressure in Pa
        :param V: molar volume [m**3/mol]
        :param T: temperautre in K
        :param T_ref: reference temperature in K
        :param val_ref: value at reference temperature [J/mol/K]
        :param log: function used for logarithm
        :return: :math:`H` [J/mol/K]
        """
        return val_ref + self.cp_ig_integral(T_ref, T) + self.R*T*self.H_R_RT_expr(P, V, T, log=log)

    """Solving equations"""
    def coefficients(self, T, P):
        """Polynomial oefficients for cubic equation of state

        .. math::
            Z^3 c_0 + Z^2c_1 + Z*c_2 + c_3 = 0


        :return: :code:`(c_0, c_1, c_2, c_3)`
        """
        B = self.beta_expr(T, P)
        e = self.epsilon
        o = self.sigma
        q = self.q_expr(T)
        return [
            1,
            B*(o + e - 1) - 1,
            B*(B*(e*o - e - o) - e - o + q),
            -B*B*(e*o*(1-B) - q)
        ]

    def cardano_constants(self, T, P):
        """

        :param T: temperature [T]
        :param P: pressure [Pa]
        :return: cardano constants p, q
        :rtype: tuple
        """
        _, b_2, b_1, b_0 = self.coefficients(T, P)
        return (
            b_1 - b_2*b_2/3,                             # p
            b_0 - b_1*b_2/3. + 2.*b_2*b_2*b_2/27.,       # q
            b_2
        )

    def one_root(self, T, P):
        p, q, b_2 = self.cardano_constants(T, P)
        d = (p/3)*(p/3)*(p/3) + (q/2)*(q/2)
        return (-q/2. + d**0.5)**(1./3.) + (-q/2. - d**0.5)**(1./3.) - b_2/3.

    def num_roots(self, T, P):
        """Find number of roots

        See :cite:`Loperena2012,Deiters2002`

        :param T: temperature in K
        :param P: pressure in Pa
        :return: number of roots
        """
        p, q, _ = self.cardano_constants(T, P)
        discriminant = (p/3)*(p/3)*(p/3) + (q/2)*(q/2)
        if discriminant > 0:
            return 1

        return 3

    def check_roots(self, T_min, T_max, P_min, P_max):
        """Check to see if all conditions have one root"""
        for T in np.linspace(T_min, T_max, 10):
            for P in np.linspace(P_min, P_max, 10):
                print('T %5.1f K, P %4.2f MPa has %i roots' % (T, P*1e-6, self.num_roots(T, P)))

    def Z_vapor_RHS(self, Z, beta, q):
        """
        Compressibility of vapor :cite:`Perry`

        .. math::
            1 + \\beta - \\frac{q\\beta (Z - \\beta)}{(Z + \\epsilon\\beta)(Z + \\sigma\\beta)}
            :label: eq_cubic_res_vapor

        :return:
        """
        return 1 + beta - q * beta * (Z - beta) / (Z + self.epsilon * beta) / (Z + self.sigma * beta)

    def Z_liquid_RHS(self, Z, beta, q):
        """
        Compressibility of vapor :cite:`Perry`

        .. math::
            \\beta + (Z + \\epsilon\\beta)(Z + \\sigma\\beta)\\left(\\frac{1 + \\beta - Z}{q\\beta}\\right)
            :label: eq_cubic_res_liquid

        """
        return beta + (Z + self.epsilon * beta) * (Z + self.sigma * beta) * (1 + beta - Z) / q / beta

    def residual(self, P, V, T):
        """

        :param P: pressure in Pa
        :param V: volume in [mol/m**3]
        :param T: temperature in K
        :return: residual for cubic equation of state
        """
        Z = self.Z_expr(P, V, T)
        return Z - self.Z_vapor_RHS(Z, self.beta_expr(T, P), self.q_expr(T))

    def iterate_to_solve_Z(self, T, P, phase):
        """

        :param T: temperature in K
        :param P: pressure in Pa
        :param phase: phase [vapor or liquid]
        :type phase: str
        :return: compressibility factor
        """
        beta = self.beta_expr(T, P)
        q = self.q_expr(T)
        if phase == 'liquid':
            Z_nm1 = beta
            func = self.Z_liquid_RHS
        elif phase == 'vapor':
            Z_nm1 = 1.
            func = self.Z_vapor_RHS
        else:
            raise Exception('Phase {} not found!'.format(phase))
        Z_n = func(Z_nm1, beta, q)
        while percent_difference(Z_n, Z_nm1) > self.tol:
            Z_nm1 = Z_n
            Z_n = func(Z_nm1, beta, q)

        return Z_n

    def iterate_Z_vapor(self, T, P):
        return self.iterate_to_solve_Z(T, P, 'vapor')

    def iterate_Z_liquid(self, T, P):
        return self.iterate_to_solve_Z(T, P, 'liquid')

    def plot_Z_vs_P(self, T, P_min, P_max, phase='vapor', symbol='o', ax=None, **kwargs):
        """Plot compressibility as a function of pressure

        :param T: temperature [K]
        :type T: float
        :param P_min: minimum pressure for plotting [Pa]
        :type P_min: float
        :param P_max: maximum pressure for plotting [Pa]
        :type P_max: float
        :param phase: phase type (liquid or vapor), defaults to vapor
        :type phase: str
        :param symbol: marker symbol, defaults to 'o'
        :type symbol: str
        :param ax: matplotlib axes for plotting, defaults to None
        :type ax: plt.axis
        :param kwargs: keyword arguments for plotting
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_ylabel('$Z$')
            ax.set_xlabel('$P$ [Pa]')

        P = np.linspace(P_min, P_max)
        Z = [self.iterate_to_solve_Z(T, pressure, phase) for pressure in P]
        ax.plot(P, Z, symbol, markerfacecolor='None', **kwargs)


class RedlichKwong(Cubic):
    """
        Redlich-Kwong Equation of State :cite:`RK`
    """

    def __init__(self, **kwargs):
        Cubic.__init__(self, sigma=1, epsilon=0, Omega=0.08664, Psi=0.42748, **kwargs)

    def alpha_expr(self, T_r):
        return 1/T_r**0.5

    def d_ln_alpha_d_ln_Tr(self, T_r):
        return -1./2./T_r/T_r


class SoaveRedlichKwong(RedlichKwong):
    """
        Soave Redlich-Kwong Equation of State :cite:`SRK`

    :param f_w: empirical expression used in :math:`\\alpha` [dimensionless?]
    :type f_w: float, derived
    """

    def __init__(self, **kwargs):
        RedlichKwong.__init__(self, **kwargs)
        self.f_w = 0.480 + 1.574*self.w - 0.176*self.w*self.w

    def alpha_expr(self, T_r):
        val = 1. + self.f_w*(1 - T_r**0.5)
        return val*val

    def d_ln_alpha_d_ln_Tr(self, T_r):
        T_r_sqrt = T_r**0.5
        return (
            self.f_w/(self.f_w*(T_r_sqrt - 1) - 1)/T_r
        )


class PengRobinson(RedlichKwong):
    """
    Peng-Robinson Equation of State :cite:`peng-robinson`

    :param f_w: empirical expression used in :math:`\\alpha` [dimensionless?]
    :type f_w: float, derived
    """

    def __init__(self, **kwargs):
        Cubic.__init__(self, sigma=1 + np.sqrt(2), epsilon=1 - np.sqrt(2), Omega=0.07780, Psi=0.45724, **kwargs)
        self.f_w = 0.37464 + 1.54226*self.w - 0.26992*self.w*self.w
