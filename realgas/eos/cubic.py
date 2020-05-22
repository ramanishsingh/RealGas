r"""

Theory
------

The generic cubic equation of state is :cite:`Perry`

.. math::
    P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}

where :math:`\epsilon` and :math:`\sigma` are pure numbers (the same for all substances),
and :math:`a(T)` and :math:`b` are given by the following equations

.. math::
    a(T) = \Psi \frac{\alpha(T_r)R^2T_{\text{c}}^2}{P_{\text{c}}}
    :label: a_expr

.. math::
    b = \Omega\frac{RT_\text{c}}{P_\text{c}}

The compressibility factor can be calculated by solving the following equation


.. math::
    Z - \left(1 + \beta - \frac{q\beta (Z - \beta)}{(Z + \epsilon\beta)(Z + \sigma\beta)}\right) = 0
    :label: eq_cubic_residual

where

.. math::
    \beta = \Omega\frac{P_\text{r}}{T_\text{r}}
    :label: beta_expr

and

.. math::
    q = \frac{\Psi \alpha(T_\text{r})}{\Omega T_\text{r}}
    :label: q_expr


An iterative routine to calculate :math:`Z` using Equation :eq:`eq_cubic_residual` is implemented,
following :cite:`Perry`, where

.. math::
    Z^\text{new} = 1 + \beta - \frac{q\beta (Z^\text{old} - \beta)}{(Z^\text{old} + \epsilon\beta)(Z^\text{old} + \sigma\beta)}
    :label: eq_iteration

that is continued until the following is true

.. math::
    \|\frac{Z^\text{new} - Z^\text{old}}{Z^\text{new} + Z^\text{old}}\times 200\| < \text{tol}
    :label: tolerance

Residual Molar Properties
-------------------------

.. math::
    \frac{H^\text{R}}{RT} = Z - 1 + \left[\frac{\text{d} \ln \alpha(T_\text{r})}{\text{d} \ln T_\text{r}} - 1\right]qI
    :label: H_R_RT_cubic

.. math::
    \frac{G^\text{R}}{RT} = Z - 1 + \ln(Z-\beta) - qI
    :label: G_R_RT_cubic

.. math::
    \frac{S^\text{R}}{R} = \ln(Z-\beta) + \frac{\text{d} \ln \alpha(T_\text{r})}{\text{d} \ln T_\text{r}} qI
    :label: S_R_R_cubic

where the following functions have been defined

.. math::
    I = \frac{1}{\sigma - \epsilon}\ln{\left(
        \frac{Z + \sigma \beta}{Z + \epsilon \beta}
    \right)}
    :label: I_expr

.. math::
    1 + \beta - \frac{q\beta (Z - \beta)}{(Z + \epsilon\beta)(Z + \sigma\beta)}
    :label: eq_cubic_rhs


Fugacity Coefficients
---------------------
The pure-component fugacity coefficient is *defined* as

.. math::
    \frac{G_i^\text{R}}{RT} = \ln\hat{\phi}_i

So that, from Equation :eq:`G_R_RT_cubic`,

.. math::
    \ln\hat{\phi}_i = Z_i - 1 + \ln(Z_i-\beta_i) - q_iI_i
    :label: ln_hat_phi_i


Mixtures
--------
.. warning::
    Not implemented yet!

.. todo::
    Implement mixtures with cubic equations of state

"""

import matplotlib.pyplot as plt
import numpy as np
from chem_util.chem_constants import gas_constant as R
from chem_util.math import percent_difference

from ..critical_constants import CriticalConstants


class Cubic(CriticalConstants):
    """Generic Cubic Equation of State

    :param sigma: Parameter defined by specific equation of state, :math:`\\sigma`
    :param epsilon: Parameter defined by specific equation of state, :math:`\\epsilon`
    :param Omega: Parameter defined by specific equation of state, :math:`\\Omega`
    :param Psi: Parameter defined by specific equation of state, :math:`\\Psi`
    :param tol: tolerance for iteration (see Equation :eq:`tolerance`), set to 0.01
    :param log: function for computing natural log, defaults to :code:`numpy.log`
    :param exp: function for computing exponential, defaults to :code:`numpy.exp`

    """

    def __init__(self, sigma: float, epsilon: float, Omega: float, Psi: float,
                 dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 log: callable = np.log, exp: callable = np.exp, name: str = 'cubic',
                 **kwargs):
        """

        :param kwargs: for customized critica constants
        """
        CriticalConstants.__init__(self, dippr_no, compound_name, cas_number, **kwargs)
        self.sigma = sigma
        self.epsilon = epsilon
        self.Omega = Omega
        self.Psi = Psi
        self.b = self.Omega * R * self.T_c / self.P_c
        self.tol = 0.01
        self.log = log
        self.exp = exp
        self.name = name

    """Expressions"""

    def alpha_expr(self, T_r):
        """An empirical expression, specific to a particular form of the equation of state

        :param T_r: reduced temperature (T/Tc), dimensionless
        :return: :math:`\\alpha (T_r)` see Table ??
        """
        raise NotImplementedError

    def beta_expr(self, T, P):
        """

        :param T:  temperature in K
        :param P: pressure in Pa
        :return: :math:`\\beta` (see Equation :eq:`beta_expr`)
        """
        P_r = P / self.P_c
        T_r = T / self.T_c
        return self.Omega * P_r / T_r

    def q_expr(self, T):
        """


        :param T: temperautre in K
        :returns: :math:`q` (see Equation :eq:`q_expr`)
        """
        T_r = T / self.T_c
        return self.Psi * self.alpha_expr(T_r) / self.Omega / T_r

    def a_expr(self, T):
        """

        :param T: temperautre in K
        :returns: :math:`a(T)` (see Equation :eq:`a_expr`)
        """
        return self.Psi * self.alpha_expr(T / self.T_c) * R * R * self.T_c * self.T_c / self.P_c

    def Z_expr(self, P, V, T):
        return P * V / R / T

    def I_expr(self, P, V, T):
        """

        :param P: pressure in Pa
        :param V: molar volume in m**3/mol
        :param T: temperature in K
        :returns: :math:`I` (see Equation :eq:`I_expr`)
        """
        Z = self.Z_expr(P, V, T)
        B = self.beta_expr(T, P)
        return self.log(
            (Z + self.sigma * B) / (Z + self.epsilon * B)
        ) / (self.sigma - self.epsilon)

    def d_ln_alpha_d_ln_Tr(self, T_r):
        """

        :param T_r: reduced temperature [dimensionless]
        :return: Expression for :math:`\\frac{\\mathrm{d} \\ln \\alpha(T_\\mathrm{r})}{\\mathrm{d} \\ln T_\\mathrm{r}}`
        """
        raise NotImplementedError

    def H_R_RT_expr(self, P, V, T):
        r"""Dimensionless residual enthalpy


        :param P: pressure in Pa
        :param V: molar volume in m**3/mol
        :param T: temperature in K
        :returns: :math:`\frac{H^\text{R}}{RT}`, see Equation :eq:`H_R_RT_cubic`
        """
        return self.Z_expr(P, V, T) - 1 \
               + (self.d_ln_alpha_d_ln_Tr(T / self.T_c) - 1) * self.q_expr(T) * self.I_expr(P, V, T)

    def G_R_RT_expr(self, P, V, T):
        r"""Dimensionless residual gibbs

        :param P: pressure in Pa
        :param V: molar volume in m**3/mol
        :param T: temperature in K
        :returns: :math:`\frac{G^\text{R}}{RT}`, see Equation :eq:`G_R_RT_cubic`
        """
        Z = self.Z_expr(P, V, T)
        return Z - 1 - self.log(Z - self.beta_expr(T, P)) - self.q_expr(T) * self.I_expr(P, V, T)

    def S_R_R_expr(self, P, V, T):
        r"""Dimensionless residual entropy

        :param P: pressure in Pa
        :param V: molar volume in m**3/mol
        :param T: temperature in K
        :returns: :math:`\frac{S^\text{R}}{R}`, see Equation :eq:`S_R_R_cubic`
        """
        Z = self.Z_expr(P, V, T)
        return self.log(Z - self.beta_expr(T, P)) \
               + self.d_ln_alpha_d_ln_Tr(T / self.T_c) * self.q_expr(T) * self.I_expr(P, V, T)

    def ln_hat_phi_k_expr(self, P, V, T):
        r"""

        :param P: pressure in Pa
        :param T: temperature in K
        :param V: molar volume in m**3/mol
        :returns: :math:`\ln{\hat{\phi}_k`, see Equation :math:`ln_hat_phi_i`
        """
        return self.G_R_RT_expr(P, V, T)

    def hat_phi_i_expr(self, *args):
        r"""expression for fugacity coefficient
        :returns: :math:`\exp\left({\ln{\hat{\phi}_i}}\right)`
        """
        return self.exp(self.ln_hat_phi_k_expr(*args))

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
            B * (o + e - 1) - 1,
            B * (B * (e * o - e - o) - e - o + q),
            -B * B * (e * o * (1 - B) - q)
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
            b_1 - b_2 * b_2 / 3,  # p
            b_0 - b_1 * b_2 / 3. + 2. * b_2 * b_2 * b_2 / 27.,  # q
            b_2
        )

    def one_root(self, T, P):
        p, q, b_2 = self.cardano_constants(T, P)
        d = (p / 3) * (p / 3) * (p / 3) + (q / 2) * (q / 2)
        return (-q / 2. + d ** 0.5) ** (1. / 3.) + (-q / 2. - d ** 0.5) ** (1. / 3.) - b_2 / 3.

    def num_roots(self, T, P):
        """Find number of roots

        See :cite:`Loperena2012,Deiters2002`

        :param T: temperature in K
        :param P: pressure in Pa
        :return: number of roots
        """
        p, q, _ = self.cardano_constants(T, P)
        discriminant = (p / 3) * (p / 3) * (p / 3) + (q / 2) * (q / 2)
        if discriminant > 0:
            return 1

        return 3

    def print_roots(self, T, P):
        """Check to see if all conditions have one root"""
        print('%s with %s EOS at T %5.1f K, P %4.2f MPa has %i roots' % (
            self.compound_name, self.name, T, P * 1e-6, self.num_roots(T, P)))

    def Z_right_hand_side(self, Z, beta, q):
        """
        Estimate of compressibility of vapor :cite:`Perry`, used for iterative methods

        :return: RHS of residual, see Equation :eq:`eq_cubic_rhs`
        """
        return 1 + beta - q * beta * (Z - beta) / (Z + self.epsilon * beta) / (Z + self.sigma * beta)

    def residual(self, P, V, T):
        """

        :param P: pressure in Pa
        :param V: volume in [mol/m**3]
        :param T: temperature in K
        :return: residual for cubic equation of state (Equation :eq:`eq_cubic_residual`)
        """
        Z = self.Z_expr(P, V, T)
        return Z - self.Z_right_hand_side(Z, self.beta_expr(T, P), self.q_expr(T))

    def iterate_to_solve_Z(self, T, P) -> float:
        """
        Iterate to compute :math:`Z` using Equation :eq:`eq_iteration`
        util termination condition (Equation :eq:`tolerance` is met)

        :param T: temperature in K
        :param P: pressure in Pa
        :return: compressibility factor
        """
        beta = self.beta_expr(T, P)
        q = self.q_expr(T)
        Z_nm1 = 1.  # initial guess
        Z_n = self.Z_right_hand_side(Z_nm1, beta, q)
        while abs(percent_difference(Z_n, Z_nm1)) > self.tol:
            Z_nm1 = Z_n
            Z_n = self.Z_right_hand_side(Z_nm1, beta, q)

        return Z_n

    def plot_Z_vs_P(self, T: float, P_min: float, P_max: float, symbol='o', ax: plt.axis = None, fig: plt.figure = None,
                    **kwargs):
        """Plot compressibility as a function of pressure

        :param T: temperature [K]
        :param P_min: minimum pressure for plotting [Pa]
        :param P_max: maximum pressure for plotting [Pa]
        :param symbol: marker symbol, defaults to 'o'
        :param ax: matplotlib axes for plotting, defaults to None
        :param kwargs: keyword arguments for plotting
        """
        if ax is None:
            if fig is None:
                fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_ylabel('$Z$')
            ax.set_xlabel('$P$ [Pa]')

        P = np.linspace(P_min, P_max)
        Z = [self.iterate_to_solve_Z(T, pressure) for pressure in P]
        ax.plot(P, Z, symbol, markerfacecolor='None', **kwargs)
        return fig, ax


class RedlichKwong(Cubic):
    r"""Redlich-Kwong Equation of State :cite:`RK`

    This Equation of state has the following parameters :cite:`Smith2005`

    ===========================       ==============================
            Symbol                                  Value
    ===========================       ==============================
    :math:`\alpha(T_\text{r})`          :math:`1/\sqrt{T_\text{r}}`
    :math:`\sigma`                               1
    :math:`\epsilon`                             0
    :math:`\Omega`                                  0.08664
    :math:`\Psi`                                    0.42748
    ===========================       ==============================


    >>> from realgas.eos.cubic import RedlichKwong
    >>> model = RedlichKwong(compound_name='Propane')
    >>> model.sigma
    1
    >>> model.epsilon
    0
    >>> model.Omega
    0.08664
    >>> model.Psi
    0.42748

    """

    def __init__(self, **kwargs):
        Cubic.__init__(self, sigma=1, epsilon=0, Omega=0.08664, Psi=0.42748, name='RK', **kwargs)

    def alpha_expr(self, T_r):
        return 1 / T_r ** 0.5

    def d_ln_alpha_d_ln_Tr(self, T_r):
        return -1. / 2. / T_r / T_r


class SoaveRedlichKwong(RedlichKwong):
    r"""Soave Redlich-Kwong Equation of State :cite:`SRK`

    This equation of state has the following parameters

    ===========================       ====================================================
            Symbol                                  Value
    ===========================       ====================================================
    :math:`\alpha(T_\text{r})`         :math:`\left[1 + f_w(1-\sqrt{T_\text{r}})\right]^2`
    :math:`\sigma`                               :math:`1`
    :math:`\epsilon`                             :math:`0`
    :math:`\Omega`                                :math:`0.08664`
    :math:`\Psi`                                    :math:`0.42748`
    ===========================       ====================================================

    where

    .. math::
        f_w = 0.480 + 1.574\omega - 0.176\omega^2
        :label: f_w_SRK

    >>> from realgas.eos.cubic import SoaveRedlichKwong
    >>> model = SoaveRedlichKwong(compound_name='Water')
    >>> model.Omega
    0.08664
    >>> model.sigma
    1
    >>> model.epsilon
    0
    >>> model.Psi
    0.42748
    >>> model.f_w_rule(0.)
    0.48
    >>> model.f_w_rule(1.)
    1.878

    :param f_w: empirical expression used in :math:`\\alpha` [dimensionless], see Equation :eq:`f_w_SRK`
    :type f_w: float, derived
    """

    def __init__(self, **kwargs):
        Cubic.__init__(self, sigma=1, epsilon=0, Omega=0.08664, Psi=0.42748, name='SRK', **kwargs)
        self.f_w = self.f_w_rule(self.w)

    def f_w_rule(self, w):
        r"""

        :param w: accentric factor
        :returns: :math:`f_w`, see Equation :eq:`f_w_SRK`
        """
        return 0.480 + 1.574 * w - 0.176 * w * w

    def alpha_expr(self, T_r):
        val = 1. + self.f_w * (1 - T_r ** 0.5)
        return val * val

    def d_ln_alpha_d_ln_Tr(self, T_r):
        T_r_sqrt = T_r ** 0.5
        return (
                self.f_w / (self.f_w * (T_r_sqrt - 1) - 1) / T_r
        )


class PengRobinson(SoaveRedlichKwong):
    r"""Peng-Robinson Equation of State :cite:`peng-robinson`

    This equation of state has the following parameters

    ===========================       ====================================================
            Symbol                                  Value
    ===========================       ====================================================
    :math:`\alpha(T_\text{r})`         :math:`\left[1 + f_w(1-\sqrt{T_\text{r}})\right]^2`
    :math:`\sigma`                               :math:`1 + \sqrt{2}`
    :math:`\epsilon`                             :math:`1 - \sqrt{2}`
    :math:`\Omega`                                :math:`0.07780`
    :math:`\Psi`                                    :math:`0.45724`
    ===========================       ====================================================

    where

    .. math::
        f_w = 0.480 + 1.574\omega - 0.176\omega^2
        :label: f_w_PR

    >>> from realgas.eos.cubic import PengRobinson
    >>> model = PengRobinson(compound_name='Water')
    >>> model.Omega
    0.0778
    >>> model.sigma
    2.4142135
    >>> model.epsilon
    -0.414213
    >>> model.Psi
    0.45724
    >>> model.f_w_rule(0.)
    0.37464
    >>> model.f_w_rule(1.)
    1.64698

    :param f_w: empirical expression used in :math:`\\alpha` [dimensionless?]
    :type f_w: float, derived
    """

    def __init__(self, **kwargs):
        Cubic.__init__(self, sigma=1 + np.sqrt(2), epsilon=1 - np.sqrt(2), Omega=0.07780, Psi=0.45724, name='PR',
                       **kwargs)
        self.f_w = self.f_w_rule(self.w)

    def f_w_rule(self, w):
        r"""

        :param w: accentric factor
        :returns: :math:`f_w`, see Equation :eq:`f_w_PR`
        """
        return 0.37464 + 1.54226 * w - 0.26992 * w * w
