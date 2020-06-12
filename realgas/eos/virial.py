r"""

.. todo::
    merge docs with those in :code:`realgas.partial_molar_properties`.
    There, the definitions of residual properties are displayed and here wee need only
    write simplified forms for the specific equation of state.

Theory
------

The second order virial equation of state is :cite:`Perry`

.. math::
    Z = 1 + B\rho = 1 + \frac{BP}{RT}
    :label: z_virial

Where the composition dependency of :math:`B` is given by the *exact* mixing rule

.. math::
    B = \sum_i \sum_j y_i y_j B_{ij}
    :label: B_mix_expr

where :math:`B_{ij}=B_{ji}`, and :math:`B_{ii}` and :math:`B_{jj}` are virial coefficients for the pure species

In this package, the useful correlation

.. math::
    \frac{BP_\text{c}}{RT_\text{c}} = B^0 + \omega B^1

or

.. math::
    B = \frac{RT_\text{c}}{P_\text{c}}\left(B^0 + \omega B^1\right)
    :label: B_expr

So that, combining Equations :eq:`z_virial` and :eq:`B_expr`,
the compressibility can be calculated from dimensionless quantities as

.. math::
    Z = 1 + \left(B^0 + \omega B^1\right)\frac{P_\text{r}}{T_\text{r}}
    :label: z_dimensionless

where :cite:`Smith2005`

.. math::
    B^0 = 0.083 - \frac{0.422}{T_\text{r}^{1.6}}
    :label: B0_expr

.. math::
    B^1 = 0.139 - \frac{0.172}{T_\text{r}^{4.2}}
    :label: B1_expr

so that

so that the following derivatives can be computed as

.. math::
    \frac{dB^0}{dT_\text{r}} = \frac{0.675}{T_\text{r}^{2.6}}
    :label: dB0dTr

.. math::
    \frac{dB^1}{dT_\text{r}} = \frac{0.722}{T_\text{r}^{5.2}}
    :label: dB1dTr

Which allow the :math:`H^\text{R}`, :math:`S^\text{R}`, and :math:`G^\text{R}` to be readily computed as follows :cite:`Perry`

.. math::
    \frac{G^\text{R}}{RT} = \left(B^0 + \omega B^1\right)\frac{P_\text{r}}{T_\text{r}}
    :label: G_R_RT_virial

.. math::
    \frac{H^\mathrm{R}}{RT} = P_\mathrm{r}\left[
        \frac{B^0}{T_\mathrm{r}} - \frac{\mathrm{d} B^0}{\mathrm{d}T_\mathrm{r}} + \omega\left(\frac{B^1}{T_\mathrm{r}} - \frac{\mathrm{d}B^1}{\mathrm{d}T_\mathrm{r}}\right)
    \right]
    :label: H_R_RT_virial

.. math::
    \frac{S^\text{R}}{R} = -P_\text{r}\left(\frac{d B^0}{d T_\text{r}} - \omega\frac{d B^1}{d T_\text{r}}\right)
    :label: S_R_R_virial


The cross coefficients are calculated as

.. math::
    B_{ij} = \frac{RT_{\text{c}ij}}{P_{\text{c}ij}}\left(B^0 + \omega_{ij}B^1\right)
    :label: B_ij_expr

so that the cross derivatives can be computed as

.. math::
    \frac{dB_{ij}}{dT} = \frac{RT_{\text{c}ij}}{P_{\text{c}ij}}\left(\frac{dB^0}{d T} + \omega_{ij}\frac{dB^1}{dT}\right)

    \frac{dB_{ij}}{dT} = \frac{R}{P_{\text{c}ij}}\left(\frac{dB^0}{d T_{\text{r}ij}} + \omega_{ij}\frac{dB^1}{dT_{\text{r}ij}}\right)

or, equivalently,

.. math::
    \frac{dB_{ij}}{dT_{\text{r}ij}} = \frac{RT_{\text{c}ij}}{P_{\text{c}ij}}\left(\frac{dB^0}{d T_{\text{r}ij}} + \omega_{ij}\frac{dB^1}{dT_{\text{r}ij}}\right)
    :label: dB_dTrij


Fugacity Coefficients
---------------------

The Fugacity coefficients are calculated as

.. math::
    \ln\hat{\phi}_i=\left(\frac{\partial (nG^\text{R}/R/T)}{\partial n_i}\right)_{P,T,n_j}

For the virial equation of state, this becomes :cite:`VanNess1982`

.. math::
    \ln\hat{\phi}_k = \frac{P}{RT}\left[B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)\right]
    :label: ln_phi_i_virial

where *both* :math:`i` and :math:`j` indices run over all species

.. math::
    \delta_{ik} = 2B_{ik} - B_{ii} - B_{kk} = \delta_{ki}
    :label: d_ik

and

.. math::
    \delta_{ii} = 0

Residual Partial Molar Properties
---------------------------------
The partial molar residual free energy of component *k* is

.. math::
    \frac{\bar{G}^\text{R}_k}{RT} = \ln{\hat{\phi}_k} = \frac{P}{RT}\left[B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)\right]
    :label: barGiR


The partial molar residual enthalpy of component *k* is


.. math::
    \begin{align}
        \frac{\bar{H}^\mathrm{R}_k}{RT} &= -T \left(\frac{\partial \ln\hat{\phi}_k}{\partial T}\right)_{P,y} \\
                             &= - \frac{T}{T_c} \left(\frac{\partial \ln\hat{\phi}_k}{\partial T_r}\right)_{P,y} \\
                             &= - \frac{T}{T_c}\left\{
                                \frac{P}{RT}\left[
                                    \frac{\partial B_{kk}}{\partial T_r}
                                     + \frac{1}{2}\sum_i\sum_j y_iy_j\left(
                                                    2\frac{\partial \delta_{ik}}{\partial T_r}
                                                    - \frac{\partial \delta_{ij}}{\partial T_r}
                                                    \right)
                                    \right]
                                - \frac{T_c\ln\hat{\phi}_i}{T}
                             \right\}
    \end{align}

where :math:`\frac{\partial \delta_{ij}}{\partial T_r}` is given by :eq:`eq_d_dij`
so that we obtain

.. math::
    \frac{\bar{H}^\mathrm{R}_k}{RT} = - \frac{P}{RT_\text{c}}\left[
                                    \frac{\partial B_{kk}}{\partial T_r}
                                     + \frac{1}{2}\sum_i\sum_j y_iy_j\left(
                                                    2\frac{\partial \delta_{ik}}{\partial T_r}
                                                    - \frac{\partial \delta_{ij}}{\partial T_r}
                                                    \right)
                                    \right]
                                + \ln\hat{\phi}_i
    :label: barHiR

Since

.. math::
    G^\mathrm{R} = H^\mathrm{R} - T S^\mathrm{R}

In terms of partial molar properties, then

.. math::
    \begin{align}
        \bar{S}_i^\mathrm{R} &= \frac{\bar{H}_i^\mathrm{R} - \bar{G}_i^\mathrm{R}}{T} \\
        \frac{\bar{S}_i^\mathrm{R}}{R} &= \frac{\bar{H}_i^\mathrm{R}}{RT} - \frac{\bar{G}_i^\mathrm{R}}{RT} \\
    \end{align}

By comparing Equation :eq:`barGiR` and :eq:`barHiR` it is observed that

.. math::
    \frac{\bar{S}_k^\mathrm{R}}{R} = - \frac{P}{RT_\text{c}}\left[
                                    \frac{\partial B_{kk}}{\partial T_r}
                                     + \frac{1}{2}\sum_i\sum_j y_iy_j\left(
                                                    2\frac{\partial \delta_{ik}}{\partial T_r}
                                                    - \frac{\partial \delta_{ij}}{\partial T_r}
                                                    \right)
                                    \right]
    :label: barSiR

where :math:`\frac{\partial \delta_{ij}}{\partial T_r}` is given by :eq:`eq_d_dij`



The partial molar residual volume of component *i* is calculated as

.. math::
    \begin{align}
        \frac{\bar{V}_k^\mathrm{R}}{RT} &= \left(\frac{\partial \ln\hat{\phi}_i}{\partial P}\right)_{T,y}\\
            &= \frac{\partial}{\partial P}\left\{\frac{P}{RT}\left[B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)\right]\right\}
    \end{align}

which simplifies to

.. math::
    \frac{\bar{V}_k^\mathrm{R}}{RT} = \frac{1}{RT}\left[B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)\right]
    :label: barViR

from which we obtain the intuitive result of

.. math::
    \bar{V}_k^\mathrm{R} = B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)


"""

import typing

import matplotlib.pyplot as plt
import numpy as np
from chem_util.chem_constants import gas_constant as R

from ..critical_constants import CriticalConstants
from ..input import RealMixture


class Virial:
    """

    :param pow: function for computing power, defaults to numpy.power
    :type pow: callable, optional
    :param exp: function for computing logarithm, defaults to numpy.exp
    :type exp: callable, optional
    """

    def __init__(self, pow: callable = np.power, exp: callable = np.exp):
        self.pow = pow
        self.exp = exp

    def B0_expr(self, T_r):
        """
        :param T_r: Reduced temperature
        :returns: Equation :eq:`B0_expr`
        """
        return 0.083 - 0.422 * self.pow(T_r, -1.6)

    def B1_expr(self, T_r):
        """

        :param T_r: reduced temperature
        :return: Equation eq:`B1_expr`
        """
        return 0.139 - 0.172 * self.pow(T_r, -4.2)

    def d_B0_d_Tr_expr(self, T_r):
        """
        :param T_r: reduced temperature
        :returns: Equation :eq:`dB0dTr`
        """
        return 0.6752 * self.pow(T_r, -2.6)

    def d_B1_d_Tr_expr(self, T_r):
        """
        :param T_r: reduced temperature
        :returns: Equation :eq:`dB1dTr`
        """
        return 0.7224 * self.pow(T_r, -5.2)

    def B_expr(self, T_r, w, T_c, P_c):
        """

        :param T_r: reduced temperature
        :param w: accentric factor
        :param T_c: critical temperautre [K]
        :param P_c: critical pressure [Pa]
        :return: Equation :eq:`B_expr`
        """
        return R * T_c / P_c * (
                self.B0_expr(T_r) + w * self.B1_expr(T_r)
        )

    def ln_hat_phi_k_expr(self, *args):
        raise NotImplementedError

    def hat_phi_i_expr(self, *args):
        r"""expression for fugacity coefficient
        :returns: :math:`\exp\left({\ln{\hat{\phi}_i}}\right)`
        """
        return self.exp(self.ln_hat_phi_k_expr(*args))


class SecondVirial(CriticalConstants, Virial):
    """Virial equation of state for one component. See :cite:`Perry,Smith2005`

    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 pow: callable = np.power, **kwargs):
        """

        :param kwargs: for customized critica constants
        """
        CriticalConstants.__init__(self, dippr_no, compound_name, cas_number, **kwargs)
        Virial.__init__(self, pow=pow)
        self.name = 'secondvirial'

    def H_R_RT_expr(self, P, T):
        """

        :param P: pressure in Pa
        :param T: temperature in K
        :return: Equation :eq:`H_R_RT_virial`
        """
        T_r = T / self.T_c
        return P / self.P_c * (
                self.B0_expr(T_r) / T_r - self.d_B0_d_Tr_expr(T_r) + self.w * (
                self.B1_expr(T_r) / T_r - self.d_B1_d_Tr_expr(T_r))
        )

    def G_R_RT_expr(self, P, T):
        """

        :param P: pressure in Pa
        :param T: Temperautre in K
        :return: Equation :eq:`G_R_RT_virial`
        """
        T_r = T / self.T_c
        P_r = P / self.P_c
        return (self.B0_expr(T_r) + self.w * self.B1_expr(T_r)) * P_r / T_r

    def S_R_R_expr(self, P, T):
        """

        :param P: pressure in Pa
        :param T: temperature in K
        :return: Equation :eq:`S_R_R_virial`
        """
        T_r = T / self.T_c
        P_r = P / self.P_c
        return -P_r * (self.d_B0_d_Tr_expr(T_r) + self.w * self.d_B1_d_Tr_expr(T_r))

    def ln_hat_phi_k_expr(self, P, T):
        r"""logarithm of fugacity coefficient

        .. note::
            single-component version

        In this case, Equation :eq:`ln_phi_i_virial` simplifies to

        .. math::
            \ln\hat{\phi}_i = \frac{PB}{RT}

        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        """
        T_r = T / self.T_c
        return P * self.B_expr(T_r, self.w, self.T_c, self.P_c) / R / T

    def calc_Z_from_units(self, P, T):
        """

        :param P: pressure in Pa
        :param T: temperature in K
        :return: Equation :eq:`z_virial`
        """
        return 1. + self.B_expr(T / self.T_c, self.w, self.T_c, self.P_c) * P / R / T

    def calc_Z_from_dimensionless(self, P_r, T_r):
        """

        :param P_r: reduced pressure, dimensionless
        :param T_r: reduced temperature, dimensionless
        :return: Equation :eq:`z_dimensionless`
        """
        return 1 + (self.B0_expr(T_r) + self.w * self.B1_expr(T_r)) * P_r / T_r

    def plot_Z_vs_P(self, T, P_min, P_max, symbol='o', ax=None, **kwargs):
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
        Z = [self.calc_Z_from_dimensionless(pressure / self.P_c, T / self.T_c) for pressure in P]
        ax.plot(P, Z, symbol, markerfacecolor='None', **kwargs)


class MixingRule(Virial):
    r"""Van der Walls mixing rule

    combining rules from :cite:`Prausnitz1986`


    .. math::
        \omega_{ij} = \frac{\omega_i + \omega_j}{2}
        :label: omega_combine

    .. math::
        T_{\text{c}ij} = \sqrt{T_{\text{c}i}T_{\text{c}j}}(1-k_{ij})
        :label: Tc_combine

    .. math::
        P_{\text{c}ij} = \frac{Z_{\text{c}ij}RT_{\text{c}ij}}{V_{\text{c}ij}}
        :label: Pc_combine

    .. math::
        Z_{\text{c}ij} = \frac{Z_{\text{c}i} + Z_{\text{c}j}}{2}
        :label: Zc_combine

    .. math::
        V_{\text{c}ij} = \left(\frac{V_{\text{c}i}^{1/3} + V_{\text{c}j}^{1/3}}{2}\right)^3
        :label: vc_combine


    """

    def __init__(self, pow: callable = np.power, exp: callable = np.exp):
        Virial.__init__(self, pow, exp)

    def w_ij_rule(self, w_i, w_j):
        r"""

        :param w_i: accentric factor of component *i*
        :param w_j: accentric factor of component *j*
        :return: Equation :eq:`omega_combine`
        """
        return (w_i + w_j) / 2.

    def T_cij_rule(self, T_ci, T_cj, k_ij):
        """

        :param T_ci: critical temperature of component i [K]
        :param T_cj: ** j [K]
        :param k_ij: k_ij parameter
        :return: Equation :eq:`Tc_combine`
        """
        return pow(T_ci * T_cj, 0.5) * (1. - k_ij)

    def P_cij_rule(self, Z_ci, V_ci, T_ci, Z_cj, V_cj, T_cj, k_ij):
        """

        :return: Equation :eq:`Pc_combine`
        """
        Z_cij = self.Z_cij_rule(Z_ci, Z_cj)
        T_cij = self.T_cij_rule(T_ci, T_cj, k_ij=k_ij)
        V_cij = self.V_cij_rule(V_ci, V_cj)
        return Z_cij * R * T_cij / V_cij

    def Z_cij_rule(self, Z_ci, Z_cj):
        """

        :param Z_ci: critical compressibility factor of component *i*
        :param Z_cj: critical compressibility factor of component *j*
        :return: Equation :eq:`Zc_combine`
        """
        return (Z_ci + Z_cj) / 2.

    def V_cij_rule(self, V_ci, V_cj):
        """

        :param V_ci: critical molar volume of component *i* [m**3/mol]
        :param V_cj: critical molar volume of component *j* [m**3/mol]
        :return: Equation eq:`Vc_combine`
        """
        return self.pow(
            (self.pow(V_ci, 1 / 3) + self.pow(V_cj, 1 / 3)) / 2., 3
        )


class SecondVirialMixture(RealMixture, MixingRule):
    r"""Second virial with mixing rule from :class:`.MixingRule`

    .. note::
        can only input both custom critical properties or both from DIPPR--cant have mixed at the moment

    :param pow: function to calculate power, defaults to numpy.power
    :param exp: function to calculate exponential,d efaults to numpy.exp
    :param kwargs: key-word arguments to pass to :class:`.RealMixture`

    """

    def __init__(self,
                 pow: callable = np.power, exp: callable = np.exp,
                 **kwargs):
        MixingRule.__init__(self, pow, exp)
        RealMixture.__init__(self, **kwargs)
        RealMixture.setup(self)

        for i in range(self.num_components):
            I = CriticalConstants(**self.get_point_input(i))
            self.set_point_input(i, **I.__dict__)
        self.check_params()

    def get_w_Tc_Pc(self, i: int, j=None):
        """Returns critical constants for calculation based off of whetner i = j or not

        :returns: (:math:`w`, :math:`T_c`, :math:`P_c`)
        :rtype: tuple
        """
        if j is None or i == j:
            return self.ws[i], self.T_cs[i], self.P_cs[i]

        return (
            self.w_ij_rule(self.ws[i], self.ws[j]),
            self.T_cij_rule(self.T_cs[i], self.T_cs[j], k_ij=self.k_ij[i][j]),
            self.P_cij_rule(
                self.Z_cs[i], self.V_cs[i], self.T_cs[i],
                self.Z_cs[j], self.V_cs[j], self.T_cs[j],
                k_ij=self.k_ij[i][j]
            )
        )

    def B_ij_expr(self, i: int, j: int, T):
        r"""

        :param i: index of first component
        :param j: index of second component
        :param T: temperature [K]
        :returns: Equation :eq:`B_ij_expr`
        """
        w, T_c, P_c = self.get_w_Tc_Pc(i, j)
        T_r = T / T_c
        return self.B_expr(T_r, w, T_c, P_c)

    def B_mix_expr(self, y_k: typing.List[typing.Union[float, typing.Any]], T):
        """

        :param y_k: mole fractions of each component :math:`k`
        :param T: temperature in K
        :returns: Equation :eq:`B_mix_expr`
        """
        return sum(
            y_k[i] * y_k[j] * self.B_ij_expr(i, j, T) for i in range(self.num_components) for j in
            range(self.num_components)
        )

    def calc_Z_from_units(self, y_k: typing.List, P, T):
        """

        :param y_k: mole fractions of each component :math:`k`
        :param P: pressure in Pa
        :param T: temperature in K
        :return: Equation :eq:`z_virial`
        """
        return 1. + self.B_mix_expr(y_k, T) * P / R / T

    def d_ik_expr(self, i: int, k: int, T):
        """

        :param i: index of component *i*
        :param k: index of component *k*
        :param T: temperature [K]
        :return: Equation :eq:`d_ik`
        """
        if i == k:
            return 0.

        return 2. * self.B_ij_expr(i, k, T) - self.B_ij_expr(i, i, T) - self.B_ij_expr(k, k, T)

    def ln_hat_phi_k_expr(self, k: int, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r"""logarithm of fugacity coefficient

        .. note::
            The order of :code:`ys` corresponds to the order of components made during initialization

        :param k: index of component *k*
        :param ys: mole fractions ordered as component order
        :param P: pressure in Pa
        :param T: temperature in K
        :returns: Equation :eq:`ln_phi_i_virial`
        """
        return P / R / T * (
                self.B_ij_expr(k, k, T)
                + 1. / 2. * sum(
            ys[i] * ys[j] * (2. * self.d_ik_expr(i, k, T) - self.d_ik_expr(i, j, T))
            for i in range(self.num_components) for j in range(self.num_components)
        )
        )

    def fugacity_i_expr(self, cas_i: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r"""Fugacity of component i in mixture :math:`f_i=\hat{\phi}_i y_i P`

        :param cas_i: cas number for component of interest
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        """
        assert cas_i in self.cas_numbers, 'Cas number not found!'
        i = self.cas_numbers.index(cas_i)
        return self.hat_phi_i_expr(i, ys, P, T) * ys[i] * P

    def bar_GiR_RT(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r"""Dimensionless residual partial molar free energy of component :math:`i`

        :param cas_k: cas number of component k
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        :returns: Equation :eq:`barGiR`
        """
        assert cas_k in self.cas_numbers, 'Cas number not found'
        k = self.cas_numbers.index(cas_k)
        return self.ln_hat_phi_k_expr(k, ys, P, T)

    def d_Bij_d_Trij(self, i: int, j: int, T):
        """

        :param i: index for component *i*
        :param j: index for componen t *j*
        :param T: temperature in K
        :return: Equation :eq:`dB_dTrij`
        """
        w, T_c, P_c = self.get_w_Tc_Pc(i, j)
        T_r = T / T_c
        return R * T_c / P_c * (self.d_B0_d_Tr_expr(T_r) + w * self.d_B1_d_Tr_expr(T_r))

    def d_dij_d_Tr(self, i, j, T):
        r"""
        .. math::
            \frac{\partial \delta_{ij}}{\partial T_{\text{r}ij}} = 2\frac{\partial B_{ij}}{\partial T_{\text{r}ij}}-\frac{\partial B_{ii}}{\partial T_{\text{r}ij}}-\frac{\partial B_{jj}}{\partial T_{\text{r}ij}}
            :label: eq_d_dij

        .. todo::
            test this with symbolic differentiation of d_ik expression

        :param i: index for component *i*
        :param j: index for componen t *j*
        :param T: temperature [K]
        """
        if i == j:
            return 0.

        return 2. * self.d_Bij_d_Trij(i, j, T) \
               - self.d_Bij_d_Trij(i, i, T) \
               - self.d_Bij_d_Trij(j, j, T)

    def bar_HiR_RT(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r"""Dimensionless residual partial molar enthalpy of component :math:`i`

        :param cas_k: cas number for component of interest
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        """
        assert cas_k in self.cas_numbers, 'Cas number not found'
        k = self.cas_numbers.index(cas_k)
        T_c = self.T_cs[k]
        return self.ln_hat_phi_k_expr(k, ys, P, T) - P / R / T_c * (
                self.d_Bij_d_Trij(k, k, T)
                + 1. / 2. * sum(
            ys[i] * ys[j] * (2. * self.d_dij_d_Tr(i, k, T) - self.d_dij_d_Tr(i, j, T))
            for i in range(self.num_components) for j in range(self.num_components)
        )
        )

    def bar_HiR_star(self, T_star, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r""" Returns the scaled entropy useful for computations :math:`\bar{H}_i^{\text{R},\star}`,
        as defined in Equation :eq:`barHiRstar_definition`

        :param cas_k: cas number of component *k*
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        """
        return T_star * self.bar_HiR_RT(cas_k, ys, P, T)

    def bar_SiR_R(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r"""Dimensionless residual partial molar entropy of component :math:`i`

        :param cas_k: cas number for component of interest
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        :returns: Equation :eq:`barSiR`
        """
        assert cas_k in self.cas_numbers, 'Cas number not found'
        k = self.cas_numbers.index(cas_k)
        T_c = self.T_cs[k]
        return -P / R / T_c * (
                self.d_Bij_d_Trij(k, k, T)
                + sum(
            1. / 2. * ys[i] * ys[j] * (2. * self.d_dij_d_Tr(i, k, T) - self.d_dij_d_Tr(i, j, T))
            for i in range(self.num_components) for j in range(self.num_components)
        )
        )

    def bar_ViR_RT(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        r""" residual Partial molar volume for component *i*

        .. note::
            This expression does not depend on :math:`P`

        :param cas_k: cas number for component of interest
        :param ys: mole fractions
        :param P: pressure in Pa
        :param T: temperature in K
        :returns: Equation :eq:`barViR`
        """
        assert cas_k in self.cas_numbers, 'Cas number not found'
        k = self.cas_numbers.index(cas_k)
        return 1. / R / T(
            self.B_ij_expr(k, k, T)
            + 1. / 2. * sum(
                ys[i] * ys[j] * (2. * self.d_ik_expr(i, k, T) - self.d_ik_expr(i, j, T))
                for i in range(self.num_components) for j in range(self.num_components)
            )
        )

    def M_R_dimensionless(self, method: callable, ys: typing.List[typing.Union[float, typing.Any]], P: float, T: float):
        """Residual property of :math:`X` for mixture.

        Similar to Equation :eq:`residual_molar` but in dimensionless form

        :param method: function to compute partial molar property of compound
        :type method: callable
        """
        return sum(
            ys[i] * method(self.cas_numbers[i], ys, P, T) for i in range(self.num_components)
        )

    def S_R_R(self, *args):
        r"""Residual entropy of mixture :math:`S^\mathrm{R}/R`"""
        return self.M_R_dimensionless(self.bar_SiR_R, *args)

    def H_R_RT(self, *args):
        r"""Residual enthalpy of mixture :math:`H^\mathrm{R}/R/T`"""
        return self.M_R_dimensionless(self.bar_HiR_RT, *args)

    def G_R_RT(self, *args):
        r"""Residual free energy of mixture :math:`G^\mathrm{R}/R/T`"""
        return self.M_R_dimensionless(self.bar_GiR_RT, *args)

    def plot_residual_HSG(self, P, T, ax=None, fig=None) -> typing.Tuple[plt.figure, plt.subplot]:
        """Plot dimensionless residual properties as a function of mole fraction

        :param P: pressure in Pa
        :param T: Temperature in K
        :param ax: matplotlib ax, defaults to None
        :param fig: matplotlib figure, defautls to None
        """

        assert self.num_components == 2, 'Plotting only implemented for binary mixtures'

        if ax is None:
            if fig is None:
                fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_ylabel('Dimensionless')
            ax.set_xlabel('Gas-Phase Mole fraction %s' % self.compound_names[0])

        y_1 = np.linspace(0., 1.)
        HR_RT = np.array(list(self.H_R_RT([i, 1. - i], P, T) for i in y_1))
        GR_RT = np.array(list(self.G_R_RT([i, 1. - i], P, T) for i in y_1))
        SR_R = np.array(list(self.S_R_R([i, 1. - i], P, T) for i in y_1))
        ax.plot(y_1, HR_RT, '-', label='$H^\\mathrm{R}/(RT)$')
        ax.plot(y_1, GR_RT, '--', label='$G^\\mathrm{R}/(RT)$')
        ax.plot(y_1, SR_R, '-.', label='$S^\\mathrm{R}/R$')
        ax.set_title('Residual properties at T = %5.2f K and P = %4.3e Pa' % (T, P))
        ax.legend()
        return fig, ax
