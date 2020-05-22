r"""

The second order virial equation of state is :cite:`Perry`

.. math::
    Z = 1 + B\rho = 1 + \frac{BP}{RT}
    :label: z_virial

Where the composition dependency of :math:`B` is given by the *exact* mixing rule

.. math::
    B = \sum_i \sum_j y_i y_j B_{ij}

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

so that the cross derivatives can be computed as

.. math::
    \frac{dB_{ij}}{dT} = \frac{RT_{\text{c}ij}}{P_{\text{c}ij}}\left(\frac{dB^0}{d T} + \omega_{ij}\frac{dB^1}{dT}\right)

or, equivalently,

.. math::
    \frac{dB_{ij}}{dT} = \frac{R}{P_{\text{c}ij}}\left(\frac{dB^0}{d T_{\text{r}ij}} + \omega_{ij}\frac{dB^1}{dT_{\text{r}ij}}\right)

Combining Rules
---------------
The following combining rules are used


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

Fugacity Coefficients
---------------------

The Fugacity coefficients are calculated as

.. math::
    \ln\hat{\phi}_i=\left(\frac{\partial (nG^\text{R}/R/T)}{\partial n_i}\right)_{P,T,n_j}

For the virial equation of state, this becomes :cite:`VanNess1982`

.. math::
    \ln\hat{\phi}_k = \frac{P}{RT}\left[B_{kk} + \frac{1}{2}\sum_i\sum_jy_iy_j\left(2\delta_{ik} - \delta_{ij}\right)\right]
    :label: ln_phi_i_virial

where

.. math::
    \delta_{ik} = 2B_{ik} - B_{ii} - B_{kk} = \delta_{ki}

and

.. math::
    \delta_{ii} = 0

"""

from ..critical_constants import CriticalConstants
from chem_util.chem_constants import gas_constant
import numpy as np
import matplotlib.pyplot as plt


class Virial:
    """

    :param R: gas constant, set to SI units
    :type R: float, hard-coded
    :param pow: function for computing power, defaults to numpy.power
    :type pow: callable, optional
    :param exp: function for computing logarithm, defaults to numpy.exp
    :type exp: callable, optional
    """

    def __init__(self, pow: callable = np.power, exp: callable = np.exp):
        self.pow = pow
        self.exp = exp
        self.R = gas_constant

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
        return self.R * T_c / P_c * (
                self.B0_expr(T_r) + w * self.B1_expr(T_r)
        )

    def ln_hat_phi_i_expr(self, *args):
        raise NotImplementedError

    def hat_phi_i_expr(self, *args):
        r"""expression for fugacity coefficient
        :returns: :math:`\exp{\ln{\hat{\phi}_i}}
        """
        return self.exp(self.ln_hat_phi_i_expr(*args))


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
        :return: Equation :eq:`S_R_RT_virial`
        """
        T_r = T / self.T_c
        P_r = P / self.P_c
        return -P_r * (self.d_B0_d_Tr_expr(T_r) + self.w * self.d_B1_d_Tr_expr(T_r))

    def ln_hat_phi_i_expr(self, P, T):
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
        return P * self.B_expr(T_r, self.w, self.T_c, self.P_c) / self.R / T

    def calc_Z_from_units(self, P, T):
        """

        :param P: pressure in Pa
        :param T: temperature in K
        :return: Equation :eq:`z_virial`
        """
        return 1. + self.B_expr(T / self.T_c, self.w, self.T_c, self.P_c) * P / self.R / T

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


class BinarySecondVirial(CriticalConstants, Virial):
    r"""Second virial with combining rules from :cite:`Prausnitz1986`


    .. todo::
        Fix docs here

    .. math::
        w_{ij} = \frac{w_i + w_j}{2}

    .. math::
        T_{\mathrm{c},ij} = \sqrt{T_{\mathrm{c},i}T_{\mathrm{c},j}}(1-k_{ij})
        :label: eq_Tcij

    .. math::
        P_{\mathrm{c},ij} = \frac{Z_{\mathrm{c},ij}RT_{\mathrm{c},ij}}{V_{\mathrm{c},ij}}

    where

    .. math::
        \begin{align}
            Z_{\mathrm{c},ij} &= \frac{Z_{\mathrm{c},i} + Z_{\mathrm{c},j}}{2} \\\
            V_{\mathrm{c},ij} &= \left(\frac{V_{\mathrm{c},i}^{1/3} + V_{\mathrm{c},j}^{1/3}}{2}\right)^{3}
        \end{align}

    :param k_ij: equation of state mixing rule in calculation of critical temperautre, see Equation :eq:`eq_Tcij`. When :math:`i=j` and for chemical similar species, :math:`k_{ij}=0`. Otherwise, it is a small (usually) positive number evaluated from minimal :math:`PVT` data or, in the absence of data, set equal to zero.
    :type k_ij: float
    :param cas_pairs: pairs of cas registry numbers, derived from cas_numbers calculated
    :type cas_pairs: list(tuple(str))
    :param other_cas: map from one cas number to other
    :type other_cas: dict

    """

    def __init__(self,
                 i_kwargs=None, j_kwargs=None,
                 k_ij=0., pow: callable = np.power, exp: callable = np.exp):
        """

        :param i_kwargs: kwargs for component i passed to critical constants
        :param j_kwargs: kwargs for component j passed to critical constants
        """
        if i_kwargs is None:
            i_kwargs = {}
        if j_kwargs is None:
            j_kwargs = {}
        self.i = CriticalConstants(**i_kwargs)
        self.j = CriticalConstants(**j_kwargs)
        Virial.__init__(self, pow, exp)
        assert self.i.cas_number != self.j.cas_number, 'Errors anticipated when cas numbers are equal if other properties are not'

        self.k_ij = k_ij
        self.w_ij = (self.i.w + self.j.w) / 2.
        self.T_c_ij = pow(self.i.T_c * self.j.T_c, 0.5) * (1. - self.k_ij)
        self.Z_c_ij = (self.i.Z_c + self.j.Z_c) / 2.
        self.V_c_ij = pow(
            (pow(self.i.V_c, 1 / 3) + pow(self.j.V_c, 1 / 3)) / 2., 3
        )
        self.P_c_ij = self.Z_c_ij * self.R * self.T_c_ij / self.V_c_ij
        self.other_cas = {
            self.i.cas_number: self.j.cas_number,
            self.j.cas_number: self.i.cas_number,
        }
        self.cas_pairs = []
        for x in self.other_cas.keys():
            for y in self.other_cas.values():
                self.cas_pairs.append((x, y))

    def get_w_Tc_Pc(self, cas_i, cas_j=None):
        """Returns critical constants for calculation based off of whetner i = j or not

        :returns: (:math:`w`, :math:`T_c`, :math:`P_c`)
        :rtype: tuple
        """
        if cas_i == self.i.cas_number and (cas_i == cas_j or cas_j is None):
            return self.i.w, self.i.T_c, self.i.P_c
        if cas_i == self.j.cas_number and (cas_i == cas_j or cas_j is None):
            return self.j.w, self.j.T_c, self.j.P_c
        if cas_i == self.other_cas[cas_j]:
            return self.w_ij, self.T_c_ij, self.P_c_ij

    def B_ij_expr(self, cas_1, cas_2, T):
        r"""
        Returns :math:`B_{ij}` considering that :code:`i=j` or :code:`i!=j`

        If :code:`i=j`, find the component :math:`k` for which :code:`k=i=j` (all cas no's equal).
        Then return :math:`B_{kk}` where

        .. math::
            B_{kk}=\frac{RT_{\mathrm{c},k}}{P_{\mathrm{c},k}}\left[
            B^0\left(\frac{T}{T_{\mathrm{c},k}}\right) + \omega_k\left(\frac{T}{T_{\mathrm{c},k}}\right)
            \right]

        Otherwise, if :code:`i!=j`, return :math:`B_{ij}`, where

        .. math::
            B_{ij}=\frac{RT_{\mathrm{c},ij}}{P_{\mathrm{c},ij}}\left[
            B^0\left(\frac{T}{T_{\mathrm{c},ij}}\right) + \omega_{ij}\left(\frac{T}{T_{\mathrm{c},ij}}\right)
            \right]

        This is implemented in a simplified fashion usintg :meth:`BinarySecondVirial.get_w_Tc_Pc`
        and then calling the generic expression for :math:`B`


        :param cas_1: cas number of first component
        :type cas_1: str
        :param cas_2: cas number of second component
        :type cas_2: str
        :param T: temperature [K]
        """
        w, T_c, P_c = self.get_w_Tc_Pc(cas_1, cas_2)
        T_r = T / T_c
        return self.B_expr(T_r, w, T_c, P_c)

    def B_mix_expr(self, y_k, T):
        """

        :param y_k: mole fractions of each component :math:`k`
        :type y_k: dict[:attr:`cas_numbers`, float]
        :param T: temperature in K
        :return: :math:`B` [m^3/mol], where :math:`Z = 1 + BP/R/T`
        """
        return sum(
            y_k[i] * y_k[j] * self.B_ij_expr(i, j, T) for i, j in self.cas_pairs
        )

    def calc_Z(self, y_k, P, T):
        """

        :param y_k: mole fractions of each component :math:`k`
        :type y_k: dict
        :param P: pressure in Pa
        :param T: temperature in K
        :return: :math:`Z` [mol/m^3], where :math:`Z = 1 + BP/R/T`
        """
        return 1. + self.B_mix_expr(y_k, T) * P / self.R / T

    def d_ij_expr(self, T):
        """
        .. math::
            \\delta_{ij} = 2B_{ij} - B_{ii} - B_{jj}
            :label: eq_dij

        :param T: temperature [K]
        :return: :math:`\\delta_{ij}` [m**3/mol]
        """
        return 2. * self.B_ij_expr(self.i.cas_number, self.j.cas_number, T) \
               - self.B_ij_expr(self.i.cas_number, self.i.cas_number, T) \
               - self.B_ij_expr(self.j.cas_number, self.j.cas_number, T)

    def ln_hat_phi_i_expr(self, cas_i, y_i, P, T):
        r"""logarithm of fugacity coefficient

        .. math::
            \ln\hat{\phi}_i = \frac{P}{RT}\left[B_{ii} + (1-y_i)^2\delta_{ij}\right]

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float

        where :math:`\delta_{ij}` is given by Equation :eq:`eq_dij`
        """
        return P * (
                self.B_ij_expr(cas_i, cas_i, T) - (1. - y_i) * (1. - y_i) * self.d_ij_expr(T)
        ) / self.R / T

    def fugacity_i_expr(self, cas_i, y_i, P, T):
        r"""Fugacity of component i in mixture :math:`f_i=\hat{\phi}_i y_i P`

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        """
        return self.hat_phi_i_expr(cas_i, y_i, P, T) * y_i * P

    def bar_GiR_RT(self, *args):
        r"""Dimensionless residual partial molar free energy of component :math:`i`

        .. math::
            \frac{\bar{G}^\mathrm{R}_i}{RT} = \ln\hat{\phi}_i
            :label: eq_GiR

        """
        return self.ln_hat_phi_i_expr(*args)

    def d_Bij_d_Tr(self, cas_i, cas_j, T):
        w, T_c, P_c = self.get_w_Tc_Pc(cas_i, cas_j)
        T_r = T / T_c
        return self.R * T_c / P_c * (self.d_B0_d_Tr_expr(T_r) + w * self.d_B1_d_Tr_expr(T_r))

    def d_dij_d_Tr(self, T):
        r"""
        .. math::
            \frac{\partial \delta_{ij}}{\partial T_r} = 2\frac{\partial B_{ij}}{\partial T_r}-\frac{\partial B_{ii}}{\partial T_r}-\frac{\partial B_{jj}}{\partial T_r}
            :label: eq_d_dij

        :param T: temperature [K]
        """
        return 2. * self.d_Bij_d_Tr(self.i.cas_number, self.j.cas_number, T) \
               - self.d_Bij_d_Tr(self.i.cas_number, self.i.cas_number, T) \
               - self.d_Bij_d_Tr(self.j.cas_number, self.j.cas_number, T)

    def Tstar_d_lnphi_dTstar(self, cas_i, y_i, P, T):
        r"""Returns

        .. math::
            \begin{align}
                T^\star \frac{\partial \ln{\hat{\phi}_i}}{\partial T^\star}
                    &=  \frac{T T_\text{ref}}{T_\text{ref}}\frac{\partial \ln{\hat{\phi}_i}}{\partial T}
                    = \frac{T}{T_\text{c}}\frac{\partial \ln{\hat{\phi}_i}}{\partial T_\text{r}} \\
                &= \frac{T}{T_\text{c}}\left[
                                        \frac{P}{RT}\left(
                                            \frac{\partial B_{ii}}{\partial T_r} + (1-y_i)^2\frac{\partial \delta_{ij}}{\partial T_r}
                                            \right)
                                        - \frac{T_c\ln\hat{\phi}_i}{T}
                \right] \\
                &= \frac{P}{R T_\text{c}}\left(
                                            \frac{\partial B_{ii}}{\partial T_r} + (1-y_i)^2\frac{\partial \delta_{ij}}{\partial T_r}
                                            \right)
                                        - \ln\hat{\phi}_i \\
                &= -\frac{\bar{H}_i^\text{R}}{RT}
            \end{align}

        where :math:`\frac{\partial \delta_{ij}}{\partial T_r}` is given by :eq:`eq_d_dij`

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        """
        return -self.bar_HiR_RT(cas_i, y_i, P, T)

    def bar_HiR_RT(self, cas_i, y_i, P, T):
        r"""Dimensionless residual partial molar enthalpy of component :math:`i`

        .. math::
            \begin{align}
                \frac{\bar{H}^\mathrm{R}_i}{RT} &= -T \left(\frac{\partial \ln\hat{\phi}_i}{\partial T}\right)_{P,y} \\
                                     &= - \frac{T}{T_c} \left(\frac{\partial \ln\hat{\phi}_i}{\partial T_r}\right)_{P,y} \\
                                     &= - \frac{T}{T_c}\left(
                                        \frac{P}{RT}\left[
                                            \frac{\partial B_{ii}}{\partial T_r} + (1-y_i)^2\frac{\partial \delta_{ij}}{\partial T_r}
                                            \right]
                                        - \frac{T_c\ln\hat{\phi}_i}{T}
                                     \right)
            \end{align}

        where :math:`\frac{\partial \delta_{ij}}{\partial T_r}` is given by :eq:`eq_d_dij`
        so that we obtain

        .. math::
            \frac{\bar{H}^\mathrm{R}_i}{RT} = -\frac{P}{RT_c}\left[
                        \frac{\partial B_{ii}}{\partial T_r} + (1-y_i)^2\frac{\partial \delta_{ij}}{\partial T_r}
                                            \right]
                                        + \ln\hat{\phi}_i
            :label: eq_HiR

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        """
        w, T_c, P_c = self.get_w_Tc_Pc(cas_i)
        return -P / self.R / T_c * (self.d_Bij_d_Tr(cas_i, cas_i, T) + (1. - y_i) * (1. - y_i) * self.d_dij_d_Tr(
            T)) + self.ln_hat_phi_i_expr(cas_i, y_i, P, T)

    def bar_SiR_R(self, cas_i, y_i, P, T):
        r"""Dimensionless residual partial molar entropy of component :math:`i`

        Since

        .. math::
            G^\mathrm{R} = H^\mathrm{R} - T S^\mathrm{R}

        In terms of partial molar properties, then

        .. math::
            \begin{align}
                \bar{S}_i^\mathrm{R} &= \frac{\bar{H}_i^\mathrm{R} - \bar{G}_i^\mathrm{R}}{T} \\
                \frac{\bar{S}_i^\mathrm{R}}{R} &= \frac{\bar{H}_i^\mathrm{R}}{RT} - \frac{\bar{G}_i^\mathrm{R}}{RT} \\
            \end{align}

        By comparing Equation :eq:`eq_GiR` and :eq:`eq_HiR` it is observed that

        .. math::
            \frac{\bar{S}_i^\mathrm{R}}{R} = -\frac{P}{RT_c}\left[
                        \frac{\partial B_{ii}}{\partial T_r} + (1-y_i)^2\frac{\partial \delta_{ij}}{\partial T_r}
                                            \right]

        where :math:`\frac{\partial \delta_{ij}}{\partial T_r}` is given by :eq:`eq_d_dij`

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        """
        w, T_c, P_c = self.get_w_Tc_Pc(cas_i)
        return -P / self.R / T_c * (
                self.d_Bij_d_Tr(cas_i, cas_i, T) + (1. - y_i) * (1. - y_i) * self.d_dij_d_Tr(T)
        )

    def bar_ViR_RT(self, cas_i, y_i, P, T):
        r""" residual Partial molar volume for component *i*

        .. math::
            \begin{align}
                \frac{\bar{V}_i^\mathrm{R}}{RT} &= \left(\frac{\partial \ln\hat{\phi}_i}{\partial P}\right)_{T,y}\\
                &= \frac{B_{ii} + (1-y_i)^2\delta_{ij}}{RT}
            \end{align}

        .. note::
            This expression does not depend on :math:`P`

        :param cas_i: cas number for component of interest
        :type cas_i: str
        :param y_i: mole fraction of component of interest
        :type y_i: float
        :param P: pressure in Pa
        :type P: float
        :param T: temperature in K
        :type T: float
        :return: :math:`\bar{V}_i^\mathrm{R}/R/T`
        """
        return (self.B_ij_expr(cas_i, cas_i, T) + (1 - y_i) * (1 - y_i) * self.d_ij_expr(T)) / self.R / T

    def X_R_dimensionless(self, method: callable, cas_i: str, y_i: float, P: float, T: float):
        """Residual property of :math:`X` for mixture.

        :param method: function to compute partial molar property of compound
        :type method: callable
        """
        y_j = 1. - y_i
        cas_j = self.other_cas[cas_i]
        return (
                y_i * method(cas_i, y_i, P, T) + y_j * method(cas_j, y_j, P, T)
        )

    def S_R_R(self, *args):
        """Residual entropy of mixture :math:`S^\mathrm{R}`"""
        return self.X_R_dimensionless(self.bar_SiR_R, *args)

    def H_R_RT(self, *args):
        """Residual enthalpy of mixture :math:`H^\mathrm{R}`"""
        return self.X_R_dimensionless(self.bar_HiR_RT, *args)

    def G_R_RT(self, *args):
        """Residual free energy of mixture :math:`G^\mathrm{R}`"""
        return self.X_R_dimensionless(self.bar_GiR_RT, *args)

    def plot_residual_HSG(self, P, T, ax=None, fig=None):
        if ax is None:
            if fig is None:
                fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_ylabel('Dimensionless')
            ax.set_xlabel('Gas-Phase Mole fraction %s' % self.i.compound_name)

        y_1 = np.linspace(0., 1.)
        HR_RT = np.array(list(self.H_R_RT(self.i.cas_number, i, P, T) for i in y_1))
        GR_RT = np.array(list(self.G_R_RT(self.i.cas_number, i, P, T) for i in y_1))
        SR_R = np.array(list(self.S_R_R(self.i.cas_number, i, P, T) for i in y_1))
        ax.plot(y_1, HR_RT, '-', label='$H^\\mathrm{R}/(RT)$')
        ax.plot(y_1, GR_RT, '--', label='$G^\\mathrm{R}/(RT)$')
        ax.plot(y_1, SR_R, '-.', label='$S^\\mathrm{R}/R$')
        ax.set_title('Residual properties at T = %5.2f K and P = %4.3e Pa' % (T, P))
        ax.legend()
        return fig, ax
