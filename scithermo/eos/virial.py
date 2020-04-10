from scithermo.chem_constants import R_si_units
from scithermo.critical_constants import CriticalConstants
import numpy as np
import matplotlib.pyplot as plt


class Virial:
    """
    :param R: gas constant, set to SI units
    :type R: float, hard-coded
    :param pow: function for computing power, defaults to :ref:`np.power`
    :type pow: callable, optional
    """
    def __init__(self, pow: callable=np.power):
        self.R = R_si_units
        self.pow = pow

    def B0_expr(self, T_r):
        return 0.083 - 0.422*self.pow(T_r, -1.6)

    def B1_expr(self, T_r):
        return 0.139 - 0.172*self.pow(T_r, -4.2)

    def d_B0_d_Tr_expr(self, T_r):
        """
        .. math::
            \\frac{\\mathrm{d} B^0}{\\mathrm{d}T_\\mathrm{r}}
        """
        return 0.6752 * self.pow(T_r, -2.6)

    def d_B1_d_Tr_expr(self, T_r):
        """
        .. math::
            \\frac{\\mathrm{d} B^1}{\\mathrm{d}T_\\mathrm{r}}
        """
        return 0.7224 * self.pow(T_r, -5.2)

    def B_expr(self, T_r, w, T_c, P_c):
        return self.R*T_c/P_c*(
            self.B0_expr(T_r) + w*self.B1_expr(T_r)
        )

    def calc_Z_from_units(self, P, T, w, T_c, P_c):
        return 1. + self.B_expr(T/T_c, w, T_c, P_c)*P/self.R/T

    def calc_Z_from_dimensionless(self, w, P_r, T_r):
        return 1 + (self.B0_expr(T_r) + w*self.B1_expr(T_r))*P_r/T_r


class SecondVirial(CriticalConstants, Virial):
    """Virial equation of state for one component. See :cite:`Perry,Smith2005`

    """
    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 pow: callable=np.power):
        """

        :param kwargs: used in :ref:`CpIdealGas`
        """
        CriticalConstants.__init__(self, dippr_no, compound_name, cas_number)
        Virial.__init__(self, pow=pow)

    def d_B_d_Tr(self, T):
        T_r = T / self.T_c
        return self.d_B0_d_Tr_expr(T_r) + self.w*self.d_B1_d_Tr_expr(T_r)

    def H_R_RT_expr(self, P, T):
        """Dimensionless residual enthalpy

        .. math::
            \\frac{H^\\mathrm{R}}{RT} = P_\\mathrm{r}\\left[
                \\frac{B^0}{T_\\mathrm{r}} - \\frac{\\mathrm{d} B^0}{\\mathrm{d}T_\\mathrm{r}} + \\omega\\left(\\frac{B^1}{T_\\mathrm{r}} - \\frac{\\mathrm{d}B^1}{\\mathrm{d}T_\\mathrm{r}}\\right)
            \\right]

        :return: Expression for residual enthalpy (divided by RT) -- dimensionless
        """
        T_r = T/self.T_c
        return P/self.P_c*(
            self.B0_expr(T_r)/T_r - self.d_B0_d_Tr_expr(T_r) + self.w*(self.B1_expr(T_r)/T_r - self.d_B1_d_Tr_expr(T_r) )
        )

    def G_R_RT_expr(self, P, T):
        """Dimensionless residual gibbs

        .. math::
            \\frac{G^\\mathrm{R}}{RT} = (B^0 + \\omega B^1)\\frac{P_\\mathrm{r}}{T_\\mathrm{r}}

        :return: Expression for residual gibbs free (divided by RT) -- dimensionless
        """
        T_r = T / self.T_c
        P_r = P / self.P_c
        return (self.B0_expr(T_r) + self.w * self.B1_expr(T_r)) * P_r / T_r

    def S_R_R_expr(self, P, T):
        """Dimensionless residual entropy

        .. math::
            \\frac{S^\\mathrm{R}}{R} = -P_\\mathrm{r}\\left(
                \\frac{\\mathrm{d} B^0}{\\mathrm{d}T_\\mathrm{r}} + \\omega\\frac{\\mathrm{d}B^1}{\\mathrm{d}T_\\mathrm{r}}
            \\right)

        :return: Expression for residual entropy (divided by R) -- dimensionless
        """
        T_r = T / self.T_c
        P_r = P / self.P_c
        return -P_r * (self.d_B0_d_Tr_expr(T_r) + self.w*self.d_B1_d_Tr_expr(T_r))

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
        Z = [self.calc_Z_from_dimensionless(self.w, pressure/self.P_c, T/self.T_c) for pressure in P]
        ax.plot(P, Z, symbol, markerfacecolor='None', **kwargs)


class BinarySecondVirial(CriticalConstants, Virial):
    r"""Second virial with combining rules from :cite:`Prausnitz1986`

    .. math::
        \begin{align}
            w_{ij} &= \frac{w_i + w_j}{2} \\
            T_{\mathrm{c},ij} &= \sqrt{T_{\mathrm{c},c1}T_{\mathrm{c},c2}}(1-k_{ij}) :label:eq_Tcij\\
            P_{\mathrm{c},ij} &= \frac{Z_{\mathrm{c},ij}RT_{\mathrm{c},ij}}{V_{\mathrm{c},ij}} \\
        \end{align}

    where

    .. math::
        \begin{align}
            Z_{\mathrm{c},ij} &= \frac{Z_{\mathrm{c},c1} + Z_{\mathrm{c},c2}}{2} \\\
            V_{\mathrm{c},ij} &= \left(\frac{V_{\mathrm{c},c1}^{1/3} + V_{\mathrm{c},c2}^{1/3}}{2}\right)^{1/3}
        \end{align}

    :param k_ij: equation of state mixing rule in calculation of critical temperautre, see Equation :eq:`eq_Tcij`. When :math:`c1=c2` and for chemical similar species, :math:`k_{ij}=0`. Otherwise, it is a small (usually) positive number evaluated from minimal :math:`PVT` data or, in the absence of data, set equal to zero.
    :type k_ij: float
    :param cas_pairs: pairs of cas registry numbers, derived from cas_numbers calculated
    :type cas_pairs: list(tuple(str))
    :param other_cas: map from one cas number to other
    :type other_cas: dict

    """
    def __init__(self,
                 dippr_no_i: str = None, compound_name_i: str = None, cas_number_i: str = None,
                 dippr_no_j: str = None, compound_name_j: str = None, cas_number_j: str = None,
                 k_ij = 0., pow: callable=np.power):
        self.c1 = CriticalConstants(dippr_no=dippr_no_i, compound_name=compound_name_i, cas_number=cas_number_i)
        self.c2 = CriticalConstants(dippr_no=dippr_no_j, compound_name=compound_name_j, cas_number=cas_number_j)
        Virial.__init__(self, pow)

        self.k_ij = k_ij
        self.w_ij = (self.c1.w + self.c2.w) / 2.
        self.T_c_ij = pow(self.c1.T_c * self.c2.T_c, 0.5) * (1. - self.k_ij)
        self.Z_c_ij = (self.c1.Z_c + self.c2.Z_c) / 2.
        self.V_c_ij = pow(
            (pow(self.c1.V_c, 1 / 3) + pow(self.c2.V_c, 1 / 3)) / 2., 1 / 3
        )
        self.P_c_ij = self.Z_c_ij*self.R*self.T_c_ij/self.V_c_ij
        self.other_cas = {
            self.c1.cas_number: self.c2.cas_number,
            self.c2.cas_number: self.c1.cas_number,
        }
        self.cas_pairs = []
        for x in self.other_cas.keys():
            for y in self.other_cas.values():
                self.cas_pairs.append((x, y))

    def get_w_Tc_Pc(self, cas_i, cas_j=None):
        """Returns critical constants for calculation based off of whetner c1 = c2 or not

        :returns: (:math:`w`, :math:`T_c`, :math:`P_c`)
        :rtype: tuple
        """
        if cas_i == self.c1.cas_number and (cas_i == cas_j or cas_j is None):
            return self.c1.w, self.c1.T_c, self.c1.P_c
        if cas_i == self.c2.cas_number and (cas_i == cas_j or cas_j is None):
            return self.c2.w, self.c2.T_c, self.c2.P_c
        if cas_i == self.other_cas[cas_j]:
            return self.w_ij, self.T_c_ij, self.P_c_ij

    def B_ij_expr(self, cas_1, cas_2, T):
        r"""
        Returns :math:`B_{ij}` considering that :code:`c1=c2` or :code:`c1!=c2`

        If :code:`c1=c2`, find the component :math:`k` for which :code:`k=c1=c2` (all cas no's equal).
        Then return :math:`B_{kk}` where

        .. math::
            B_{kk}=\frac{RT_{\mathrm{c},k}}{P_{\mathrm{c},k}}\left[
            B^0\left(\frac{T}{T_{\mathrm{c},k}}\right) + \omega_k\left(\frac{T}{T_{\mathrm{c},k}}\right)
            \right]

        Otherwise, if :code:`c1!=c2`, return :math:`B_{ij}`, where

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
        :type y_k: dict
        :param T: temperature in K
        :return: :math:`B` [m^3/mol], where :math:`Z = 1 + BP/R/T`
        """
        return sum(
           y_k[i]*y_k[j]*self.B_ij_expr(i, j, T) for i, j in self.cas_pairs
        )

    def calc_Z(self, y_k, P, T):
        """

        :param y_k: mole fractions of each component :math:`k`
        :type y_k: dict
        :param P: pressure in Pa
        :param T: temperature in K
        :return: :math:`Z` [mol/m^3], where :math:`Z = 1 + BP/R/T`
        """
        return 1. + self.B_mix_expr(y_k, T)*P/self.R/T

    def d_ij_expr(self, T):
        """
        .. math::
            \\delta_{ij} = 2B_{ij} - B_{ii} - B_{jj}
            :label: eq_dij

        :param T: temperature [K]
        :return: :math:`\\delta_{ij}` [m**3/mol]
        """
        return 2.*self.B_ij_expr(self.c1.cas_number, self.c2.cas_number, T) \
                - self.B_ij_expr(self.c1.cas_number, self.c1.cas_number, T) \
                - self.B_ij_expr(self.c2.cas_number, self.c2.cas_number, T)

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
        return P*(
                self.B_ij_expr(cas_i, cas_i, T) - (1.-y_i)*(1.-y_i)*self.d_ij_expr(T)
        )/self.R/T

    def bar_GiR_RT(self, *args):
        r"""Dimensionless residual partial molar free energy of component :math:`c1`

        .. math::
            \frac{\bar{G}^\mathrm{R}_i}{RT} = \ln\hat{\phi}_i
            :label: eq_GiR

        """
        return self.ln_hat_phi_i_expr(*args)

    def d_Bij_d_Tr(self, cas_i, cas_j, T):
        w, T_c, P_c = self.get_w_Tc_Pc(cas_i, cas_j)
        T_r = T / T_c
        return self.R*T_c/P_c*(self.d_B0_d_Tr_expr(T_r) + w*self.d_B1_d_Tr_expr(T_r))

    def d_dij_d_Tr(self, T):
        r"""
        .. math::
            \frac{\partial \delta_{ij}}{\partial T_r} = 2\frac{\partial B_{ij}}{\partial T_r}-\frac{\partial B_{ii}}{\partial T_r}-\frac{\partial B_{jj}}{\partial T_r}
            :label: eq_d_dij

        :param T: temperature [K]
        """
        return 2.*self.d_Bij_d_Tr(self.c1.cas_number, self.c2.cas_number, T) \
               - self.d_Bij_d_Tr(self.c1.cas_number, self.c1.cas_number, T) \
               - self.d_Bij_d_Tr(self.c2.cas_number, self.c2.cas_number, T)

    def bar_HiR_RT(self, cas_i, y_i, P, T):
        r"""Dimensionless residual partial molar enthalpy of component :math:`c1`

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
            \frac{\bar{H}^\mathrm{R}_i}{RT} = \frac{P}{RT_c}\left[
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
        return P/self.R/T_c*(self.d_Bij_d_Tr(cas_i, cas_i, T) + (1.-y_i)*(1.-y_i)*self.d_dij_d_Tr(T)) + self.ln_hat_phi_i_expr(cas_i, y_i, P, T)

    def bar_SiR_R(self, cas_i, y_i, P, T):
        r"""Dimensionless residual partial molar entropy of component :math:`c1`

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
            \frac{\bar{S}_i^\mathrm{R}}{R} = \frac{P}{RT_c}\left[
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
        return P/self.R/T_c*(
            self.d_Bij_d_Tr(cas_i, cas_i, T) + (1.-y_i)*(1.-y_i)*self.d_dij_d_Tr(T)
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
        return (self.B_ij_expr(cas_i, cas_i, T) + (1-y_i)*(1-y_i)*self.d_ij_expr(T))/self.R/T

    def X_R_dimensionless(self, method: callable, cas_i: str, y_i: float, P: float, T: float):
        """Residual property of :math:`X` for mixture.

        :param method: function to compute partial molar property of compound
        :type method: callable
        """
        y_j = 1.-y_i
        cas_j = self.other_cas[cas_i]
        return (
                y_i*method(cas_i, y_i, P, T) + (1.-y_i)*method(cas_j, y_j, P, T)
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

    def plot_residual_HSG(self, P, T, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_ylabel('Dimensionless')
            ax.set_xlabel('Gas-Phase Mole fraction %s' % self.c1.compound_name)

        y_1 = np.linspace(0., 1.)
        HR_RT = np.array(list(self.H_R_RT(self.c1.cas_number, i, P, T) for i in y_1))
        GR_RT = np.array(list(self.G_R_RT(self.c1.cas_number, i, P, T) for i in y_1))
        SR_R = np.array(list(self.S_R_R(self.c1.cas_number, i, P, T) for i in y_1))
        ax.plot(y_1, HR_RT, '-', label='$H^\\mathrm{R}/(RT)$')
        ax.plot(y_1, GR_RT, '--', label='$G^\\mathrm{R}/(RT)$')
        ax.plot(y_1, SR_R, '-.', label='$S^\\mathrm{R}/R$')
        ax.set_title('Residual properties at T = %5.2f K and P = %4.3e Pa' % (T, P))
        ax.legend()
        return ax
