from thermodynamic_properties.chem_constants import R_si_units
from thermodynamic_properties.critial_constants import CriticalConstants
from thermodynamic_properties.util import percent_difference
import numpy as np
from thermodynamic_properties.cp_ig import CpIdealGas
import matplotlib.pyplot as plt


class SecondVirial(CpIdealGas, CriticalConstants):
    """Virial equation of state for one component. See :cite:`Perry,Smith2005`

    :param R: gas constant, set to SI units
    :type R: float, hard-coded
    :param pow: function for computing power, defaults to :ref:`np.power`
    :type pow: callable, optional
    """
    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 pow: callable=np.power, **kwargs):
        """

        :param kwargs: used in :ref:`CpIdealGas`
        """
        CriticalConstants.__init__(self, dippr_no, compound_name, cas_number)
        CpIdealGas.__init__(self, dippr_no, compound_name, cas_number, **kwargs)
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

    def calc_Z(self, P, T):
        P_r = P/self.P_c
        T_r = T/self.T_c
        return 1 + (self.B0_expr(T_r) + self.w*self.B1_expr(T_r))*P_r/T_r

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
        Z = [self.calc_Z(pressure, T) for pressure in P]
        ax.plot(P, Z, symbol, markerfacecolor='None', **kwargs)
