from thermodynamic_properties.chem_constants import R_si_units
import math
from typing import List
from thermodynamic_properties.critial_constants import CriticalConstants


class PengRobinson(CriticalConstants):
    """

    Peng-Robinson Equation of State :cite:`peng-robinson`

    :param R: gas constant, set to SI units
    :type R: float, hard-coded
    :param a: PR EOS Parameter, set in Equation :eq:`eq_a` from critical properties
    :type a: callable
    :param b: PR EOS Parameter, set in Equation :eq:`eq_b` from critical properties
    :type b: float
    """
    def __init__(self, **kwargs):
        CriticalConstants.__init__(self, **kwargs)
        self.R = R_si_units
        self.b = self.b_rule()
        self.a = lambda T: self.get_a_expr(T)

    def kappa_rule(self):
        """
        .. math::
            0.37464 + 1.54226\\omega - 0.26992\\omega^2
            :label: eq_kappa


        :return: expression for :math:`\\kappa`
        """
        return 0.37464 + 1.54226*self.w - 0.26992*self.w*self.w

    def get_a_expr(self, T, my_pow: callable=pow) -> callable:
        """
        .. math::
            a = \\left(0.45724\\frac{R^2T_c^2}{P_c}\\right)\\left[1 + \\kappa  - \\kappa\\left(\\frac{T}{T_c}\\right)^{1/2}\\right]^2
            :label: eq_a

        where :math:`\\kappa` is calculated in Equation :eq:`eq_kappa`

        .. note::
            This parameter depends on temperautre

        :param T: temperature [K]
        :param my_pow: function to perform power :math:`f(x,p) = x^p`, defaults to :py:ref:`pow`
        :type my_pow: callable
        """
        ki = self.kappa_rule()
        return 0.45724*self.R*self.R*self.T_c*self.T_c/self.P_c*my_pow(
            1. + ki - ki*my_pow(T/self.T_c, 0.5), 2
        )

    def b_rule(self):
        """
        .. math::
            b = 0.07780\\frac{RT_c}{P_c}
            :label: eq_b

        :return: :math:`b`
        """
        return 0.07780*self.R*self.T_c / self.P_c

    def residual(self, P, T, V):
        """Residual for PR-EOS (cubic equation of state)

        .. math::
            Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
            :label: eq_pr_eos_residual

        where

        .. math::
            \\begin{align}
                A &= \\frac{a(T)P}{R^2T^2}\\\\
                B &= \\frac{bP}{RT} \\\\
                Z &= \\frac{PV}{RT}
            \\end{align}

        :param P: pressure [Pa]
        :param T: temperature [K]
        :param V: molar volume [m^3/mol]

        """
        R = self.R
        A = self.a(T)*P/R/R/T/T
        B = self.b*P/R/T
        Z = P*V/R/T
        Z_squared = Z*Z
        B_squared = B*B

        return (
                Z * Z_squared
                - (1 - B) * Z_squared
                + Z * (A - 3. * B_squared - 2. * B)
                - (A * B - B_squared - B * B_squared) == 0.
        )



class PengRobinsonFactory:
    """

    .. note:: currently neglects the :math:`k_{ij}` mixing parameter

    :param data: Peng Robinson parameter class for each component
    :type data: dict[:attr:`components`, :py:class:`thermodynamic_properties.equation_of_state.PengRobinson`]
    :param components: names of components
    :type components: list
    """
    def __init__(self, components: List[str], args: List[dict]):
        """

        :param components: list of components
        :param args: args for each component
        """
        self.R = R_si_units
        self.unary = {
            key: PengRobinson(**val) for key, val in zip(components, args)
        }
        self.components = components
        self._sqrt_2 = self.sqrt(2)

    def sqrt(self, val):
        return math.sqrt(val)

    def a_ij_expr(self, i, j, T):
        """Mixing rule for :math:`a_{ij}`

        .. math::

            a_{ij} = \sqrt{a_ia_j}

        .. note::
            Assumes :math:`k_{ij}=0`

        :param i: component compound_name
        :param j: component compound_name
        :param T: temperature [K]
        :return: :math:`a_{ij}`
        """
        return self.sqrt(self.unary[i].a(T)*self.unary[j].a(T))

    def a_expr(self, Y_i_k):
        sum = 0.
        for i in Y_i_k.keys():
            for j in Y_i_k.keys():
                sum += Y_i_k[i] * Y_i_k[j] * self.a_ij_expr(i, j)
        return sum

    def b_expr(self, Y_i_k):
        return sum(
            self.b_i[key] * val for key, val in Y_i_k.items()
        )

    def A_expr(self, Y_i_k, P):
        return self.a_expr(Y_i_k) * P / self.R / self.R / self.T / self.T

    def B_expr(self, Y_i_k, P):
        return self.b_expr(Y_i_k) * P / self.R / self.T

    def C_k_expr(self, i, Y_i_k, P):
        """Constant in fugacity coefficient expression

        .. math::

            \\frac{A}{2\\sqrt{2}B}\\left(\\frac{b_i}{b}-\\frac{2\\sum_j Y_{j,k} a_{i,j}}{a}\\right)

        """
        A = self.A_expr(Y_i_k, P)
        B = self.B_expr(Y_i_k, P)
        return A/2./self._sqrt_2/B*(
            self.b_i[i]/self.b_expr(Y_i_k)
            - 2.*sum(Y_i_k[j]*self.a_ij_expr(i, j) for j in self.components)/self.a_expr(Y_i_k)
        )

    def phi_i_k_expr_brute_force(self, i, Z, Y_i_k, P):
        B = self.B_expr(Y_i_k, P)
        return pyo.exp(
            self.b_i[i]/self.b_expr(Y_i_k)*(Z - 1.)
            - pyo.log(Z - B)
            + self.C_k_expr(i, Y_i_k, P)*pyo.log(
                (Z + (1 + self._sqrt_2)*B) / (Z - (self._sqrt_2 - 1)*B)
            )
        )

    def phi_i_k_expr(self, i, Z, Y_i_k, P):
        B = self.B_expr(Y_i_k, P)
        return pyo.exp(
            self.b_i[i] / self.b_expr(Y_i_k) * (Z - 1.)
            + self.C_k_expr(i, Y_i_k, P)*pyo.log(
                (Z + (1 + self._sqrt_2) * B) / (Z - (self._sqrt_2 - 1) * B)
            )
        ) / (Z - B)

    def phi_i_r_expr(self, i, *args):
        Y_i_k = self.get_Y_i_r(*args)
        P = self.P_r_expr(*args)
        Z = self.Z_r_expr(*args)
        return self.phi_i_k_expr(i, Z, Y_i_k, P)

    def phi_i_p_expr(self, i, *args):
        Y_i_k = self.get_Y_i_p(*args)
        P = self.P_p_expr(*args)
        Z = self.Z_p_expr(*args)
        return self.phi_i_k_expr(i, Z, Y_i_k, P)

    @staticmethod
    def residual(Z, A, B):
        """Peng-Robinson EOS

        Ind. Eng. Chem. Fundam. Vol 15 1976
        """
        Z_squared = Z*Z
        B_squared = B*B

        return (
                 Z * Z_squared
                 - (1 - B) * Z_squared
                 + Z * (A - 3. * B_squared - 2. * B)
                 - (A * B - B_squared - B * B_squared) == 0.
        )


