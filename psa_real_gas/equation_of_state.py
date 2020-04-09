from .chem_constants import R_si_units


class PengRobinsonUnary:
    """
    :param T_c: Critical temperature [K]
    :type T_c: float
    :param P_c: Critical pressure [Pa]
    :type P_c: float
    :param w: Accentric factor [dimensionless]
    :type w: float
    """
    def __init__(self, **kwargs):
        self.R = 8.314
        self.T_c = kwargs.pop('T_c')
        self.P_c = kwargs.pop('P_c')
        self.w = kwargs.pop('w')

    def kappa_rule(self):
        """
        :return: expression for :math:`\\kappa`
        :rtype: float
        """
        return 0.37464 + 1.54226*self.w - 0.26992*self.w*self.w

    def a_expr(self, T):
        """
        .. math::
            \\left(0.45724\\right)

        :param T: temperature [K]
        :return: :math:`a`
        """

        Tic = self.T_i_c[component]
        Pic = self.P_i_c[component]
        R = self.R
        T_reduced = T / Tic
        ki = self.kappa_i[component]
        return 0.45724*self.R*self.R*Tic*Tic/Pic*pow(1. + ki - ki*pow(T_reduced, 0.5), 2)

class PengRobinson:
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def add_parameters(self, T_i_c=None, P_i_c=None, w_i=None, d_ij=None, **kwargs):
        """

        :param T_i_c: Critical temperatures [K]
        :type T_i_c: dict
        :param P_i_c: Critical pressures [Pa]
        :type P_i_c: dict
        :param w_i: Accentric factor, :math:`\omega_i` [dimensionless]
        :type w_i: dict
        :param kwargs: dictionary for other parameters
        """
        self.name += '_peng_robinson'

        # additional input parameters for PR
        self.T_i_c = T_i_c
        self.P_i_c = P_i_c
        self.w_i = w_i
        self.d_ij = d_ij

        self.input_parameters += ['T_i_c', 'P_i_c', 'w_i', 'd_ij']

        # additional dimensionless parameters
        self.kappa_i = {
            key: 0.37464 + 1.54226*val - 0.26992*val*val for key, val in self.w_i.items()
        }
        self.a_i = {
            key: self.a_i_rule(key) for key in self.components
        }
        self.b_i = {
            key: 0.07780*self.R*self.T_i_c[key]/self.P_i_c[key] for key in self.components
        }
        self.dimensionless_parameters += ['kappa_i', 'a_i', 'b_i']

        import math
        # todo: move this to __init__ method?
        self._sqrt_2 = math.sqrt(2)

    def a_i_rule(self, component):

    def a_ij_expr(self, i, j):
        """Mixing rule for :math:`a_{ij}`

        .. math::

            a_{ij} = \sqrt{a_ia_j}(1-\\delta_{ij})

        :param i: component name
        :param j: component name
        :return: mixing rule
        :rtype: float
        """
        return pyo.sqrt(self.a_i[i]*self.a_i[j])*(1. - self.d_ij[i, j])

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

    def retentate_gas_law_rule(self, *args):
        if self.on_inlet_retentate_boundary(*args):
            return self.enforce_P_feed_boundary(*args)

        Y_i_k = self.get_Y_i_r(*args)
        P = self.P_r_expr(*args)
        Z = self.Z_r_expr(*args)
        A = self.A_expr(Y_i_k, P)
        B = self.B_expr(Y_i_k, P)
        return self.residual(Z, A, B)

    def permeate_gas_law_rule(self, *args):
        if self.on_top_permeate_boundary(*args):
            return self.enforce_P_top_permeate_boundary(*args)

        Y_i_k = self.get_Y_i_p(*args)
        P = self.P_p_expr(*args)
        Z = self.Z_p_expr(*args)
        A = self.A_expr(Y_i_k, P)
        B = self.B_expr(Y_i_k, P)
        return self.residual(Z, A, B)


class PR_EOS(PengRobinson, Model):
    def __init__(self, **kwargs):
        Model.__init__(self, **kwargs)

    def add_parameters(self, T_i_c=None, P_i_c=None, w_i=None, **kwargs):
        Model.add_parameters(self, **kwargs)
        PengRobinson.add_parameters(self, T_i_c=T_i_c, P_i_c=P_i_c, w_i=w_i, **kwargs)

    def add_variables(self):
        self.Z_feed = pyo.Var(initialize=1, doc='compressibility factor')

    def add_equations(self):
        self.equation_of_state = pyo.Constraint(expr=self.equation_of_state_rule())

    def equation_of_state_rule(self):
        A = self.A_expr(self.Y_i_feed, self.P_feed)
        B = self.B_expr(self.Y_i_feed, self.P_feed)
        return self.residual(self.Z_feed, A, B)
