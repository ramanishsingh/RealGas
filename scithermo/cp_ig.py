import numpy as np
import matplotlib.pyplot as plt
from scithermo.util import percent_difference


class CpIdealGas:
    """Heat Capacity :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` [J/mol/K] at Constant Pressure of Inorganic and Organic Compounds in the
    Ideal Gas State Fit to Hyperbolic Functions :cite:`DIPPR`

    .. math::
        C_{\\mathrm{p}}^{\\mathrm{IG}} = C_1 + C_2\\left[\\frac{C_3/T}{\\sinh{(C_3/T)}}\\right] + C_4 \\left[\\frac{C_5/T}{\\cosh{(C_5/T)}}\\right]^2
        :label: cp_ig


    where :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` is in J/mol/K and :math:`T` is in K.

    Computing integrals of Equation :eq:`cp_ig` is challenging.
    Instead, the function is fit to a polynomial within a range of interest (:attr:`T_min_fit`, :attr:`T_max_fit`)
    Then, this polynomial is integrated symbolically using :ref:`sympy` and lambdified into a function
    that can be called at runtime in an external computational environment by calling
    :meth:`scithermo.cp_ig.CpIdealGas.cp_ig_integral` (see Equation :eq:`cp_int`)


    :param dippr_no: dippr_no of compound by DIPPR table, defaults to None
    :type dippr_no: str, optional
    :param compound_name: name of chemical compound, defaults to None
    :type compound_name: str, optional
    :param cas_number: CAS registry number for chemical compound, defaults to None
    :type cas_number: str, optional
    :param MW: molecular weight in g/mol
    :type MW: float, derived from input
    :param Tmin: minimum temperature of validity for relationship [K]
    :type T_min: float, derived from input
    :param Tmax: maximum temperature of validity [K]
    :type T_max: float, derived from input
    :param C1: parameter in Equation :eq:`cp_ig`
    :type C1: float, derived from input
    :param C2: parameter in Equation :eq:`cp_ig`
    :type C2: float, derived from input
    :param C3: parameter in Equation :eq:`cp_ig`
    :type C3: float, derived from input
    :param C4: parameter in Equation :eq:`cp_ig`
    :type C4: float, derived from input
    :param C5: parameter in Equation :eq:`cp_ig`
    :type C5: float, derived from input
    :param units: units for :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}`, set to J/mol/K
    :type units: str



    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 T_min_fit: float = None, T_max_fit: float = None, n_points_fit: int = 1000,
                 poly_order: int=2):
        """
        :param n_points_fit: number of points for fitting polynomial and plotting, defaults to 1000
        :type n_points_fit: int, optional
        :param poly_order: order of polynomial for fitting,d efaults to 2
        :type poly_order: int, optional
        :param T_min_fit: minimum temperature for fitting, defaults to Tmin
        :param T_max_fit: maximum temperature for fitting, defaults to Tmax
        """
        from scithermo import os, ROOT_DIR
        file = os.path.join(ROOT_DIR, 'cp_ig.csv')
        my_header = [
            'Cmpd. no.', 'Name', 'Formula', 'CAS no.', 'Mol. wt. [g/mol]',
            'C1 x 1e-5', 'C2 x 1e-5', 'C3 x 1e-3', 'C4 x 1e-5', 'C5',
            'Tmin [K]', 'Cp at Tmin x 1e-5', 'Tmax [K]', 'Cp at Tmax x 1e-5'
        ]
        self.dippr_no = dippr_no
        self.compound_name = compound_name
        self.cas_number = cas_number

        found_compound = False
        self.units = 'J/mol/K'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            assert header == my_header, 'Wrong header!'
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == self.dippr_no or vals[1] == self.compound_name or vals[3] == self.cas_number:
                    assert not found_compound, 'Input compound found twice in table!'
                    found_compound = True
                    # found
                    (self.dippr_no, self.compound_name, self.formula, self.cas_number, self.MW,
                     self._C1em5, self._C2em5, self._C3em3, self._C4e5, self.C5,
                     self.T_min, self._Cp_Tmin_em5, self.T_max, self._Cp_Tmax_em5) = vals

        assert found_compound, 'No compound was found in table! for {}, {}, {}'.format(self.dippr_no, self.compound_name, self.cas_number)
        self.MW = float(self.MW)
        self.T_min = float(self.T_min)
        self.T_max = float(self.T_max)
        self.C1 = float(self._C1em5) * 1e2
        self.C2 = float(self._C2em5) * 1e2
        self.C3 = float(self._C3em3) * 1e3
        self.C4 = float(self._C4e5) * 1e2
        self.C5 = float(self.C5)
        self.Cp_Tmin = float(self._Cp_Tmin_em5) * 1e2
        self.Cp_Tmax = float(self._Cp_Tmax_em5) * 1e2

        if T_min_fit is None:
            T_min_fit = self.T_min
        if T_max_fit is None:
            T_max_fit = self.T_max

        self.T_fit = np.linspace(T_min_fit, T_max_fit, n_points_fit)
        self.T_min_fit = T_min_fit
        self.T_max_fit = T_max_fit

        assert self.T_min <= T_min_fit < self.T_max, 'Minimum temp not in correct range: {} !<= {} !<= {}'.format(self.T_min, T_min_fit, self.T_max)
        assert self.T_min < T_max_fit <= self.T_max, 'Max temp not in correct range: {} !< {} !<= {}'.format(self.T_min, T_max_fit, self.T_max)

        self.Cp_poly = np.poly1d(
            np.polyfit(self.T_fit, self.eval(self.T_fit), poly_order)
        )
        self.anti_derivative = np.polyint(self.Cp_poly)

        if self.get_numerical_percent_difference() > 0.001:
            self.plot()
            plt.show()
            raise Exception('Error in integration is too large! Try using a smaller temperature range for fitting '
                            'and/or increasing the number of fitting points and polynomial degree')

    def eval(self, T, f_sinh=np.sinh, f_cosh=np.cosh):
        """

        :param T: temperature in K
        :param f_sinh: function for hyperbolic sine, defaults to :ref:`np.sinh`
        :type f_sinh: callable
        :param f_cosh: function for hyperbolic cosine, defaults to :ref:`np.cosh`
        :type f_cosh: callable
        :return: :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` J/mol/K (see equation :eq:`cp_ig`)
        """
        C3_T = self.C3 / T
        C5_T = self.C5 / T
        return self.C1 + self.C2 * (C3_T/f_sinh(C3_T)) + self.C4 * (C5_T/f_cosh(C5_T))*(C5_T/f_cosh(C5_T))

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        T_all = np.linspace(self.T_min, self.T_max, 1000)
        vals = self.eval(T_all)
        approx_vals = self.Cp_poly(self.T_fit)
        ax.plot(T_all, vals, '-', label='Hyperbolic Functions')
        ax.plot(self.T_fit, approx_vals, '--', markerfacecolor='None', label='Polynomial')
        ax.legend()
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('CpIg [%s]' % self.units)

    def cp_ig_integral(self, T_a, T_b):
        """

        .. math::
            \\int_{T_a}^{T_b}C_{\\mathrm{p}}^{\\mathrm{IG}}(T^\\star) \\mathrm{d}T^\\star
            :label: cp_int

        :param T_a: start temperature in K
        :param T_b: finish temperature in K
        :return: integral
        """
        return self.anti_derivative(T_b) - self.anti_derivative(T_a)

    def numerical_integration(self, T_a, T_b):
        """Numerical integration using scipy"""
        import scipy.integrate as si
        return si.quad(self.eval, T_a, T_b)

    def get_numerical_percent_difference(self):
        """Calculate the percent difference with numerical integration obtained by :ref:`scipy`"""
        integral_poly_fit = self.cp_ig_integral(self.T_min_fit, self.T_max_fit)
        integral_numerical, err_numerical = self.numerical_integration(self.T_min_fit, self.T_max_fit)
        return percent_difference(integral_poly_fit, integral_numerical)
