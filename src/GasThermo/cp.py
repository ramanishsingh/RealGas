import numpy as np
import logging
import matplotlib.pyplot as plt
from chem_util.math import percent_difference
from src.GasThermo.chem_constants import R_si_units


class CpIdealGas:
    r"""Heat Capacity :math:`C_{\mathrm{p}}^{\mathrm{IG}}` [J/mol/K] at Constant Pressure of Inorganic and Organic Compounds in the
    Ideal Gas State Fit to Hyperbolic Functions :cite:`DIPPR`

    .. math::
        C_{\mathrm{p}}^{\mathrm{IG}} = C_1 + C_2\left[\frac{C_3/T}{\sinh{(C_3/T)}}\right] + C_4 \left[\frac{C_5/T}{\cosh{(C_5/T)}}\right]^2
        :label: cp_ig

    where :math:`C_{\mathrm{p}}^{\mathrm{IG}}` is in J/mol/K and :math:`T` is in K.


    Computing integrals of Equation :eq:`cp_ig` is challenging.
    Instead, the function is fit to a polynomial within a range of interest
    so that it can be integrated by using an antiderivative that is a polynomial.


    :param dippr_no: dippr_no of compound by DIPPR table, defaults to None
    :type dippr_no: str, optional
    :param compound_name: name of chemical compound, defaults to None
    :type compound_name: str, optional
    :param cas_number: CAS registry number for chemical compound, defaults to None
    :type cas_number: str, optional
    :param MW: molecular weight in g/mol
    :type MW: float, derived from input
    :param T_min: minimum temperature of validity for relationship [K]
    :type T_min: float, derived from input
    :param T_max: maximum temperature of validity [K]
    :type T_max: float, derived from input
    :param T_min_fit: minimum temperature for fitting, defaults to Tmin
    :param T_max_fit: maximum temperature for fitting, defaults to Tmax
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
    :param Cp_units: units for :math:`C_{\mathrm{p}}^{\mathrm{IG}}`, defaults to J/mol/K
    :type Cp_units: str, optional
    :param T_units: units for :math:`T`, defaults to K
    :type T_units: str, optional
    :param n_points_fit: number of points for fitting polynomial and plotting, defaults to 1000
    :type n_points_fit: int, optional
    :param poly_order: order of polynomial for fitting, defaults to 2
    :type poly_order: int, optional

    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 T_min_fit: float = None, T_max_fit: float = None, n_points_fit: int = 1000,
                 poly_order: int = 2, T_units='K', Cp_units='J/mol/K'):
        from src.GasThermo import os, ROOT_DIR
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
        self.Cp_units = Cp_units
        self.T_units = T_units
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

        assert found_compound, 'No compound was found in table! for {}, {}, {}'.format(self.dippr_no,
                                                                                       self.compound_name,
                                                                                       self.cas_number)
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

        assert self.T_min <= T_min_fit < self.T_max, 'Minimum temp not in correct range: {} !<= {} !<= {}'.format(
            self.T_min, T_min_fit, self.T_max)
        assert self.T_min < T_max_fit <= self.T_max, 'Max temp not in correct range: {} !< {} !<= {}'.format(self.T_min,
                                                                                                             T_max_fit,
                                                                                                             self.T_max)
        self.poly_order = poly_order

        self.Cp_poly = np.poly1d(
            np.polyfit(self.T_fit, self.eval(self.T_fit), self.poly_order)
        )
        self.anti_derivative = np.polyint(self.Cp_poly)

        if self.get_numerical_percent_difference() > 0.001:
            fig, ax = self.plot()
            fig.savefig('error_cp.png')
            raise Exception('Error in integration is too large! Try using a smaller temperature range for fitting '
                            'and/or increasing the number of fitting points and polynomial degree')

    def eval(self, T, f_sinh=np.sinh, f_cosh=np.cosh):
        """

        :param T: temperature in K
        :param f_sinh: function for hyperbolic sine, defaults to :code:`np.sinh`
        :type f_sinh: callable
        :param f_cosh: function for hyperbolic cosine, defaults to :code:`np.cosh`
        :type f_cosh: callable
        :return: :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` J/mol/K (see equation :eq:`cp_ig`)
        """
        C3_T = self.C3 / T
        C5_T = self.C5 / T
        return self.C1 + self.C2 * (C3_T / f_sinh(C3_T)) + self.C4 * (C5_T / f_cosh(C5_T)) * (C5_T / f_cosh(C5_T))

    def plot(self, fig= None, ax=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        T_all = np.linspace(self.T_min, self.T_max, 1000)
        vals = self.eval(T_all)
        approx_vals = self.Cp_poly(self.T_fit)
        ax.plot(T_all, vals, '-', label='Hyperbolic Functions')
        ax.plot(self.T_fit, approx_vals, '--', markerfacecolor='None', label='Polynomial')
        ax.legend()
        ax.set_xlabel('Temperature [%s]' % self.T_units)
        ax.set_ylabel('CpIg [%s]' % self.Cp_units)
        return fig, ax

    def cp_integral(self, T_a, T_b):
        r"""

        .. math::
            \int_{T_a}^{T_b}C_{\mathrm{p}}^{\mathrm{IG}}(T^\prime) \mathrm{d}T^\prime
            :label: cp_int

        :param T_a: start temperature in K
        :param T_b: finish temperature in K
        :return: integral
        """
        return self.anti_derivative(T_b) - self.anti_derivative(T_a)

    def numerical_integration(self, T_a, T_b) -> tuple:
        """Numerical integration using scipy"""
        import scipy.integrate as si
        return si.quad(self.eval, T_a, T_b)

    def get_numerical_percent_difference(self):
        """Calculate the percent difference with numerical integration obtained by :ref:`scipy`"""
        integral_poly_fit = self.cp_integral(self.T_min_fit, self.T_max_fit)
        integral_numerical, err_numerical = self.numerical_integration(self.T_min_fit, self.T_max_fit)
        return percent_difference(integral_poly_fit, integral_numerical)


class CpStar(CpIdealGas):
    r"""Dimensionless Heat Capacity at Constant Pressure of Inorganic and Organic Compounds in the
    Ideal Gas State Fit to Hyperbolic Functions :cite:`DIPPR`

    The dimensionless form is obtained by introducing the following variables

    .. math::
        \begin{align}
            C_{\mathrm{p}}^\star &= \frac{C_\mathrm{p}^\mathrm{IG}}{\text{R}} \\
            T^\star &= \frac{T}{T_\text{ref}}
        \end{align}

    where R is the gas constant in units of J/mol/K,
    and :math:`T_\text{ref}` is a reference temperature [K] input by the user (see :attr:`T_ref`)

    The heat capacity in dimensionless form becomes

    .. math::
        C_{\mathrm{p}}^\star = C_1^\star + C_2^\star \left[\frac{C_3^\star/T^\star}{\sinh{(C_3^\star/T^\star)}}\right] + C_4^\star \left[\frac{C_5^\star/T^\star}{\cosh{(C_5^\star/T^\star)}}\right]^2
        :label: cp_ig_star

    where

    .. math::
        \begin{align}
            C_1^\star &= \frac{C_1}{\text{R}} \\
            C_2^\star &= \frac{C_2}{\text{R}} \\
            C_3^\star &= \frac{C_3}{T_\text{ref}} \\
            C_4^\star &= \frac{C_4}{\text{R}} \\
            C_5^\star &= \frac{C_5}{T_\text{ref}}
        \end{align}


    :param T_ref: reference temperature [K] for dimensionless computations
    :type T_ref: float

    """

    def __init__(self, T_ref: float, **kwargs):
        """
        :param kwargs: for CpIdealGas
        """
        CpIdealGas.__init__(self, Cp_units='dimensionless', T_units='dimensionless', **kwargs)
        self.T_ref = T_ref
        self.R = R_si_units

        self.T_min = self.T_min / self.T_ref
        self.T_max = self.T_max / self.T_ref
        self.C1 = self.C1 / self.R
        self.C2 = self.C2 / self.R
        self.C3 = self.C3 / self.T_ref
        self.C4 = self.C4 / self.R
        self.C5 = self.C5 / self.T_ref
        self.Cp_Tmin = self.Cp_Tmin / self.R
        self.Cp_Tmax = self.Cp_Tmin / self.R

        self.T_fit = self.T_fit / self.T_ref
        self.T_min_fit = self.T_min_fit / self.T_ref
        self.T_max_fit = self.T_max_fit / self.T_ref

        del self.Cp_poly
        del self.anti_derivative
        self.Cp_poly = np.poly1d(
            np.polyfit(self.T_fit, self.eval(self.T_fit), self.poly_order)
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
        :param f_sinh: function for hyperbolic sine, defaults to :code:`np.sinh`
        :type f_sinh: callable
        :param f_cosh: function for hyperbolic cosine, defaults to :code:`np.cosh`
        :type f_cosh: callable
        :return: :math:`C_{\\mathrm{p}}^{\\star}` [dimensionless] (see equation :eq:`cp_ig_star`)
        """
        return CpIdealGas.eval(self, T, f_sinh, f_cosh)


class CpRawData:
    """From raw data for Cp(T)
    * fit to polynomial of temperature
    * fit polynomial to antiderivative


    :param T_min_fit: minimum temperature for fitting function [K]
    :type T_min_fit: float, optional
    :param T_max_fit: maximum temperature for fitting function [K]
    :type T_max_fit: float, optional
    :param poly_order: order of polynomial for fitting, defaults to 2
    :type poly_order: int, optional
    :param T_raw: raw temperatures in K
    :type T_raw: list
    :param Cp_raw: raw heat capacities in J/K/mol
    :type Cp_raw: list
    :param Cp_units: units for :math:`C_{\mathrm{p}}`, defaults to J/mol/K
    :type Cp_units: str, optional
    :param T_units: units for :math:`T`, defaults to K
    :type T_units: str, optional

    """

    def __init__(self, T_raw: list, Cp_raw: list, T_min_fit: float = None, T_max_fit: float = None,
                 poly_order: int = 2, T_units='K', Cp_units='J/mol/K'):
        self.poly_order = poly_order
        self.T_min_fit = T_min_fit
        self.T_max_fit = T_max_fit
        self.T_raw = T_raw
        self.Cp_raw = Cp_raw
        self.Cp_units = Cp_units
        self.T_units = T_units

        assert len(T_raw) == len(Cp_raw), 'Inconsistent input data'

        if self.T_min_fit is None:
            self.T_min_fit = min(T_raw)
        if self.T_max_fit is None:
            self.T_max_fit = max(T_raw)

        indices_for_fit = [i for i in range(len(T_raw)) if self.T_min_fit <= T_raw[i] <= self.T_max_fit]
        assert len(indices_for_fit) > 1, 'Not enough indices to fit within temperature range of raw data({},{})'.format(min(T_raw), max(T_raw))
        self.T_fit = [T_raw[i] for i in indices_for_fit]
        self.Cp_fit = [Cp_raw[i] for i in indices_for_fit]

        self.Cp_poly = np.poly1d(
            np.polyfit(self.T_fit, self.Cp_fit, self.poly_order)
        )
        self.anti_derivative = np.polyint(self.Cp_poly)
        self.R2 = self.get_R2()
        logging.debug('R2 is {}'.format(self.R2))

        if not (0.99 <= self.R2 <= 1.):
            fig, ax = self.plot()
            fig.savefig('cp_error.png')
            plt.show()
            raise Exception('Fit is too poor (R2={} not in (0.99,1)) too large! Try using a smaller temperature range for fitting '
                            'and/or increasing the number of fitting points and polynomial degree'.format(self.R2))

    def eval(self, T):
        return self.Cp_poly(T)

    def get_R2(self):
        y_mean = np.mean(self.Cp_fit)
        SS_tot = sum((y_i - y_mean)*(y_i - y_mean) for y_i in self.Cp_fit)
        SS_res = sum((y_i - f_i)*(y_i-f_i) for y_i, f_i in zip(self.Cp_fit, self.Cp_poly(self.T_fit)))
        return 1. - SS_res/SS_tot

    def plot(self, fig=None, ax=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        approx_vals = self.Cp_poly(self.T_fit)
        ax.plot(self.T_raw, self.Cp_raw, 'o', markerfacecolor='None', label='Raw Data')
        ax.plot(self.T_fit, approx_vals, '--', markerfacecolor='None', label='Polynomial')
        ax.legend()
        ax.set_xlabel('Temperature [%s]' % self.T_units)
        ax.set_ylabel('Cp [%s]' % self.Cp_units)
        return fig, ax

    def get_max_percent_difference(self):
        """Get largest percent difference"""
        y_fit = self.Cp_poly(self.T_fit)
        pd = [percent_difference(i, j) for i, j in zip(y_fit, self.Cp_fit)]
        return max(pd)


class CpStarRawData(CpRawData):
    """From raw data for Cp

        * convert to dimensionless units
        * fit cp from CpRawData

    :param T_ref: reference temperature in K
    :type T_ref: float
    :param R: gas constant, defaults to R_si_units
    :type R: float, optional
    :param T_min_fit: minimum temperature for fitting function [K]
    :type T_min_fit: float, optional
    :param T_max_fit: maximum temperature for fitting function [K]
    :type T_max_fit: float, optional
    :param poly_order: order of polynomial for fitting, defaults to 2
    :type poly_order: int, optional
    :param T_raw: raw temperatures in K
    :type T_raw: list
    :param Cp_raw: raw heat capacities in J/K/mol
    :type Cp_raw: list

    """

    def __init__(self, T_raw: list, Cp_raw: list, T_ref: float, R: float=R_si_units,
                 T_min_fit: float = None, T_max_fit: float = None, **kwargs):
        T_raw = [i/T_ref for i in T_raw]
        Cp_raw = [i/R for i in Cp_raw]
        if T_min_fit is not None:
            T_min_fit = T_min_fit / T_ref

        if T_max_fit is not None:
            T_max_fit = T_max_fit / T_ref

        kwargs['Cp_units'] = 'dimensionless'
        kwargs['T_units'] = 'dimensionless'
        kwargs = {'T_min_fit': T_min_fit, 'T_max_fit': T_max_fit, **kwargs}

        CpRawData.__init__(self, T_raw, Cp_raw, **kwargs)
