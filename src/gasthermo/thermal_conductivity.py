import numpy as np
import matplotlib.pyplot as plt
from chem_util.math import percent_difference
from . import os, ROOT_DIR


class ThermalConductivity:
    r"""Thermal Conductivity of Inorganic and Organic Substances [W/m/K] :cite:`DIPPR`

    .. math::
        k = \frac{C_1T^{C_2}}{1 + C_3/T + C_4 / T^2}
        :label: k_vapor

    where :math:`k` is the thermal conductivity in W/m/K and :math:`T` is in K.
    Thermal conductivites are either at 1 atm or the vapor pressure, whichever is lower.


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
    :param C1: parameter in Equation :eq:`k_vapor`
    :type C1: float, derived from input
    :param C2: parameter in Equation :eq:`k_vapor`
    :type C2: float, derived from input
    :param C3: parameter in Equation :eq:`k_vapor`
    :type C3: float, derived from input
    :param C4: parameter in Equation :eq:`k_vapor`
    :type C4: float, derived from input
    :param units: units for :math:`k`, set to W/m/K
    :type units: str
    :param T_min_fit: minimum temperature for fit, defaults to T_min
    :param T_max_fit: maximum temperature for fit, defaults to T_max

    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 T_min_fit: float = None, T_max_fit: float = None, n_points_fit: int = 1000,
                 poly_order: int = 2):
        file = os.path.join(ROOT_DIR, 'thermal_conductivity.csv')
        my_header = [
            'Cmpd. no.', 'Name', 'Formula', 'CAS no.', 'Mol. wt. [g/mol]',
            'C1', 'C2', 'C3', 'C4',
            'Tmin [K]', 'Val at Tmin', 'Tmax [K]', 'Val at Tmax'
        ]
        self.dippr_no = dippr_no
        self.compound_name = compound_name
        self.cas_number = cas_number

        found_compound = False
        self.units = 'W/m/K'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            assert header == my_header, 'Wrong header!'
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == self.dippr_no or vals[1] == self.compound_name or vals[3] == self.cas_number:
                    assert not found_compound, 'Input compound found twice in table!'
                    found_compound = True
                    (self.dippr_no, self.compound_name, self.formula, self.cas_number) = vals[:4]
                    (self.MW, self.C1, self.C2, self.C3, self.C4,
                     self.T_min, self.K_T_min, self.T_max, self.K_T_max) = map(float, vals[4:])

        assert found_compound, 'No compound was found in table! for {}, {}, {}'.format(
            self.dippr_no, self.compound_name, self.cas_number
        )

        # test that values are consistent
        if percent_difference(self.K_T_min, self.eval(self.T_min)) > 0.1:
            raise Exception('Inconsistent data for {}, Percent difference is {}'.format(
            self.compound_name, percent_difference(self.K_T_min, self.eval(self.T_min)))
            )
        if percent_difference(self.K_T_max, self.eval(self.T_max)) > 0.1:
            raise Exception('Inconsistent data for {}, Percent difference is {}'.format(
                self.compound_name, percent_difference(self.K_T_max, self.eval(self.T_max)))
            )

        self.T_min_fit = T_min_fit
        self.T_max_fit = T_max_fit
        if self.T_min_fit is None:
            self.T_min_fit = self.T_min
        if self.T_max_fit is None:
            self.T_max_fit = self.T_max

        self.T_fit = np.linspace(self.T_min_fit, self.T_max_fit, n_points_fit)

        assert self.T_min <= self.T_min_fit < self.T_max, 'Minimum temp not in correct range: {} !<= {} !<= {}'.format(
            self.T_min, self.T_min_fit, self.T_max
        )
        assert self.T_min < self.T_max_fit <= self.T_max, 'Max temp not in correct range: {} !< {} !<= {}'.format(
            self.T_min, self.T_max_fit, self.T_max
        )

        self.poly_order = poly_order

        y = self.eval(self.T_fit)
        self.coeffs = np.polyfit(self.T_fit, y, self.poly_order)

        ssxm, ssxym, ssyxm, ssym = np.cov(self.T_fit, y, bias=1).flat
        self.R2 = ssxym*ssxym / (ssxm*ssym)
        assert self.R2 > 0.95, 'Poor fit to polynomial, R2 is only %3.2f' % self.R2

        self.K = np.poly1d(self.coeffs)

    def eval(self, T):
        """

        :param T: temperature in K
        :return: :math:`k` W/m/K (see equation :eq:`k_vapor`)
        """
        return self.C1*T**self.C2/(
            1. + self.C3/T + self.C4/T/T
        )

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        T_all = np.linspace(self.T_min, self.T_max, 1000)
        vals = self.eval(T_all)
        approx_vals = self.K(self.T_fit)
        ax.plot(T_all, vals, '-', label='DIPPR')
        ax.plot(self.T_fit, approx_vals, '--', markerfacecolor='None', label='Polynomial')
        ax.legend()
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('K [%s]' % self.units)


class ThermalConductivityMixture:
    """Viscosity of vapor mixture using Wilke mixing rule

    :param name_to_cas: mapping of chemical name to cas registry number
    :type name_to_cas: dict[:attr:`components`, str]
    :param mixing_rule: mixing rule for calculation of viscosity, defaults to :code:`Simple`
    :type mixing_rule: str, optional
    :param pure: pure component viscosity info, obtained rom  :class:`src.gasthermo.vapor_viscosity.Viscosity`
    :type pure: dict[:attr:`components`, Viscosity]
    """
    def __init__(self, name_to_cas: dict, mixing_rule='Simple'):
        self.pure = {
            key: ThermalConductivity(cas_number=val) for key, val in name_to_cas.items()
        }
        self.mixing_rule = mixing_rule

    def eval_simple(self, y_i, T):
        """Calculate thermal conductivity using simple relationship

        :param y_i: mole fraction of each component i
        :type y_i: dict[:attr:`component`, float]
        """
        return sum(y_i[i]*val.eval(T) for i, val in self.pure.items())

    def eval_HR(self, y_i, T):
        """Weights based off of sqrt of molecular weights"""
        return sum(y_i[i]*val.eval(T)*val.MW**(1/2) for i, val in self.pure.items()) / sum(
            y_i[i] * val.MW ** (1 / 2) for i, val in self.pure.items()
        )

    def eval(self, y_i, T):
        if self.mixing_rule == 'Simple':
            return self.eval_simple(y_i, T)
        if self.mixing_rule == 'Herning Zipperer':
            return self.eval_HR(y_i, T)

        raise Exception('Mixing rule {} not found'.format(self.mixing_rule))
