import math


class CpIdealGas:
    """Heat Capacity :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` [J/kmol/K] at Constant Pressure of Inorganic and Organic Compounds in the
    Ideal Gas State Fit to Hyperbolic Functions :cite:`DIPPR`

    .. math::
        C_{\\mathrm{p}}^{\\mathrm{IG}} = C_1 + C_2\\left[\\frac{C_3/T}{\\sinh{(C_3/T)}\\right] + C_4 \\left[\\frac{C_5/T}{\\cosh{(C_5/T)}\\right]^2
        :label: cp_ig

    where :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` is in J/kmol/K and :math:`T` is in K.

    :param dippr_no: dippr_no of compound by DIPPR table, defaults to None
    :type dippr_no: str, optional
    :param compound_name: name of chemical compound, defaults to None
    :type compound_name: str, optional
    :param cas_number: CAS registry number for chemical compound, defaults to None
    :type cas_number: str, optional
    :param MW: molecular weight in g/mol
    :type MW: float, derived from input
    :param Tmin: minimum temperature of validity for relationship [K]
    :type Tmin: float, derived from input
    :param Tmax: maximum temperature of validity [K]
    :type Tmax: float, derived from input
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
    :param units: units for :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}`, set to J/kmol/K
    :type units: str

    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None):
        from thermodynamic_properties import os, ROOT_DIR
        file = os.path.join(ROOT_DIR, 'cp_ig.csv')
        my_header = [
            'Cmpd. no.', 'Name', 'Formula', 'CAS no.', 'Mol. wt. [g/mol]',
            'C1 x 1e-5', 'C2 x 1e-5', 'C3 x 1e-3', 'C4 x 1e-5', 'C5',
            'Tmin [K]', 'Cp at Tmin x 1e-5', 'Tmax [K]', 'Cp at Tmax x 1e-5'
        ]
        self.dippr_no = dippr_no
        self.compund_name = compound_name
        self.cas_number = cas_number

        found_compound = False
        self.units = 'J/kmol/K'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            assert header == my_header, 'Wrong header!'
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == self.dippr_no or vals[1] == self.compund_name or vals[3] == self.cas_number:
                    assert not found_compound, 'Input compound found twice in table!'
                    found_compound = True
                    # found
                    (self.dippr_no, self.compund_name, self.formula, self.cas_number, self.MW,
                     self._C1em5, self._C2em5, self._C3em3, self._C4e5, self.C5,
                     self.Tmin, self._Cp_Tmin_em5, self.Tmax, self._Cp_Tmax_em5) = vals

        assert found_compound, 'No compound was found in table! for {}, {}, {}'.format(self.dippr_no, self.compund_name, self.cas_number)
        self.MW = float(self.MW)
        self.Tmin = float(self.Tmin)
        self.Tmax = float(self.Tmax)
        self.C1 = float(self._C1em5) * 1e5
        self.C2 = float(self._C2em5) * 1e5
        self.C3 = float(self._C3em3) * 1e3
        self.C4 = float(self._C4e5) * 1e5
        self.C5 = float(self.C5)
        self.Cp_Tmin = float(self._Cp_Tmin_em5) * 1e5
        self.Cp_Tmax = float(self._Cp_Tmax_em5) * 1e5

        # if verbose:
        #     print('Setting heat capacity constants for %s (taken from Perrys):' % compound_name)
        #     for key, val in self.constants.items():
        #         print('            %s:' % key, val)
        #     print('             value at 300K [J/kmol/K]=', self.eval(300.))

    def eval(self, T):
        """

        .. math::
            C_{\\mathrm{p}}^{\\mathrm{IG}} = C_1 + C_2\\left[\\frac{C_3/T}{\\sinh{(C_3/T)}\\right] + C_4 \\left[\\frac{C_5/T}{\\cosh{(C_5/T)}\\right]^2
            :label: cp_ig

        :param T: temperature in K
        :return: :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` J/kmol/K
        """
        C3_T = self.C3 / T
        C5_T = self.C5 / T
        return self.C1 + self.C2 * (C3_T/math.sinh(C3_T)) + self.C4 * (C5_T/math.cosh(C5_T))*(C5_T/math.cosh(C5_T))

    # def integral(self, T):
    #     return sum(self.constants[i] * pow(T, i) / i for i in range(1, self.num_constants + 1))
    #
    # def integral_dT(self, T_ref, T):
    #     """
    #     .. math::
    #
    #         \\int_{Tref}^T CpL dT
    #
    #     """
    #     return self.integral(T) - self.integral(T_ref)
