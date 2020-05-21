from chem_util.math import percent_difference
from chem_util.chem_constants import gas_constant
from . import os, ROOT_DIR


class CriticalConstants:
    """
    Get critical constants of a compound

    If critical constants are not passed in, reads from DIPPR table

    :param dippr_no: dippr_no of compound by DIPPR table, defaults to None
    :type dippr_no: str, optional
    :param compound_name: name of chemical compound, defaults to None
    :type compound_name: str, optional
    :param cas_number: CAS registry number for chemical compound, defaults to None
    :type cas_number: str, optional
    :param MW: molecular weight in g/mol
    :type MW: float, derived from input
    :param T_c: critical temperature [K]
    :type T_c: float, derived from input
    :param P_c: critical pressure [Pa]
    :type P_c: float, derived from input
    :param V_c: critical molar volume [m^3/mol]
    :type V_c: float, derived from input
    :param Z_c: critical compressibility factor [dimensionless]
    :type Z_c: float, derived from input
    :param w: accentric factor [dimensionless]
    :type w: float, derived from input
    :param tol: tolerance for percent difference in Zc calulcated and tabulated, set to 0.5
    :type tol: float, hard-coded
    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None, **kwargs):
        file = os.path.join(ROOT_DIR, 'critical_constants.csv')
        my_header = [
            'Cmpd. no.', 'Name', 'Formula', 'CAS no.', 'Mol. wt. [g/mol]',
            'Tc [K]', 'Pc [MPa]', 'Vc [m3/kmol]', 'Zc', 'Acentric factor'
        ]
        self.R = gas_constant
        self.dippr_no = dippr_no
        self.compound_name = compound_name
        self.cas_number = cas_number

        self.MW = kwargs.pop('MW', None)
        self.P_c = kwargs.pop('P_c', None)
        self.V_c = kwargs.pop('V_c', None)
        self.Z_c = kwargs.pop('Z_c', None)
        self.T_c = kwargs.pop('T_c', None)
        self.w = kwargs.pop('w', None)

        self.tol = 0.5

        if self.MW is None and self.P_c is None and self.V_c is None and self.Z_c is None and self.T_c is None and self.w is None:
            # if havent input critical compounds, get from DIPPR table
            found_compound = False
            with open(file, 'r') as f:
                header = next(f).rstrip('\n').split(',')
                assert header == my_header, 'Wrong header!'
                for line in f:
                    vals = line.rstrip('\n').split(',')
                    if vals[0] == self.dippr_no or vals[1] == self.compound_name or vals[3] == self.cas_number:
                        assert not found_compound, 'Input compound found twice in table!'
                        found_compound = True
                        # found
                        (self.dippr_no, self.compound_name, self.formula, self.cas_number, *floating_point_vals) = vals
                        self.MW, self.T_c, self.P_c, self.V_c, self.Z_c, self.w = map(float, floating_point_vals)
                        self.P_c = 1e6 * self.P_c
                        self.V_c = self.V_c / 1000.

            assert found_compound, 'No compound was found in table! for {}, {}, {}'.format(self.dippr_no, self.compound_name, self.cas_number)
            assert self.Z_c_percent_difference() < self.tol, 'Critical compressibility inconsistency!'
        else:
            assert (self.compound_name is not None and self.cas_number is not None
                    and self.MW is not None and self.P_c is not None and self.V_c is not None and self.Z_c is not None
                    and self.w is not None and self.T_c is not None), 'Inconsistent input, need to input all values'

    def calc_Z_c(self):
        """Calculate critical compressibility, for comparison to tabulated value"""
        return self.P_c * self.V_c / self.R / self.T_c

    def Z_c_percent_difference(self):
        """calculate percent difference between Z_c calculated and tabulated"""
        return percent_difference(self.calc_Z_c(), self.Z_c)
