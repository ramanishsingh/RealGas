from chem_util.chem_constants import gas_constant as R
from chem_util.math import percent_difference

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
    :param T_c: critical temperature [K]
    :param P_c: critical pressure [Pa]
    :param V_c: critical molar volume [m^3/mol]
    :param Z_c: critical compressibility factor [dimensionless]
    :param w: accentric factor [dimensionless]
    :param tol: tolerance for percent difference in Zc calulcated and tabulated, set to 0.5
    :type tol: float, hard-coded
    """

    def __init__(self, dippr_no: str = None, compound_name: str = None, cas_number: str = None,
                 MW: float = None, P_c: float = None, V_c: float = None, Z_c: float = None,
                 T_c: float = None, w: float = None):
        file = os.path.join(ROOT_DIR, 'critical_constants.csv')
        my_header = [
            'Cmpd. no.', 'Name', 'Formula', 'CAS no.', 'Mol. wt. [g/mol]',
            'Tc [K]', 'Pc [MPa]', 'Vc [m3/kmol]', 'Zc', 'Acentric factor'
        ]
        self.dippr_no = dippr_no
        self.compound_name = compound_name
        self.cas_number = cas_number
        self.MW = MW
        self.P_c = P_c
        self.V_c = V_c
        self.Z_c = Z_c
        self.T_c = T_c
        self.w = w

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
            # make dummy class from DIPPR table, test that custom parameters are at least close to DIPPR
            dippr_cls = CriticalConstants(compound_name=self.compound_name)
            assert abs(percent_difference(self.Z_c, dippr_cls.Z_c)) < 50., 'Percent difference too high for Z_c'
            assert abs(percent_difference(self.V_c, dippr_cls.V_c)) < 50., 'Percent difference too high for V_c'
            assert abs(percent_difference(self.P_c, dippr_cls.P_c)) < 50., 'Percent difference too high for P_c'
            assert abs(percent_difference(self.T_c, dippr_cls.T_c)) < 50., 'Percent difference too high for T_c'
            assert abs(percent_difference(self.w, dippr_cls.w)) < 50., 'Percent difference too high for w'
            assert abs(percent_difference(self.MW, dippr_cls.MW)) < 0.01, 'Percent difference too high for MW'
            assert (self.compound_name is not None and self.cas_number is not None
                    and self.MW is not None and self.P_c is not None and self.V_c is not None and self.Z_c is not None
                    and self.w is not None and self.T_c is not None), 'Inconsistent input, need to input all values'

    def calc_Z_c(self):
        """Calculate critical compressibility, for comparison to tabulated value"""
        return self.P_c * self.V_c / R / self.T_c

    def Z_c_percent_difference(self):
        """calculate percent difference between Z_c calculated and tabulated"""
        return percent_difference(self.calc_Z_c(), self.Z_c)
