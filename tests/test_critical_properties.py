import unittest
from scithermo.critical_constants import CriticalConstants
from scithermo.eos.cubic import PengRobinson
from scithermo.eos.virial import SecondVirial
from scithermo.chem_constants import R_si_units

from tests.test_cp_ig import compounds_to_test


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 0.5

    def test_input_custom(self):
        """Test can input custom properties"""
        kwargs = {'T_c': 374.5, 'V_c': 9.44e-5, 'Z_c': 0.284, 'w': 0.0942, 'MW': 34.081,
                  'compound_name': 'H2S', 'cas_number': '0'}
        kwargs['P_c'] = kwargs['Z_c']*R_si_units*kwargs['T_c']/kwargs['V_c']
        I = CriticalConstants(**kwargs)
        J = PengRobinson(**kwargs)
        K = SecondVirial(**kwargs)
        T, P = 300., 10e5
        self.assertAlmostEqual(J.iterate_to_solve_Z(T, P, 'vapor'), K.calc_Z_from_units(P, T),places=2)


if __name__ == '__main__':
    unittest.main()
