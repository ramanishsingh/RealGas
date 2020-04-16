import unittest
from scithermo.critical_constants import CriticalConstants
from scithermo.eos.cubic import PengRobinson
from scithermo.eos.virial import SecondVirial, BinarySecondVirial
from scithermo.chem_constants import R_si_units


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 0.5
        self.kwargs_i = {'T_c': 374.5, 'V_c': 9.44e-5, 'Z_c': 0.284, 'w': 0.0942, 'MW': 34.081,
                  'compound_name': 'H2S', 'cas_number': 'i'}
        self.kwargs_i['P_c'] = self.kwargs_i['Z_c']*R_si_units*self.kwargs_i['T_c']/self.kwargs_i['V_c']
        self.kwargs_j = {
            'T_c': 191.4, 'V_c': 0.100, 'Z_c': 0.286, 'w': 0.0115, 'MW': 16.042,
            'compound_name': 'CH4', 'cas_number': 'j'
        }
        self.kwargs_j['P_c'] = self.kwargs_j['Z_c']*R_si_units*self.kwargs_j['T_c']/self.kwargs_j['V_c']

    def test_input_custom(self):
        """Test can input custom properties"""
        I = CriticalConstants(**self.kwargs_i)
        J = PengRobinson(**self.kwargs_i)
        K = SecondVirial(**self.kwargs_i)
        T, P = 300., 10e5
        self.assertAlmostEqual(J.iterate_to_solve_Z(T, P, 'vapor'), K.calc_Z_from_units(P, T),places=2)

    def test_input_binaryvirial(self):
        unary = SecondVirial(**self.kwargs_i)
        binary = BinarySecondVirial(i_kwargs=self.kwargs_i, j_kwargs=self.kwargs_j)
        P, T = 300., 10e5
        self.assertAlmostEqual(binary.calc_Z({'j': 0., 'i': 1.}, P, T), unary.calc_Z_from_units(P, T))


if __name__ == '__main__':
    unittest.main()
