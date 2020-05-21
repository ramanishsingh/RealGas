from scithermo.critical_constants import CriticalConstants
from scithermo.eos.cubic import PengRobinson
from scithermo.eos.virial import SecondVirial, BinarySecondVirial
import numpy as np
from scithermo.chem_constants import R_si_units
from shared.util import percent_difference

tol = 0.1
kwargs_i = {'T_c': 374.5, 'V_c': 9.44e-5, 'Z_c': 0.284, 'w': 0.0942, 'MW': 34.081,
                 'compound_name': 'H2S', 'cas_number': 'i'}
kwargs_i['P_c'] = kwargs_i['Z_c'] * R_si_units * kwargs_i['T_c'] / kwargs_i['V_c']
kwargs_j = {
    'T_c': 191.4, 'V_c': 0.100, 'Z_c': 0.286, 'w': 0.0115, 'MW': 16.042,
    'compound_name': 'CH4', 'cas_number': 'j'
}
kwargs_j['P_c'] = kwargs_j['Z_c'] * R_si_units * kwargs_j['T_c'] / kwargs_j['V_c']


def test_input_custom():
    """Test can input custom properties"""
    I = CriticalConstants(**kwargs_i)
    J = PengRobinson(**kwargs_i)
    K = SecondVirial(**kwargs_i)
    T, P = 300., 10e5
    assert percent_difference(J.iterate_to_solve_Z(T, P, 'vapor'), K.calc_Z_from_units(P, T)) < 0.1, "compressibility factors from second virial and PR differ by > 0.1 %%"


def test_input_binaryvirial():
    unary = SecondVirial(**kwargs_i)
    binary = BinarySecondVirial(i_kwargs=kwargs_i, j_kwargs=kwargs_j)
    P, T = 300., 10e5
    assert np.isclose(binary.calc_Z({'j': 0., 'i': 1.}, P, T), unary.calc_Z_from_units(P, T)), 'Binary compressibitily not same as unary compressibility'

