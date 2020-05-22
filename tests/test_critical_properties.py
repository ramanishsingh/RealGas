from src.gasthermo.critical_constants import CriticalConstants
from src.gasthermo.eos.cubic import PengRobinson
from src.gasthermo.eos.virial import SecondVirial, SecondVirialMixture
import numpy as np
from chem_util.chem_constants import gas_constant
from chem_util.math import percent_difference

tol = 0.1
# H2S
i = 'H2S'
casi = '000'
T_ci = 374.5
V_ci = 9.44e-5
Z_ci = 0.284
w_i = 0.0942
MW_i = 34.081
P_ci = Z_ci * gas_constant * T_ci / V_ci
# compound_name_i = 'H2S'

# CH4
j = 'CH4'
casj = '2200'
T_cj = 191.4
V_cj = 0.100
Z_cj = 0.286
w_j = 0.0115
MW_j = 16.042
P_cj = Z_cj * gas_constant * T_cj / V_cj


def test_input_custom():
    """Test can input custom properties"""
    kwargs = dict(compound_name=i, cas_number=casi, T_c=T_ci, P_c=P_ci, V_c=V_ci, Z_c=Z_ci, w=w_i, MW=MW_i)
    I = CriticalConstants(**kwargs)
    J = PengRobinson(**kwargs)
    K = SecondVirial(**kwargs)
    T, P = 300., 10e5
    assert percent_difference(J.iterate_to_solve_Z(T, P, 'vapor'), K.calc_Z_from_units(P, T)) < 0.1, "compressibility factors from second virial and PR differ by > 0.1 %%"


def test_input_binaryvirial():
    unary = SecondVirial(compound_name=i, cas_number=casi, T_c=T_ci, P_c=P_ci, V_c=V_ci, Z_c=Z_ci, w=w_i, MW=MW_i)
    binary = SecondVirialMixture(
        compound_names=[i, j],
        cas_numbers=[casi, casj],
        T_cs=[T_ci, T_cj],
        P_cs=[P_ci, P_cj],
        V_cs=[V_ci, V_cj],
        Z_cs=[Z_ci, Z_cj],
        ws=[w_i, w_j],
        MWs=[MW_i, MW_j]
    )
    P, T = 300., 10e5
    assert np.isclose(binary.calc_Z_from_units([1., 0.], P, T), unary.calc_Z_from_units(P, T)), 'Binary compressibitily not same as unary compressibility'

