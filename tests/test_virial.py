from src.gasthermo.eos import virial
import numpy as np


def test_correlation_coeffs():
    I = virial.Virial()
    import sympy as s
    T_r = s.symbols('T_r')
    B0 = I.B0_expr(T_r)
    dB0_dTr = s.diff(B0, T_r)
    B1 = I.B1_expr(T_r)
    dB1_dTr = s.diff(B1, T_r)
    print(I.d_B1_d_Tr_expr(T_r), dB1_dTr)
    assert str(I.d_B0_d_Tr_expr(T_r)) == str(dB0_dTr), 'B0 not correct!'
    assert str(I.d_B1_d_Tr_expr(T_r)) == str(dB1_dTr), 'B1 not correct!'
    assert np.isclose(I.B0_expr(1.), 0.083-0.422)
    assert np.isclose(I.B1_expr(1.), 0.139-0.172)



