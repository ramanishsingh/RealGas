"""
.. todo::
    test that becomes constant with temperature when temperature gets small
"""


import numpy as np
import logging
logging.basicConfig(level=logging.DEBUG)
from src.GasThermo.cp import CpIdealGas, CpStar, CpRawData, CpStarRawData
from chem_util.math import percent_difference

compounds_to_test = ['Butane', 'Carbon dioxide', 'Carbon monoxide',
                     'Propane', 'Propylene', 'Water']

kwargs = dict(T_min_fit=200., T_max_fit=500.)
tol = 0.1


def test_cp_ig():
    for i in compounds_to_test:
        I = CpIdealGas(compound_name=i, **kwargs)
        assert percent_difference(I.Cp_Tmin, I.eval(I.T_min)) < tol, 'T min not consistent'
        assert percent_difference(I.Cp_Tmax, I.eval(I.T_max)) < tol, 'T max not consistent'


def test_cp_star():
    for i in compounds_to_test:
        I = CpStar(compound_name=i, T_ref=300., **kwargs)
        assert percent_difference(I.Cp_Tmin, I.eval(I.T_min)) < tol, 'T min not consistent'
        assert percent_difference(I.Cp_Tmax, I.eval(I.T_max)) < tol, 'T max not consistent'


def test_cp_star_ig():
    r"""
    test that the following is true
    .. math::
        \text{R}T_\text{ref}\int C_\mathrm{p}^\star \mathrm{d}T^{\star,\prime} = \int C_{\mathrm{p}}^{\text{IG}}\mathrm{d}T^\prime
    """
    for i in compounds_to_test:
        I = CpStar(compound_name=i, T_ref=300., **kwargs)
        J = CpIdealGas(compound_name=i, **kwargs)
        assert np.isclose(
            I.R*I.T_ref*I.cp_integral(200./I.T_ref, 400./I.T_ref), J.cp_integral(200., 400)
        ), 'Integral of dimensionless not the same as dimensional'


def test_cp_raw():
    """Test getting cp and fitting from file for units and dimensionless versions
    The code automatically fits Cp and calculates goodness of fit and will fail if fit is not good enough
    """
    from chem_util.io import read_csv
    df = read_csv('tests/Cp_raw_data/MFI.tsv', delimiter='\t')
    for key in ('T [K]', 'Cp [J/mol/K]'):
        df[key] = list(map(float, df[key]))

    I = CpRawData(df['T [K]'], df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)
    J = CpStarRawData(df['T [K]'], df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350., T_ref=200.)
    I = CpRawData(T_raw=df['T [K]'], Cp_raw=df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)
    J = CpStarRawData(T_ref=200., T_raw=df['T [K]'], Cp_raw=df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)
