"""
.. todo::
    test that becomes constant with temperature when temperature gets small
"""


import logging
logging.basicConfig(level=logging.DEBUG)
import unittest
from scithermo.cp import CpIdealGas, CpStar, CpRawData, CpStarRawData
from scithermo.util import percent_difference

compounds_to_test = ['Butane', 'Carbon dioxide', 'Carbon monoxide',
                     'Propane', 'Propylene', 'Water']


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        """Set tolerance for percent difference"""
        self.kwargs = dict(T_min_fit=200., T_max_fit=500.)
        self.tol = 0.1

    def test_cp_ig_Tmin(self):

        for i in compounds_to_test:
            I = CpIdealGas(compound_name=i, **self.kwargs)
            self.assertTrue(
                percent_difference(I.Cp_Tmin, I.eval(I.T_min)) < self.tol
            )

    def test_cp_ig_Tmax(self):

        for i in compounds_to_test:
            I = CpIdealGas(compound_name=i, **self.kwargs)
            self.assertTrue(
                percent_difference(I.Cp_Tmax, I.eval(I.T_max)) < self.tol
            )

    def test_cp_star_Tmin(self):

        for i in compounds_to_test:
            I = CpStar(compound_name=i, T_ref=300., **self.kwargs)
            self.assertTrue(
                percent_difference(I.Cp_Tmin, I.eval(I.T_min)) < self.tol
            )

    def test_cp_star_Tmax(self):

        for i in compounds_to_test:
            I = CpStar(compound_name=i, T_ref=300., **self.kwargs)
            self.assertTrue(
                percent_difference(I.Cp_Tmax, I.eval(I.T_max)) < self.tol
            )

    def test_cp_star_ig(self):
        r"""
        test that the following is true

        .. math::
            \text{R}T_\text{ref}\int C_\mathrm{p}^\star \mathrm{d}T^{\star,\prime} = \int C_{\mathrm{p}}^{\text{IG}}\mathrm{d}T^\prime

        """
        for i in compounds_to_test:
            I = CpStar(compound_name=i, T_ref=300., **self.kwargs)
            J = CpIdealGas(compound_name=i, **self.kwargs)
            self.assertAlmostEqual(
                I.R*I.T_ref*I.cp_integral(200./I.T_ref, 400./I.T_ref), J.cp_integral(200., 400)
            )

    def test_cp_raw(self):
        """Test getting cp and fitting from file for units and dimensionless versions
        The code automatically fits Cp and calculates goodness of fit and will fail if fit is not good enough
        """
        import pandas as pd
        df = pd.read_csv('Cp_raw_data/MFI.tsv', delimiter='\t')
        I = CpRawData(df['T [K]'], df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)
        J = CpStarRawData(df['T [K]'], df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350., T_ref=200.)
        I = CpRawData(T_raw=df['T [K]'], Cp_raw=df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)
        J = CpStarRawData(T_ref=200., T_raw=df['T [K]'], Cp_raw=df['Cp [J/mol/K]'], T_min_fit=200., T_max_fit=350.)


if __name__ == '__main__':
    unittest.main()
