"""
.. todo::
    test that becomes constant with temperature when temperature gets small
"""


import unittest
from scithermo.cp_ig import CpIdealGas
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


if __name__ == '__main__':
    unittest.main()
