import unittest
from thermodynamic_properties.cp_ig import CpIdealGas


def percent_difference(x, y):
    return (x - y) / (x + y) * 200.


compounds_to_test = ['Butane', 'Carbon dioxide', 'Carbon monoxide', 'Decane', 'Hydrogen', 'Nitrogen', 'Oxygen',
                     'Pentane', 'Propane', 'Propylene', 'Water']


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        """Set tolerance for percent difference"""
        self.tol = 0.1

    def test_cp_ig_Tmin(self):
        for i in map(str, range(1, 14)):
            I = CpIdealGas(dippr_no=i)
            self.assertTrue(
                percent_difference(I.Cp_Tmin, I.eval(I.Tmin)) < self.tol
            )

        for i in compounds_to_test:
            I = CpIdealGas(compound_name=i)
            self.assertTrue(
                percent_difference(I.Cp_Tmin, I.eval(I.Tmin)) < self.tol
            )

    def test_cp_ig_Tmax(self):
        """Test Cp IG for random molecule"""
        for i in map(str, range(1, 14)):
            I = CpIdealGas(dippr_no=i)
            self.assertTrue(
                percent_difference(I.Cp_Tmax, I.eval(I.Tmax)) < self.tol
            )

        for i in compounds_to_test:
            I = CpIdealGas(compound_name=i)
            self.assertTrue(
                percent_difference(I.Cp_Tmax, I.eval(I.Tmax)) < self.tol
            )


if __name__ == '__main__':
    unittest.main()
