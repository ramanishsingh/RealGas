import unittest
from thermodynamic_properties.critical_constants import CriticalConstants

from tests.test_cp_ig import compounds_to_test


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 0.5

    def test_Z_c(self):
        for i in compounds_to_test:
            I = CriticalConstants(compound_name=i)
            self.assertLess(I.Z_c_percent_difference(), self.tol)


if __name__ == '__main__':
    unittest.main()
