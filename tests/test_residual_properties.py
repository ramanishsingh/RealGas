import unittest

from scithermo.eos.cubic import PengRobinson, RedlichKwong, SoaveRedlichKwong
from scithermo.eos.virial import SecondVirial
from .test_cp_ig import compounds_to_test


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.kwargs = dict(T_min_fit=200., T_max_fit=500.)

    def test_virial(self):
        T = 300.
        P = 1e5
        for i in compounds_to_test:
            I = SecondVirial(compound_name=i, pow=pow, **self.kwargs)
            self.assertAlmostEqual(I.G_R_RT_expr(P, T), I.H_R_RT_expr(P, T)-I.S_R_R_expr(P, T))

    def test_cubic(self):
        T = 300.

        P = 1e5
        for i in compounds_to_test:
            I = PengRobinson(compound_name=i, **self.kwargs)
            Z = I.iterate_to_solve_Z(T, P, 'vapor')
            V = Z*I.R*T/P
            self.assertAlmostEqual(I.G_R_RT_expr(P, V, T), I.H_R_RT_expr(P, V, T)-I.S_R_R_expr(P, V, T))


if __name__ == '__main__':
    unittest.main()
