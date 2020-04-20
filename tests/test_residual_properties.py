import unittest

from scithermo.eos.cubic import PengRobinson, RedlichKwong, SoaveRedlichKwong
from scithermo.eos.virial import SecondVirial, BinarySecondVirial
from .test_cp_ig import compounds_to_test


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        pass

    def test_virial(self):
        T = 300.
        P = 1e5
        for i in compounds_to_test:
            I = SecondVirial(compound_name=i, pow=pow)
            self.assertAlmostEqual(I.G_R_RT_expr(P, T), I.H_R_RT_expr(P, T)-I.S_R_R_expr(P, T))

    def test_binary_virial(self):
        r"""For binary virial, test that

        .. math::
            \frac{G^\text{R}}{RT} = \frac{H^\text{R}}{RT} - \frac{S^\text{R}}{R}

        and, in the limit of pure components, test that
        the residual properties calculated by BinarySecondVirial and SecondVirial (the pure comeponent one)
        are equivalent

        """
        P, T = 10e5, 250.
        I = SecondVirial(compound_name='Water')
        J = SecondVirial(compound_name='Tetrahydrofuran')
        mixture = BinarySecondVirial(
            i_kwargs={'compound_name': 'Water'},
            j_kwargs={'compound_name': 'Tetrahydrofuran'}, k_ij=0.)
        import numpy as np
        y_1 = np.linspace(0., 1.)
        HR_RT = np.array(list(mixture.H_R_RT(mixture.i.cas_number, i, P, T) for i in y_1))
        GR_RT = np.array(list(mixture.G_R_RT(mixture.i.cas_number, i, P, T) for i in y_1))
        SR_R = np.array(list(mixture.S_R_R(mixture.i.cas_number, i, P, T) for i in y_1))
        for g, h, s in zip(GR_RT, HR_RT, SR_R):
            self.assertAlmostEqual(g, h - s)
        self.assertAlmostEqual(mixture.H_R_RT(mixture.i.cas_number, 1., P, T), I.H_R_RT_expr(P, T))
        self.assertAlmostEqual(mixture.G_R_RT(mixture.i.cas_number, 1., P, T), I.G_R_RT_expr(P, T))
        self.assertAlmostEqual(mixture.S_R_R(mixture.i.cas_number, 1., P, T), I.S_R_R_expr(P, T))
        self.assertAlmostEqual(mixture.H_R_RT(mixture.j.cas_number, 1., P, T), J.H_R_RT_expr(P, T))
        self.assertAlmostEqual(mixture.G_R_RT(mixture.j.cas_number, 1., P, T), J.G_R_RT_expr(P, T))
        self.assertAlmostEqual(mixture.S_R_R(mixture.j.cas_number, 1., P, T), J.S_R_R_expr(P, T))

    def test_cubic(self):
        T = 300.
        P = 1e5
        for i in compounds_to_test:
            I = PengRobinson(compound_name=i)
            Z = I.iterate_to_solve_Z(T, P, 'vapor')
            V = Z*I.R*T/P
            self.assertAlmostEqual(I.G_R_RT_expr(P, V, T), I.H_R_RT_expr(P, V, T)-I.S_R_R_expr(P, V, T))


if __name__ == '__main__':
    unittest.main()
