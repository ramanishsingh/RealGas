import numpy as np
from src.gasthermo.eos.cubic import PengRobinson
from src.gasthermo.eos.virial import SecondVirial, SecondVirialMixture
from tests.test_cp_ig import compounds_to_test


def test_virial():
    T = 300.
    P = 1e5
    for i in compounds_to_test:
        I = SecondVirial(compound_name=i, pow=pow)
        assert np.isclose(
            I.G_R_RT_expr(P, T), I.H_R_RT_expr(P, T) - I.S_R_R_expr(P, T)
        ), 'Residual properties from virial not same'


def test_binary_virial():
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
    mixture = SecondVirialMixture(compound_names=['Water', 'Tetrahydrofuran'], k_ij=0.)
    import numpy as np
    y_1 = np.linspace(0., 1.)
    HR_RT = np.array(list(mixture.H_R_RT([i, 1.-i], P, T) for i in y_1))
    GR_RT = np.array(list(mixture.G_R_RT([i, 1.-i], P, T) for i in y_1))
    SR_R = np.array(list(mixture.S_R_R([i, 1.-i], P, T) for i in y_1))
    for g, h, s in zip(GR_RT, HR_RT, SR_R):
        assert np.isclose(g, h - s), 'residual not thermodynamically consistent'
    assert np.isclose(mixture.H_R_RT([1., 0.], P, T), I.H_R_RT_expr(P, T)), 'HR for binary virial does not converge to unary'
    assert np.isclose(mixture.G_R_RT([1., 0.], P, T), I.G_R_RT_expr(P, T)), 'GR for binary virial does not converge to unary'
    assert np.isclose(mixture.S_R_R([1., 0.], P, T), I.S_R_R_expr(P, T)), 'Sr for binary virial for component i does not converge to unary'
    assert np.isclose(mixture.H_R_RT([0., 1.], P, T), J.H_R_RT_expr(P, T)), 'Hr for binary virial of component j does not converge to virial'
    assert np.isclose(mixture.G_R_RT([0., 1.], P, T), J.G_R_RT_expr(P, T)), 'GR for binary virial of component j does not converge to unary'
    assert np.isclose(mixture.S_R_R([0., 1.], P, T), J.S_R_R_expr(P, T)), 'SR for binary virial of component j does not converge to unary'


def test_cubic():
    T = 300.
    P = 1e5
    for i in compounds_to_test:
        I = PengRobinson(compound_name=i)
        Z = I.iterate_to_solve_Z(T, P, 'vapor')
        V = Z*I.R*T/P
        assert np.isclose(I.G_R_RT_expr(P, V, T), I.H_R_RT_expr(P, V, T)-I.S_R_R_expr(P, V, T)), 'compound not thermodynamically consistend residual properties {}'.format(i)
