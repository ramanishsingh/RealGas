import matplotlib.pyplot as plt


def Cp_ig_demo():
    from thermodynamic_properties.cp_ig import CpIdealGas
    I = CpIdealGas(compound_name='Water', T_min_fit=250., T_max_fit=1000., poly_order=3)
    I.plot()
    plt.show()


def EOS_demo():
    from thermodynamic_properties.eos.cubic import PengRobinson, RedlichKwong, SoaveRedlichKwong
    from thermodynamic_properties.eos.virial import SecondVirial
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$P$ [Pa]')
    ax.set_ylabel('$Z$')
    for cmpd, marker in zip(['Propane', 'Propylene'], ['o', 'x']):
        for cls, name in zip([PengRobinson, RedlichKwong, SoaveRedlichKwong, SecondVirial], ['PR', 'RK', 'SRK', 'Virial']):
            I = cls(compound_name=cmpd, T_min_fit=200., T_max_fit=400.)
            I.plot_Z_vs_P(300., 1e3, I.P_c, ax=ax, marker=marker, label='%s, %s' % (cmpd, name))
    ax.legend()
    plt.show()


def EOS_roots():
    """Can we use a single root trick to speed up computations?"""
    from thermodynamic_properties.eos.cubic import PengRobinson
    I = PengRobinson(compound_name='Propane', T_min_fit=250., T_max_fit=500.)
    I.check_roots(250., 300., 1e5, 5e5)


if __name__ == '__main__':
    # EOS_roots()
    EOS_demo()
    # Cp_ig_demo()