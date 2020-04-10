import matplotlib.pyplot as plt


def Cp_ig_demo():
    from thermodynamic_properties.cp_ig import CpIdealGas
    I = CpIdealGas(compound_name='Water', T_min_fit=250., T_max_fit=1000., poly_order=3)
    I.plot()
    plt.show()


def EOS_demo():
    from thermodynamic_properties.equation_of_state import PengRobinson, RedlichKwong, SoaveRedlichKwong
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$P$ [Pa]')
    ax.set_ylabel('$Z$')
    for cmpd, marker in zip(['Propane', 'Propylene'], ['o', 'x']):
        for cls, name in zip([PengRobinson, RedlichKwong, SoaveRedlichKwong], ['PR', 'RK', 'SRK']):
            I = cls(compound_name=cmpd)
            I.plot_Z_vs_P(300., 1e3, I.P_c, ax=ax, marker=marker, label='%s, %s' % (cmpd, name))
    ax.legend()
    plt.show()


if __name__ == '__main__':
    EOS_demo()
    Cp_ig_demo()