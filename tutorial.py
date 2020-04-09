import matplotlib.pyplot as plt


def EOS_demo():
    from thermodynamic_properties.equation_of_state import PengRobinson, RedlichKwong, SoaveRedlichKwong
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$P$ [MPa]')
    ax.set_ylabel('$Z$')
    for cmpd, marker in zip(['Propane', 'Propylene'], ['o', 'x']):
        for cls, name in zip([PengRobinson, RedlichKwong, SoaveRedlichKwong], ['PR', 'RK', 'SRK']):
            I = cls(compound_name=cmpd)
            I.plot_Z_vs_P(300., 1e3, I.P_c, ax=ax, marker=marker, label='%s, %s' % (cmpd, name))
    ax.legend()
    plt.show()


if __name__ == '__main__':
    EOS_demo()
