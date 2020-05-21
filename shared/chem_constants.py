"""
Chemical constants
------------------

    .. note::
        All units are in SI units

    The gas constant, :math:`R=8.314` m^3*Pa/mol/K
    Room temperature is 300~K

    :math:`N_\text{av}=6.022\\times 10^{23}`
"""

gas_constant = 8.314  # m**3*Pa/mol/K [=] m**3*kg/m/s/s/mol/K => m*m*kg/s/s/mol/K
room_temperature = 300.
N_av = 6.022e23
atmospheric_pressure = 101325   # Pa
k_B = gas_constant/N_av  # boltzmann constant, [m*m*kg/s/s/K / molec]
k_B_ergs = 1.3806e-16  # erg/molec/K


def SLPM__m3_s(Q_s, P_x, T_x=room_temperature, Z_x=1):
    r"""convert flow rate from L/s (at STP) to m**3/s at (*P*)

    If molar flow rate is held constant at temperature and pressure, the following is true

    .. math::
        n_\text{s} = n_\text{x}

    where s refers to standard conditions and x refers to any temperature and pressure.

    Then, using the molar volume of the gas, :math:`V`,
     and the volumetric flow rate of the gas, and substituting in the
    gas law

    .. math::
        \begin{align}
            \frac{Q_\text{s}}{V_\text{s}} &= \frac{Q_\text{x}}{V_\text{x}} \\
            \frac{Q_\text{s}P_\text{s}}{R T_\text{s}} &= \frac{Q_\text{x}P_\text{x}}{Z_\text{x}RT_\text{x}} \\
        \end{align}

    So that

    .. math::
        Q_\text{x} = Q_\text{s}\left(\frac{P_\text{s}}{P_\text{x}}\right)\left(\frac{T_\text{x}}{T_\text{s}}\right)Z_\text{x}

    :param Q_s: SLPM flow rate
    :param P_x: real pressure [Pa]
    :param T_x: real temperature [K]
    :param Z_x: real compressibility
    :return: Q_x, molar floar rate at P_x, T_x, Z_x in m**3/s
    """
    P_s = 101325.
    T_s = 273.15
    Q_s_SLPS = Q_s/60.  # standard liter per s
    Q_s_SM3PS = Q_s_SLPS / 1000.  # standard m**3 / s
    return Q_s_SM3PS * P_s / P_x * T_x / T_s * Z_x


Si = 28.0855
O = 15.99405


if __name__ == '__main__':
    area = 3.14*(3.71*3.71)/4/100./100.
    print(SLPM__m3_s(8.6, 10 * atmospheric_pressure) / area)
