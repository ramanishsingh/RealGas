

class AntoineCoefficients:
    """
    all obtained from Perrys Chemical Engineering Handbook 7th edition

    .. note::
        all pressures in Pa
        all temperature in K
    """
    def __init__(self, name):
        self.p_sat = None
        if name == 'propane':
            self.C = [59.078, -3492.6, -6.0669, 1.0919e-5, 2]
            self.Trange = [85.47, 369.83] # K
        elif name == 'propene' or name == 'propylene':
            self.C = [57.263, -3382.4, -5.7707, 1.0431e-5, 2]
            self.Trange = [87.89, 365.57] # K
        elif name == 'water':
            self.p_sat = 21.068  # mm Hg at 23 deg C
            self.p_sat = self.p_sat / 750.061 * 14.5038 # psia
            self.Trange = [22.9+273.15, 23.1+273.15]
        elif name == 'ethane':
            self.C = [51.857, -2598.7, -5.1283, 1.4913e-5, 2]
            self.Trange = [90.35, 305.32] # K
        elif name == 'ethene' or name == 'ethylene':
            self.C = [74.242, -2707.2, -9.8462, 2.2457e-02, 1]
            self.Trange = [104, 282.34] # K
        else:
            raise NotImplemented
        self.name = name

    def calc_p_sat(self, T):
        # temperature in K
        if T < self.Trange[0] or T > self.Trange[1]:
            print('Temperature %f is outside range of T_min = %f to T_max = %f K'%(T, self.Trange[0], self.Trange[1]))
        C1, C2, C3, C4, C5 = self.C
        return math.exp(C1 + C2/T + C3*math.log(T) + C4*math.pow(T,C5))


import math


if __name__ == "__main__":
    for component in ['propane', 'propylene']:
        I = AntoineCoefficients(component)
        print(component, I.calc_p_sat(296.))
