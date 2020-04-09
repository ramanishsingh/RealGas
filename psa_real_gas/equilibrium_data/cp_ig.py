

class CpIdealGas:
    """Heat Capacity :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` [J/kmol/K] at Constant Pressure of Inorganic and Organic Compounds in the
    Ideal Gas State Fit to Hyperbolic Functions :cite:`DIPPR`

    .. math::
        C_{\\mathrm{p}}^{\\mathrm{IG}} = C_1 + C_2\\left[\\frac{C_3/T}{\\sinh{(C_3/T)}\\right] + C_4 \\left[\\frac{C_5/T}{\\cosh{(C_5/T)}\\right]^2

    where :math:`C_{\\mathrm{p}}^{\\mathrm{IG}}` is in J/kmol/K and :math:`T` is in K.

    """
    def __init__(self, compound_name, verbose=False):
        from psa_real_gas import os, ROOT_DIR
        file = os.path.join(ROOT_DIR, 'equilibrium_data', 'heat_capacity_liquid.csv')
        self.num_constants = 5
        self.units = 'J/kmol/K'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == compound_name:
                    self.constants = {
                        key: float(vals[header.index('C%i' % key)]) for key in range(1, self.num_constants + 1)
                    }
                    self.value_T_min = float(vals[header.index('Val(Tmin)')])
                    self.T_min = float(vals[header.index('Tmin [K]')])
                    self.T_max = float(vals[header.index('Tmax [K]')])

        if verbose:
            print('Setting heat capacity constants for %s (taken from Perrys):' % compound_name)
            for key, val in self.constants.items():
                print('            %s:' % key, val)
            print('             value at 300K [J/kmol/K]=', self.eval(300.))

    def eval(self, T):
        """return heat capacity liquid in J/kmol/K"""
        return sum(self.constants[i]*pow(T, i-1) for i in range(1, self.num_constants + 1))

    def integral(self, T):
        return sum(self.constants[i]*pow(T, i)/i for i in range(1, self.num_constants + 1))

    def integral_dT(self, T_ref, T):
        """
        .. math::

            \\int_{Tref}^T CpL dT

        """
        return self.integral(T) - self.integral(T_ref)
