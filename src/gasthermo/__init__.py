r"""


.. note::
    In this file, the extra :code:`src` directory is included.
    This import path is not necessary when :code:`gasthermo` is installed from :code:`pip`

Heat Capacity
-------------
Determine the temperature dependence of the ideal gas heat capacity, :math:`C_\text{p}^\text{IG}` for water

>>> import matplotlib.pyplot as plt
>>> from src.gasthermo.cp import CpIdealGas
>>> I = CpIdealGas(compound_name='Air', T_min_fit=250., T_max_fit=800., poly_order=3)
>>> I.eval(300.), I.Cp_units
(29.00369515161452, 'J/mol/K')
>>> I.eval(300.)/I.MW
1.0015088104839267
>>> # we can then plot and visualize the resutls
>>> fig, ax = I.plot()
>>> fig.savefig('docs/source/air.png')

And we will get something that looks like the following

.. image:: air.png

and we notice that the polynomial (orange dashed lines) fits the hyperbolic function well.


Equations of State
------------------

Cubic
*****

>>> from src.gasthermo.eos.cubic import PengRobinson, RedlichKwong, SoaveRedlichKwong
>>> P = 8e5  # Pa
>>> T = 300. # K
>>> cls = PengRobinson(compound_name='Propane')
>>> cls.iterate_to_solve_Z(P=P, T=T)
0.8803217148747003
>>> cls = RedlichKwong(compound_name='Propane')
>>> cls.iterate_to_solve_Z(P=P, T=T)
0.8896324535506871
>>> cls = SoaveRedlichKwong(compound_name='Propane')
>>> Z = cls.iterate_to_solve_Z(P=P, T=T)
>>> Z
0.8852911867637512
>>> # calculate residual properties
>>> from chem_util.chem_constants import gas_constant as R
>>> V = Z*R*T/P
>>> cls.S_R_R_expr(P, V, T)
-0.2739464327960456
>>> cls.H_R_RT_expr(P, V, T)
-0.4006043220021135
>>> cls.G_R_RT_expr(P, V, T) - cls.H_R_RT_expr(P, V, T) + cls.S_R_R_expr(P, V, T)
0.0

Virial
******

>>> from src.gasthermo.eos.virial import SecondVirial
>>> cls = SecondVirial(compound_name='Propane')
>>> cls.calc_Z_from_units(P=8e5, T=300.)
0.8726051032825523


Mixtures
--------
.. note::
    currently only implemented for virial equation of state

Residual Properties
*******************

Below, an example is shown for calculating residual properties of THF/Water mixtures

>>> from src.gasthermo.eos.virial import BinarySecondVirial
>>> P, T = 1e5, 300.
>>> mixture = BinarySecondVirial(i_kwargs = {'compound_name': 'Water'}, j_kwargs = {'compound_name': 'Tetrahydrofuran'}, k_ij=0.)
>>> import matplotlib.pyplot as plt
>>> fig, ax = mixture.plot_residual_HSG(P, T)
>>> fig.savefig('docs/source/THF-WATER.png')

So that the results look like the following

.. image:: THF-WATER.png

We note that the residual properties will not always vanish
in the limit of pure components like excess properties
since the pure-components may not be perfect gases.

Other Utilities
---------------
Determine whether a single real root of the cubic equation of state can be used for
simple computational implementation.
In some regimes, the cubic equation of state only has 1 real root--in this case, the compressibility
factor can be obtained easily.

>>> from src.gasthermo.eos.cubic import PengRobinson
>>> I = PengRobinson(compound_name='Propane')
>>> I.num_roots(300., 5e5)
3
>>> I.num_roots(100., 5e5)
1

Gotchas
-------
* All units are SI units

"""

import os
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
