Definitions
===========

Residual properties for a given thermodynamic property :math:`M` are defined as

.. math::
    M = M^\text{IG} + M^\text{R}

where :math:`M^\text{IG}` is the value of the property in the ideal gas state
and :math:`M^\text{R}` is the residual value of the property.

More information on residual properties can be found in standard texts :cite:`Smith2005`

We **define** partial molar property :math:`\bar{M}_i` of species *i* in a mixture as

.. math::
    \bar{M}_i = \left(\frac{\partial(nM)}{\partial n_i}\right)_{P, T, n_j}
    :label: partial_molar_property

The mixture property is related to the partial molar property as

.. math::
    nM=\sum_i n_i\bar{M}_i

or, in terms of gas-phase mole fractions :math:`y_i`,

.. math::
    M=\sum_i y_i\bar{M}_i

The following relationships also hold

.. math::
    \bar{M}_i = \bar{M}^\text{IG} + \bar{M}^\text{R}
    :label: residual_partial_molar

.. math::
    M^\text{R} = \sum_i y_i\bar{M}^\text{R}
    :label: residual_molar

Nomenclature
------------

====  ===================        ==============================================================================
Code     Symbol                         Description
====  ===================        ==============================================================================
P      :math:`P`                 Pressure in Pa
V      :math:`V`                 Molar Volume in :math:`\text{m}^3/\text{mol}`
R      :math:`R`                 gas constant SI units (:math:`\text{m}^3\times\text{Pa}/\text{mol}/\text{K}`)
T      :math:`T`                 temperature in K
T_c    :math:`T_\text{c}`        critical temperature in K
P_c    :math:`P_\text{c}`        critical pressure in Pa
T_r    :math:`T_\text{r}`        reduced temperature (dimensionless)
V_c    :math:`V_\text{c}`        Critical volume :math:`\text{m}^3/\text{mol}`
w      :math:`\omega`            Accentric factor
y_i    :math:`y_i`               mole fraction of component *i*
--     :math:`n_i`               number of moles of component *i*
====  ===================        ==============================================================================
