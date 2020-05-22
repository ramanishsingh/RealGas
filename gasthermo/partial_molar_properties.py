r"""
Partial Molar Properties
========================

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


Ideal Gas
---------

The *Gibbs theorem* is

    A partial molar property (other than volume) of a contituent species
    in an ideal-gas mixture is equal to the corresponding molar property
    of the species as a pure ideal gas at the mixture temperature but
    at a pressure equal to its partial pressure in the mixture


The partial molar volume of an ideal gas, :math:`\bar{V}_i^\text{IG}`, is

.. math::
    \bar{V}_i^\text{IG} = \frac{RT}{P}
    :label: barViIG

The partial molar enthalpy of an ideal gas, :math:`\bar{H}_i^\text{IG}`, is

.. math::
    \bar{H}_i^\text{IG} = H_i^\text{IG}
    :label: barHiIG

which results from the enthalpy of an ideal gas being independent of pressure.
Therefore, we can compute the ideal gas partial molar enthalpy if we have
the ideal gas heat capacities, as follows

.. math::
    \bar{H}_i^\text{IG} = H_i^\text{IG} = \int_{T_\text{ref}}^T C_{\text{p},i}^\text{IG}\mathrm{d}T^\prime

where :math:`T_\text{ref}` is a reference temperature and :math:`T^\prime` is a dummy variable for integration.
Often times, we want to compute these quantities in dimensionless units,

.. math::
    \begin{align}
        \bar{H}_i^{\text{IG},\star} &= \frac{\bar{H}_i^\text{IG}}{RT_\text{ref}} \\
                                    &= \int_{1}^{T^\star} \frac{C_{\text{p},i}^\text{IG}}{R}\mathrm{d}T^\prime \\
                                    &= \int_{1}^{T^\star} C_{\text{p},i}^{\text{IG},\star}\mathrm{d}T^\prime
    \end{align}

where

.. math::
    T^\star=T/T_\text{ref}

is a dimensionless variable
and

.. math::
    C_{\text{p},i}^{\text{IG},\star} = C_{\text{p},i}^\text{IG}/R

is a dimensionless parameter that is a function of :math:`T^\star`.

We note that the thermodynamic integration reference temperature does not have to be the same
as the temperature for scaling, but we have made them the same here for simplicity.

.. todo::
    Implement ideal gas partial molar entropy??
    Implement ideal gas partial molar free energy??
    This might not be useful though because these values seem to *depend* on mixture properties

Residual
--------

The residual partial molar volume of component *i*, :math:`\bar{V}_i^\text{R}`,
can be calculated as

.. math::
    \bar{V}_i^\mathrm{R} = RT\left(\frac{\partial \ln\hat{\phi}_i}{\partial P}\right)_{T,y}\\
    :label: barViR_definition

Defining the dimensionless quantity

.. math::
    \bar{V}_i^{\mathrm{R},\star} = \frac{\bar{V}_i^\text{R}}{RT_\text{ref}}

The expression in dimensionless units can be simplified to

.. math::
    \bar{V}_i^{\mathrm{R},\star} = T^\star\left(\frac{\partial \ln{\hat{\phi}_i}}{\partial P}\right)_{T^\star, y}
    :label: barViRstar_definition

The residual partial molar enthalpy of component *i*, :math:`\bar{H}_i^\text{R}`, can be calculated as

.. math::
    \bar{H}_i^{\text{R}} = -RT^2\left(\frac{\partial \ln\hat{\phi}_i}{\partial T}\right)
    :label: barHiR_definition

Defining the dimensionless quantity

.. math::
    \begin{align}
        \bar{H}_i^{\text{R},\star} &= \frac{\bar{H}_i^\text{R}}{RT_\text{ref}} \\
         &= -\frac{RT^2}{RT_\text{ref}}\left(\frac{\partial \ln\hat{\phi}_i}{\partial T}\right) \\
         &= -\frac{RT^2}{RT_\text{ref}^2}\left(\frac{\partial \ln\hat{\phi}_i}{\partial T^\star}\right) \\
    \end{align}

The expression in dimensionless units is computed as

.. math::
    \bar{H}_i^{\text{R},\star} = -\left(T^\star\right)^2\left(\frac{\partial \ln\hat{\phi}_i}{\partial T^\star}\right) \\
    :label: barHiRstar_definition

The residual partial molar free energy of component *i*, :math:`\bar{G}_i^\text{R}`,
which *defines* the fugacity coefficient :math:`\hat{\phi}_i`, is

.. math::
    \bar{G}_i^\text{R} = RT \ln\hat{\phi}_i
    :label: barGiR_definition


With these definitions, however, we note that we need an equation of state to calculate the partial molar properties.
In this package, the second-order virial equation of state currently implements the necessary derivatives.

"""


from gasthermo.eos.virial import SecondVirialMixture
import gasthermo.input as inp
from gasthermo.cp import CpIdealGas, CpStar
import numpy as np
import typing


class IdealMixture(inp.IdealMixture):
    def __init__(self, **kwargs):
        inp.IdealMixture.__init__(self, **kwargs)


class Mixture(SecondVirialMixture):
    r"""Mixture

    .. note::
        can only input both custom critical properties or both from DIPPR--cant have mixed at the moment

    :param ideal_gas: whether or not the gas is ideal, defaults to True
    :type ideal_gas: bool, optional
    :param kwargs: key-word arguments for instantiation of :class:`gasthermo.eos.virial.SecondVirialMixture`

    """

    def __init__(self, ideal_gas=True, **kwargs):
        self.ideal_gas = ideal_gas
        SecondVirialMixture.__init__(self, **kwargs)

