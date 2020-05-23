r"""

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

which results from the enthalpy of an ideal gas being independent of pressure.
Therefore, we can compute the ideal gas partial molar enthalpy if we have
the ideal gas heat capacities, as follows

.. math::
    \bar{H}_i^\text{IG} = H_i^\text{IG} = \int_{T_\text{ref}}^T C_{\text{p},i}^\text{IG}\mathrm{d}T^\prime
    :label: barHiIG

where :math:`T_\text{ref}` is a reference temperature and :math:`T^\prime` is a dummy variable for integration.
Often times, we want to compute these quantities in dimensionless units,

.. math::
    \begin{align}
        \bar{H}_i^{\text{IG},\star} &= \frac{\bar{H}_i^\text{IG}}{RT_\text{ref}} \\
                                    &= \int_{1}^{T^\star} \frac{C_{\text{p},i}^\text{IG}}{R}\mathrm{d}T^\prime \\
                                    &=
    \end{align}

or

.. math::
    \bar{H}_i^{\text{IG},\star} = \int_{1}^{T^\star} C_{\text{p},i}^{\text{IG},\star}\mathrm{d}T^\prime
    :label: barHiIGstar

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


from realgas.eos.virial import SecondVirialMixture
import realgas.input as inp
from realgas.cp import CpIdealGas, CpStar
from chem_util.chem_constants import gas_constant as R
from realgas.input import IdealMixture
import numpy as np
import typing


class Mixture(SecondVirialMixture):
    """
    :param ideal: whether or not ideal gas, defaults to True
    :type ideal: bool, optional
    :param kwargs: key-word arguments for :class:`realgas.eos.virial.SecondVirialMixture`


    """

    def __init__(self, cp_args: typing.List[dict], ideal=True, **kwargs):
        self.ideal = ideal
        self.data = [
            CpIdealGas(**i) for i in cp_args
        ]
        self.cas_numbers = [I.cas_number for I in self.data]
        if not self.ideal:
            SecondVirialMixture.__init__(self, **kwargs)

    def bar_Hi_IG(self, cas_i: str, T, T_ref=0):
        """

        :param T_ref: reference temperature in K for enthalpy, defaults to 0
        :param T: temperautre in K
        :return: :math:`\bar{H}_i^\text{IG}`, see Equation :eq:`barHiIG`
        """
        assert cas_i in self.cas_numbers, 'Cas number not found!'
        i = self.cas_numbers.index(cas_i)
        return self.data[i].cp_integral(T_ref, T)

    def bar_Vi_IG(self, T, P):
        """

        :param T: temperature in K
        :param P: pressure in Pa
        :return: :math:`\bar{V}_i^\text{IG}`, see Equation :eq:`barViIG`
        """
        return R * T / P

    def bar_Hi_R(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        assert not self.ideal, 'No residual props for ideal gas'
        return R * T * self.bar_HiR_RT(cas_k, ys, P, T)

    def bar_Vi_R(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        assert not self.ideal, 'No residual props for ideal gas'
        return R * T * self.bar_ViR_RT(cas_k, ys, P, T)

    def bar_Hi(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T, T_ref=0):
        if self.ideal:
            return self.bar_Hi_IG(cas_k, T, T_ref=T_ref)

        return self.bar_Hi_IG(cas_k, T, T_ref=T_ref) + R * T * self.bar_HiR_RT(cas_k, ys, P, T)

    def bar_Vi(self, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        if self.ideal:
            return self.bar_Vi_IG(T, P)

        return self.bar_Vi_IG(T, P) + R * T * self.bar_Vi_R(cas_k, ys, P, T)

    def enthalpy(self, ys: typing.List[typing.Union[float, typing.Any]], P: float, T: float, T_ref=0.):
        """Residual property of :math:`X` for mixture.

        Similar to Equation :eq:`residual_molar` but in dimensionless form

        :param method: function to compute partial molar property of compound
        :type method: callable
        """
        return sum(
            ys[i] * self.bar_Hi(self.cas_numbers[i], ys, P, T, T_ref) for i in range(self.num_components)
        )


class MixtureDimensionless(SecondVirialMixture):
    """

    .. todo::
        add docs!

    :param ideal: whether or not ideal gas, defaults to True
    :type ideal: bool, optional
    :param kwargs: key-word arguments for :class:`realgas.eos.virial.SecondVirialMixture`
    """

    def __init__(self, cp_args: typing.List[dict], ideal=True, **kwargs):
        self.ideal = ideal
        self.data = [
            CpStar(**i) for i in cp_args
        ]
        self.cas_numbers = [I.cas_number for I in self.data]
        if not self.ideal:
            SecondVirialMixture.__init__(self, **kwargs)

    def bar_Hi_IG_star(self, cas_i, T_star, T_ref_star=0):
        """

        :param T_ref_star: dimensionless reference temperature in K for enthalpy, defaults to 0
        :param T_star: dimensionless temperature
        :return: :math:`\bar{H}_i^\text{IG}`, see Equation :eq:`barHiIGstar`
        """
        assert cas_i in self.cas_numbers, 'Cas number not found!'
        i = self.cas_numbers.index(cas_i)
        return self.data[i].cp_integral(T_star, T_ref_star)

    def bar_Hi_star(self, T_star, cas_k: str, ys: typing.List[typing.Union[float, typing.Any]], P, T):
        if self.ideal:
            return self.bar_Hi_IG_star(cas_k, T_star)

        return self.bar_Hi_IG_star(cas_k, T_star) + self.bar_HiR_star(T_star, cas_k, ys, P, T)

    def enthalpy_star(self, ys: typing.List[typing.Union[float, typing.Any]], P: float, T: float):
        """Residual property of :math:`X` for mixture.

        Similar to Equation :eq:`residual_molar` but in dimensionless form

        :param method: function to compute partial molar property of compound
        :type method: callable
        """
        return sum(
            ys[i] * self.bar_Hi_star(self.cas_numbers[i], ys, P, T) for i in range(self.num_components)
        )
