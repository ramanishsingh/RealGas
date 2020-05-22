"""
Input
=====

Helper classes for inputting parameters
"""


class IdealMixture:
    r"""

    :param num_components: number of components
    :type num_components: int
    :param dippr_nos: dippr numbers of components
    :type dippr_nos: typing.Optional[typing.Union[str, None]]
    :param compound_names: names of components
    :type compound_names: typing.Optional[typing.Union[str, None]]
    :param cas_numbers: cas registry numbers
    :type cas_numbers:  typing.Optional[typing.Union[str, None]]
    """
    def __init__(self, **kwargs):
        self.num_components = kwargs.pop('num_components', 0)
        self.dippr_nos = kwargs.pop('dippr_nos', None)
        self.compound_names = kwargs.pop('compound_names', None)
        self.cas_numbers = kwargs.pop('cas_numbers', None)
        self._has_been_setup = False

    def setup(self):
        """setup input parameters"""
        self._find_number_components()
        self.check_input()
        self._update_params()
        self._has_been_setup = True

    def _find_number_components(self):
        if self.dippr_nos is not None:
            self.num_components = len(self.dippr_nos)
        if self.compound_names is not None:
            self.num_components = len(self.compound_names)
        if self.cas_numbers is not None:
            self.num_components = len(self.cas_numbers)

    def _update_params(self):
        if self.dippr_nos is None:
            self.dippr_nos = [None for I in range(self.num_components)]
        if self.compound_names is None:
            self.compound_names = [None for I in range(self.num_components)]
        if self.cas_numbers is None:
            self.cas_numbers = [None for I in range(self.num_components)]

    def check_input(self):
        assert self.num_components > 0, 'No components found!'

    def get_point_input(self, i: int) -> dict:
        """

        :param i: index of point
        :return: keyword arguments for point input
        """
        return {'dippr_no': self.dippr_nos[i] if self.dippr_nos is not None else None,
                'compound_name': self.compound_names[i] if self.compound_names is not None else None,
                'cas_number': self.cas_numbers[i] if self.cas_numbers is not None else None}

    def set_point_input(self, i: int, **kwargs):
        """

        :param i: index of point
        """
        self.dippr_nos[i] = kwargs.pop('dippr_no')
        self.compound_names[i] = kwargs.pop('compound_name')
        self.cas_numbers[i] = kwargs.pop('cas_number')

    def check_params(self):
        pass


class RealMixture(IdealMixture):
    r"""

    :param MWs: molecular weights of component sin g/mol
    :type MWs: typing.Optional[typing.List[float]]
    :param P_cs: critical pressures of componets in Pa
    :type P_cs: typing.Optional[typing.List[float]]
    :param T_cs: critical temperatures of pure components in K
    :type T_cs: typing.Optional[typing.List[float]]
    :param V_cs: critical molar volumes of pure components in m**3/mol
    :type V_cs: typing.Optional[typing.List[float]]
    :param ws: accentric factors of components
    :type ws: typing.Optional[typing.List[float]]
    :param k_ij: equation of state mixing rule in calculation of critical temperautre, see Equation :eq:`Tc_combine`. When :math:`i=j` and for chemical similar species, :math:`k_{ij}=0`. Otherwise, it is a small (usually) positive number evaluated from minimal :math:`PVT` data or, in the absence of data, set equal to zero.
    :type k_ij: typing.Union[float, typing.List[float]], defaults to 0

    """

    def __init__(self, **kwargs):
        IdealMixture.__init__(self, **kwargs)
        self.MWs = kwargs.pop('MWs', None)
        self.P_cs = kwargs.pop('P_cs', None)
        self.V_cs = kwargs.pop('V_cs', None)
        self.Z_cs = kwargs.pop('Z_cs', None)
        self.T_cs = kwargs.pop('T_cs', None)
        self.ws = kwargs.pop('ws', None)
        self.k_ij = kwargs.pop('k_ij', 0.)

    def _find_number_components(self):
        IdealMixture._find_number_components(self)
        if self.dippr_nos is None and self.compound_names is None and self.cas_numbers is None:
            assert (
                    self.MWs is not None
                    and self.P_cs is not None
                    and self.T_cs is not None
                    and self.V_cs is not None
                    and self.Z_cs is not None
                    and self.ws is not None

                    ), 'Incorrect input'
            assert self.num_components == 0, 'What happended?'
        if self.num_components == 0:
            self.num_components = len(self.MWs)

        if isinstance(self.k_ij, float):
            val = self.k_ij
            self.k_ij = [[0. for i in range(self.num_components)] for j in range(self.num_components)]
            for i in range(self.num_components):
                for j in range(self.num_components):
                    if i != j:
                        self.k_ij[i][j] = val

        for i in range(self.num_components):
            for j in range(self.num_components):
                assert abs(self.k_ij[i][j]) < 1e-8, 'K[i][i] must be zero!'

    def _update_params(self):
        IdealMixture._update_params(self)
        if self.MWs is None:
            self.MWs = [None for I in range(self.num_components)]
        if self.P_cs is None:
            self.P_cs = [None for I in range(self.num_components)]
        if self.Z_cs is None:
            self.Z_cs = [None for I in range(self.num_components)]
        if self.T_cs is None:
            self.T_cs = [None for I in range(self.num_components)]
        if self.V_cs is None:
            self.V_cs = [None for I in range(self.num_components)]
        if self.ws is None:
            self.ws = [None for I in range(self.num_components)]

    def get_point_input(self, i: int) -> dict:
        """

        :param i: index of point
        :return: keyword arguments for point input
        """
        data = IdealMixture.get_point_input(self, i)
        data['MW'] = self.MWs[i] if self.MWs is not None else None
        data['P_c'] = self.P_cs[i] if self.P_cs is not None else None
        data['V_c'] = self.V_cs[i] if self.V_cs is not None else None
        data['Z_c'] = self.Z_cs[i] if self.Z_cs is not None else None
        data['T_c'] = self.T_cs[i] if self.T_cs is not None else None
        data['w'] = self.ws[i] if self.ws is not None else None
        return data

    def set_point_input(self, i: int, **kwargs):
        """

        :param i: index of point
        """
        IdealMixture.set_point_input(self, i, **kwargs)
        self.MWs[i] = kwargs.pop('MW')
        self.P_cs[i] = kwargs.pop('P_c')
        self.T_cs[i] = kwargs.pop('T_c')
        self.V_cs[i] = kwargs.pop('V_c')
        self.Z_cs[i] = kwargs.pop('Z_c')
        self.ws[i] = kwargs.pop('w')

    def check_params(self):
        IdealMixture.check_params(self)
        if len(set(self.T_cs)) > 1:
            assert len(set(self.cas_numbers)) > 1, 'Errors anticipated when cas numbers equal if other props are not'