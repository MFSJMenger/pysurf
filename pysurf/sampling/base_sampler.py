from collections import namedtuple
from abc import abstractmethod

from colt import Plugin

CrdCondition = namedtuple("CrdCondition",['crd'])
DynCondition = namedtuple("DynCondition", ['crd', 'veloc', 'state'])


class SamplerFactory(Plugin):
    """ Factory for samplers, which provide only coordinates. It is also the underlying class for
        DynSamplingFactory, which is an extension for samplers with velocities and initial states.
    """
    _plugins_storage = '_methods'
    _is_plugin_factory = True

    condition = CrdCondition
    _n_conditions = None

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases("method", {name: method.colt_user_input
                                 for name, method in cls._methods.items()})

    def get_number_of_conditions(self, n_max):            
        if self._n_conditions is None:
            return n_max
        if self._n_conditions > n_max:
            return n_max
        return self._n_conditions
            

class CrdSamplerBase(SamplerFactory):
    """Basic sampler class"""

    extend_user_input: "inherited"
    #
    _register_plugin = False
    _user_input = "inherited"

    @classmethod
    @abstractmethod
    def from_config(cls, config, start=0):
        pass

class DynSamplerFactory(CrdSamplerBase):
    """Factory to store intial condition samplers"""
    # setup questions
    _colt_user_input = 'inherited'
    extend_user_input: 'inherited'
    # setup plugin
    _plugins_storage = '_methods'
    _is_plugin_specialisation = True
    _is_plugin_factory = True
    _register_plugin = False
    #
    condition = DynCondition

    @classmethod
    def _extend_user_input(cls, questions):
        """ This class will not be inherited """
        # Update _user_input from Sampling by adding an additional question
        questions.add_questions_to_block("""
            # State on which trajectories start
            initial state = 0 :: int
        """)
        questions.generate_cases("method", {name: method.colt_user_input
                                 for name, method in cls._methods.items()})


class DynSamplerBase(DynSamplerFactory):
    """Base Class for Dynamics Conditions Sampler"""

    _register_plugin = False
    _colt_user_input = 'inherited'

    @classmethod
    @abstractmethod
    def from_config(cls, config, start=0):
        pass
