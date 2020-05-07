from collections import namedtuple
from abc import abstractmethod

from pysurf.colt import PluginBase

CrdCondition = namedtuple("CrdCondition",['crd'])
DynCondition = namedtuple("DynCondition", ['crd', 'veloc', 'state'])


class SamplerFactory(PluginBase):
    """ Factory for samplers, which provide only coordinates. It is also the underlying class for 
        DynSamplingFactory, which is an extension for samplers with velocities and initial states.
    """
    _plugins_storage = '_methods'
    _is_plugin_factory = True

    condition = CrdCondition
    
    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                 for name, method in cls._methods.items()})


class CrdSamplerBase(SamplerFactory):
    """Basic sampler class"""

    subquestions: "inherited"
    #
    _register_plugin = False
    _questions = "inherited"

    @classmethod
    @abstractmethod
    def from_config(cls, config, start=0):
        pass

class DynSamplerFactory(CrdSamplerBase):
    """Factory to store intial condition samplers"""
    # setup questions
    questions = 'inherited'
    subquestions: 'inherited'
    # setup plugin
    _plugins_storage = '_methods'
    _is_plugin_specialisation = True
    _is_plugin_factory = True
    _register_plugin = False
    #
    condition = DynCondition

    @classmethod
    def _generate_subquestions(cls, questions):
        """ This class will not be inherited """
        # Update _questions from Sampling by adding an additional question
        questions.add_questions_to_block("""
            # State on which trajectories start
            initial state = 0 :: int
        """)
        questions.generate_cases("method", {name: method.questions
                                 for name, method in cls._methods.items()})


class DynSamplerBase(DynSamplerFactory):
    """Base Class for Dynamics Conditions Sampler"""

    _register_plugin = False
    questions = 'inherited'

    @classmethod
    @abstractmethod
    def from_config(cls, config, start=0):
        pass
