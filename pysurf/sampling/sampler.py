from abc import abstractmethod
from collections import namedtuple
from copy import deepcopy
#
import numpy as np
#from pysurf.logger import get_logger
from pysurf.utils.osutils import exists_and_isfile
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.colt import Colt
from pysurf.logger import Logger, get_logger
#
from pysurf.molecule.molecule import Molecule
from pysurf.sampling.base_sampler import SamplerFactory
from .base_sampler import DynCondition, CrdCondition
from .normalmodes import Mode




class Sampler(Colt):
    """ Sampling is the header class for the sampling routines. It asks the main questions, selects
        the sampler and reads and writes the conditions to the sampling database.
    """
    _questions = """
    method = :: str

    #State whether the system is a model system
    model = False :: bool
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                 for name, method in SamplerFactory._methods.items()})

    @classmethod
    def from_config(cls, config):
        sampler = cls._get_sampler(config['method'])
        return cls(config, sampler)

    @classmethod
    def from_inputfile(cls, inputfile):
        # Generate the config
        quests = cls.generate_questions('SAMPLING', config=None)
        config = quests.check_only(inputfile)
        return cls.from_config(config)    



    def __init__(self, config, sampler, logger=None):
        """ Sampling always goes with a database, if not needed use Sampler class
        """

        self.config = config
        self.sampler = sampler

        if logger is None:
            self.logger = get_logger('sampler.log', 'sampler')
        else:
            self.logger = logger


    def get_condition(self, idx):
        return self.sampler.get_condition()

    def __iter__(self):
        self._start = 1  # skip equilibrium structure
        return self

    def __next__(self):
        cond = self.get_condition(self._start)
        if cond is not None:
            self._start += 1
            return cond
        raise StopIteration

    @staticmethod
    def _get_sampler(config):
        return SamplerFactory._methods[config.value].from_config(config)
    
    @property
    def equilibrium(self):
        return self.sampler.get_init()
