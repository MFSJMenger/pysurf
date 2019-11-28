import os
from copy import deepcopy
from collections import namedtuple, OrderedDict
#
import numpy as np
import numpy.random as random
#
from pysurf.molden import MoldenParser
#
from pysurf.logger import get_logger
from pysurf.utils.osutils import exists_and_isfile
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import DatabaseTools
from .base_sampling import InitialCondition
from .sampling import SamplingBase
#
from pysurf.molecule.molecule import Molecule
from .normalmodes import Mode 
#from .wigner import Wigner


class InitialConditions(SamplingBase):
    _questions = 'inherited'
    subquestions : 'inherited'
    _register_plugin = False

    @classmethod
    def _generate_subquestions(cls, questions):
        """ This class will not be inherited """
        # Update _questions from Sampling by adding an additional question
        questions.add_questions_to_block("""
            # State on which trajectories start
            initial state = 0 :: int 
            """)


    def __init__(self, config, logger = None):
        super().__init__(config, logger)
        self.condition = InitialCondition
    
    def get_condition(self, idx):
        if idx >= self.nconditions:
            return None
        crd = self._db.get('crd', idx)
        veloc = self._db.get('veloc', idx)
        #state = self._db.get('state', idx)
        return InitialCondition(crd, veloc, None)


    @property
    def _settings(self):
        settings = super()._settings
        
        # Add initial state
        settings['variables']['state'] = DBVariable(np.double, ('frame', 'one'))
        
        # Add velocities
        if self.model is False:
            settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'natoms', 'three'))
        if self.model is True:
            settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'nmodes'))
        
        return settings

        
    @staticmethod
    def _write_condition(db, cond):
        db.append('crd', cond.crd)
        db.append('veloc', cond.veloc)
        db.append('state', cond.state)
        db.increase

