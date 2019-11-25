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
from .wigner import WignerSampling
from .base_sampling import InitialCondition
from .sampling import Sampling
#
from pysurf.molecule.molecule import Molecule
from .normalmodes import Mode 
#
from pysurf.colt import AskQuestions
from pysurf.colt import Colt
#  used to save just InitialConditions (coordinates, velocities) for a given molecule


class InitialConditions(Sampling, Colt):
    _methods = OrderedDict({'wigner': WignerSampling})

    # Update _questions from Sampling by adding an additional question
    _questionsInitCond = """
    # State on which trajectories start
    initial state = 0 :: int 
    """
    _questions = Sampling._questions.rstrip('"') + _questionsInitCond.lstrip('"')

#    _questions_old = """
#    # database file where the initial conditions are saved or from which the initial conditions
#    # are taken if the file already exists.
#    outputfile = initconds.db
#
#
#    # Describes which sampling algorithm is used to generate the initial conditions.
#    # The default is wigner.
#    sampling = wigner :: str :: [wigner]
#
#    # State whether the system you are describing is a model system
#    model = False :: bool
#
#    # Initial state of the trajectories
#    initial state = 0 :: int
#
#    # Number of initial conditions that have to be created.
#    # The default value is 100.
#    number of initial conditions = 100 :: int
#    """

    def __init__(self, config, logger = None):
        super(InitialConditions, self).__init__(config, logger)
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

