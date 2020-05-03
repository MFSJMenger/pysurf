import os
import numpy as np

from collections import namedtuple

from pysurf.logger import Logger, get_logger
from pysurf.colt import Colt
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.utils.constants import fs2au
from pysurf.utils.strutils import split_str
from pysurf.spp import SurfacePointProvider
from .base_propagator import PropagatorFactory
from pysurf.sampling import Sampling


class RunTrajectory(Colt):
    _questions = """
    # Total propagation time in fs
    time_final [fs] = 100 :: float                                                            
    
    # Time step in fs for the propagation                                                            
    timestep [fs] = 0.5 :: float                                                                   
     
    # File with initial condition                                                                    
    initial condition = init.db :: str                                                          
   
    # Number of total states
    n_states = 2 :: int

    method = LandauZener :: str

    #Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: str    
    
    restart = True :: bool
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                 for name, method in PropagatorFactory._methods.items()})

    def __init__(self, config):
        self.logger = get_logger('prop.log', 'prop')
        self.logger.header('PROPAGATION', config)
        sampling = Sampling.from_db(config['initial condition'])

        self.nsteps = int(np.ceil(config['time_final [fs]'] / config['timestep [fs]']))
        natoms = len(sampling.atomids)

        self.logger.info('Start propagation of trajectory')

        #get propagator
        propagator = PropagatorFactory._methods[config['method'].value](config['spp'], 
                                                                        sampling,
                                                                        config['n_states'],
                                                                        restart=config['restart'],
                                                                        logger=self.logger)
        propagator.run(self.nsteps, config['timestep [fs]']*fs2au)
    
    @classmethod
    def from_inputfile(cls, inputfile):
        quests = cls.generate_questions(config=inputfile)
        config = quests.ask(inputfile)
        return cls(config)
