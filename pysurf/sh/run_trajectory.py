import logging
import os
import numpy as np
import logging

from collections import namedtuple

from pysurf.colt import Colt
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.utils.constants import fs2au
from pysurf.utils.strutils import split_str
from pysurf.sh.landau_zener import landau_zener_surfacehopping
from pysurf.sampling import Sampling

class RunTrajectory(Colt):
    _questions = """
    # Total propagation time in fs
    propagation time in fs = 100 :: float                                                            
    
    # Time step in fs for the propagation                                                            
    time step in fs = 0.5 :: float                                                                   
     
    # File with initial condition                                                                    
    initial condition = init.db :: str                                                          
    
    #Filepath for the inputfile of the Surface Point Provider
    spp inputfile = spp.inp :: str    
    """
    def __init__(self, config):
        self.logger = logging.getLogger('runtrajectory')
        self.inputfile = inputfile

        quests = self.generate_questions("SURFACE HOPPING", config=inputfile)
        self.config = quests.ask(inputfile)

        self.initcond = Sampling.from_db(self.config['initial condition'])

        self.logger.info('Start propagation of trajectory')
        landau_zener_surfacehopping(self.initcond.get_condition(1),
                                                self.initstate,
                                                self.nsteps,
                                                self.config['spp inputfile'],
                                                self.config['time step in fs']*fs2au)
    
    @classmethod
    def from_inputfile(cls,inputfile):
        quests = cls.generate_questions(config=inputfile)
        config = quests.ask(inputfile)
        return cls(config)
