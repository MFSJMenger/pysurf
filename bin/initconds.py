import os
import numpy as np
import click

from pysurf.colt import AskQuestions
from pysurf.utils.constants import fs2au
from pysurf.utils.osutils import exists_and_isfile
from pysurf.utils.strutils import split_str
from pysurf.initconds.wigner import WignerSampling
from pysurf.initconds.wigner import InitialConditions
from pysurf.logger import get_logger
from pysurf.utils.decorator import print_timelog, _total_times, timeit

class InitConds():
    def __init__(self, inputfile):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """

        question_string = """
        # database file where the initial conditions are saved or from which the initial conditions
        # are taken if the file already exists.
        outputfile = initconds.db

        # Input source for the normal modes and/or frequencies, which are used to generate the 
        # initial conditions.
        # Possible options are:
        # - molden
        # - frequencies
        from = none :: str :: [molden, frequencies]

        # Describes which sampling algorithm is used to generate the initial conditions.
        # The default is wigner.
        sampling = wigner

        # Number of initial conditions that have to be created.
        # The default value is 100.
        number of initial conditions = 100 :: int

        # If initial conditions are generated from a molden file, this subquestion asks for the 
        # molden file.
        [from(molden)]
        moldenfile = none

        # If initial conditions are generated from frequencies, here a list of frequencies has to
        # be given for the modes. The list can be in the python list format or just comma separated
        # values.
        [from(frequencies)]
        frequencies = none
        """

        self.logger = get_logger('initconds.log', 'initconds')

        quests = AskQuestions.from_string("INITIAL CONDITIONS", question_string, config=inputfile)
        
        # Read or create the inputfile and ask for missing options with the AskQuestion class
        self.config = quests.ask(inputfile)
#        self.config = quests.check_only(inputfile)
        quests.create_config_from_answers(inputfile)

        # create initconds
        if exists_and_isfile(self.config['outputfile']):
            self.logger.info('Get initial conditions from existing database: ' 
                             + str(self.config['outputfile']))
            self.initconds = InitialConditions.from_db(self.config['outputfile'])
            nconditions = self.config['number of initial conditions']
            if self.initconds.nconditions < nconditions:
                newconditions = nconditions - self.initconds.nconditions
                self.logger.info('Adding ' + str(newconditions) + ' new initial conditions...')
                self.initconds.add_initial_conditions(newconditions)
        else:
            self.logger.info('Create initial conditions and store them in '
                             + str(self.config['outputfile']))
            self.initconds = self._create_inits(self.config['outputfile'])

    def _create_inits(self, filename):
        ninits = self.config['number of initial conditions']
        if self.config['from'] == 'molden':
            moldenfile = self.config['from']['moldenfile']
            sampling = WignerSampling.from_molden(moldenfile)
            conditions = sampling.create_initial_conditions(filename, ninits)
        elif self.config['from'] == 'frequencies':
            freqs = self.config['from']['frequencies']
            freqs = split_str(freqs)
            vfloat = np.vectorize(float)
            freqs = vfloat(freqs)
            sampling = WignerSampling.from_freq(freqs)
            conditions = sampling.create_initial_conditions(filename, ninits, model=True)
        return conditions


def command_initconds(filename):
    InitConds(filename)

if __name__=="__main__":
    command_initconds("initconds.inp")
