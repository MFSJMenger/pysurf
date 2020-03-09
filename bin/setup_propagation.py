import os
import numpy as np
import click
from shutil import copy2

from pysurf.colt import Colt
from pysurf.utils.osutils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.sampling.initialconditions import InitialConditions

class SetupPropagation(Colt):
    _questions = """
    # Total propagation time in fs
    propagation time in fs = 100 :: float

    # Time step in fs for the propagation
    time step in fs = 0.5 :: float
    
    # Number of trajectories for the propagation
    number of trajectories = 100 :: int

    # Database containing all the initial conditions
    database with the initial conditions = sampling.db :: str

    # Filepath for the inputfile of the Surface Point Provider
    spp inputfile = spp.inp :: str

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy database = yes :: str :: [yes, no]

    [copy database(yes)]
    database filename = none :: str
    """
    

    def __init__(self, inputfile):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        self.inputfile = inputfile
        self.propagationfolder = 'prop'
        self.trajfolder = 'traj.'
        self.initcondname = 'initcond.db'
        self.logger = get_logger('setup_propagation.log', 'setup_propagation')
        
        quests = self.generate_questions("SURFACE HOPPING", config=inputfile)
        self.config = quests.ask(inputfile)

        if not(os.path.isdir(self.propagationfolder)):
            os.mkdir(self.propagationfolder)
        
        # Open DB of initial conditions once, so that it is available
        if exists_and_isfile(self.config['database with the initial conditions']):
            self.initconds = InitialConditions.from_db(self.config['database with the initial conditions'])
        else:
            self.logger.error('Database with the initial conditions not found')

        for i in range(self.config['number of trajectories']):
            self._setup_trajectory_folder(i)

    def _setup_trajectory_folder(self, number):
        foldername = self.propagationfolder + '/' + self.trajfolder + str(number)
        if os.path.isdir(foldername):
            self.logger.info('Folder {0} already exists. Skipping trajectory'.format(foldername))
            return

        os.mkdir(foldername)

        copy2(self.inputfile, foldername)
        if exists_and_isfile(self.config['spp inputfile']):
            copy2(self.config['spp inputfile'], foldername)
        else:
            self.logger.error('spp inputfile not found!')

        if self.config['copy database'] == 'yes':
            if exists_and_isfile(self.config['copy database']['database filename']):
                copy2(self.config['copy database']['database filename'], foldername)
            else:
                self.logger.error('Database file should be copied but cannot be found!')

        self._setup_initial_condition_for_trajectory(number)

    def _setup_initial_condition_for_trajectory(self, number):
        foldername = self.propagationfolder + '/' + self.trajfolder + str(number)
        cond = self.initconds.get_condition(number)
        initcondname = os.path.join(foldername,self.initcondname)
        self.initconds.export_condition(initcondname, number)


@click.command()
@click.option('-f', 'filename', default='propagation.inp')
def command_setup_propagation(filename):
    SetupPropagation(filename)

if __name__=="__main__":
    command_setup_propagation()
