import os
import numpy as np
from shutil import copy2

from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.utils.osutils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.spp import SurfacePointProvider
from pysurf.sh import RunTrajectory

class SetupPropagation(Colt):
    _questions = """
    # Number of trajectories for the propagation
    n_traj = 100 :: int

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file

    #Filepath for the inputfile of the Propagation
    prop = prop.inp

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy_db = yes :: str :: [yes, no]

    #[copy_db(yes)]
    #db_file = none :: str
    """
    

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        self.config = config
        self.propagationfolder = 'prop'
        self.trajfolder = 'traj_'
        self.logger = get_logger('setup_propagation.log', 'setup_propagation')
        
        if not(os.path.isdir(self.propagationfolder)):
            os.mkdir(self.propagationfolder)
        
        # Open DB of initial conditions once, so that it is available
        self.sampling = Sampling.from_db(self.config['sampling_db'])

        RunTrajectory.generate_input(config['prop'], config=config['prop'])
        SurfacePointProvider.generate_input(config['spp'], config=config['spp'])

        for i in range(self.config['n_traj']):
            self._setup_trajectory_folder(i)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    @classmethod
    def from_inputfile(cls, inputfile):
        quests = cls.generate_questions(config=inputfile)
        config = quests.ask(inputfile)
        return cls.from_config(config)

    def _setup_trajectory_folder(self, number):
        foldername = os.path.join(self.propagationfolder, self.trajfolder + '{:04d}'.format(number))
        if os.path.isdir(foldername):
            self.logger.info('Folder {0} already exists. Skipping trajectory'.format(foldername))
            return

        os.mkdir(foldername)

        copy2(self.config['prop'], foldername)
        copy2(self.config['spp'], foldername)

        if self.config['copy_db'] == 'yes':
            copy2(self.config['copy_db']['db_file'], foldername)

        initname = os.path.join(foldername, 'init.db'.format(number))
        self.sampling.export_condition(initname, number)


if __name__=="__main__":
    SetupPropagation.from_inputfile('setup_propagation.inp')
