import os
from shutil import copy2 as copy
#
from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.setup import SetupBase


class SetupPropagation(SetupBase):

    folder = 'prop'
    subfolder = 'traj'

    _questions = """
    # Total propagation time in fs
    propagation time in fs = 100 :: float

    # Time step in fs for the propagation
    time step in fs = 0.5 :: float
    
    # Number of trajectories for the propagation
    n_traj = 100 :: int

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy_db = yes :: str :: [yes, no]

    [copy_db(yes)]
    db_file = none :: str
    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('setup_propagation.log', 'setup_propagation')
        SetupBase.__init__(self, logger)

        n_traj = config['n_traj']

        # Open DB of initial conditions once, so that it is available
        sampling = Sampling.from_db(config['sampling_db'])
        #
        self.setup_folders(range(n_traj), config, sampling)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling):
        copy(self.inputfile, foldername)
        copy(config['spp'], foldername)

        if config['copy_db'] == 'yes':
            copy(config['copy_db']['db_file'], foldername)

        initname = os.path.join(foldername, 'init.db')
        sampling.export_condition(initname, number)


if __name__=="__main__":
    SetupPropagation.from_commandline()
