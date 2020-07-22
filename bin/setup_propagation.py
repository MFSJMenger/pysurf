import os
from shutil import copy2 as copy
#
from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.setup import SetupBase
from pysurf.spp import SurfacePointProvider
from pysurf.dynamics import RunTrajectory
from pysurf.utils import exists_and_isfile


class SetupPropagation(SetupBase):

    folder = 'prop'
    subfolder = 'traj'

    _questions = """
    # Number of trajectories for the propagation
    n_traj = -1 :: int

    # Database containing all the initial conditions
    sampling_db = :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file
    
    # initial excited state for the trajectory
    initial state =  :: int

    #Filepath for the inputfile of the Propagation
    prop = prop.inp :: file

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy_db = none :: str
    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('setup_propagation.log', 'setup_propagation')
        SetupBase.__init__(self, logger)


        # Open DB of initial conditions once, so that it is available
        sampling = Sampling.from_db(config['sampling_db'])

        #Make sure that inputfile for the SPP exists and is complete
        
        if exists_and_isfile(config['spp']): lconfig = config['spp']
        else: lconfig = None
        SurfacePointProvider.generate_input(config['spp'], config=lconfig)

        #Make sure that inputfile for RunTrajectory exists and is complete
        if exists_and_isfile(config['prop']): lconfig = config['prop']
        else: lconfig = None
        RunTrajectory.generate_input(config['prop'], config=lconfig)

        #
        if config['n_traj'] == -1:
            ntraj = len(sampling._db)
        else:
            ntraj = config['n_traj']
        if sampling.nconditions < ntraj:
            logger.error(f"Too few initial conditions in {config['sampling_db']}")

        self.setup_folders(range(ntraj), config, sampling)


    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling):
        copy(config['prop'], foldername)
        copy(config['spp'], foldername)

        if config['copy_db'] != 'none':
            copy(config['copy_db'], foldername)

        initname = os.path.join(foldername, 'init.db')
        #setup new database 
        new_sampling = Sampling.create_db(initname, sampling.info['variables'], sampling.info['dimensions'], sampling.system, sampling.modes, model=sampling.model, sp=True)
        #copy condition to new db
        condition = sampling.get_condition(number)
        new_sampling.write_condition(condition, 0)
        new_sampling.set('currstate', config['initial state'], 0)


if __name__=="__main__":
    SetupPropagation.from_commandline()
