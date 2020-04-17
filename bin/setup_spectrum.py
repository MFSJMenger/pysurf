import os
from shutil import copy2 as copy

from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.setup import SetupBase


class SetupSpectrum(SetupBase):

    folder = 'spectrum'
    subfolder = 'condition'

    _questions = """
    # Number of trajectories for the propagation
    n_traj = 100 :: int

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file
    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('setup_spectrum.log', 'setup_spectrum')
        SetupBase.__init__(self, logger)
        #
        sampling = Sampling.from_db(config['sampling_db'])
        #
        self.setup_folders(range(config['n_traj']), config, sampling, skip_existing=True)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling):
        copy(inputfile, foldername)
        copy(config['spp'], foldername)
        initname = os.path.join(foldername, 'init_{:04d}.db'.format(number))
        sampling.export_condition(initname, number)


if __name__=="__main__":
    SetupSpectrum.from_commandline()
