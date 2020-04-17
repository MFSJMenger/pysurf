import os
from shutil import copy2 as copy
#
from pysurf.logger import get_logger
from pysurf.sampling.sampling import SamplingBase
from pysurf.setup import SetupBase


class SetupSampling(SetupBase):

    folder = 'sampling'
    subfolder = 'init'

    _questions = """
    # Number of Calculations
    njobs = 100 :: int

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: existing_file

    file = ex.xyz :: existing_file

    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('setup_sampling.log', 'setup_sampling')
        SetupBase.__init__(self, logger)
        # Open DB of initial conditions once, so that it is available
        initconds = SamplingBase.from_db(config['sampling_db'])
        #
        self.setup_folders(range(config['njobs']), config, initconds)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, initconds):
        copy(config['spp'], foldername)
        copy(config['file'], foldername)
        # 
        xyzfile = os.path.join(foldername, 'geom.xyz')
        initconds.write_xyz(xyzfile, number)


if __name__=="__main__":
    SetupSampling.from_commandline()
