import os
import numpy as np
from shutil import copy2

from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.utils.osutils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.sampling.sampling import Sampling

class SetupSampling(Colt):

    _questions = """
    # Number of Calculations
    njobs = 100 :: int

    # Name of the folder
    name = prop :: str

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
        self.propagationfolder = config['name']
        self.initcondname = config['sampling_db']
        self.sppinp = config['spp']
        self.refgeom = config['file']
        self.logger = get_logger('setup_sampling.log', 'setup_sampling')
        
        if not(os.path.isdir(self.propagationfolder)):
            os.mkdir(self.propagationfolder)
        
        # Open DB of initial conditions once, so that it is available
        self.initconds = Sampling.from_db(config['sampling_db'])

        for i in range(config['njobs']):
            self._setup_trajectory_folder(i)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def _setup_trajectory_folder(self, number):
        foldername = os.path.join(self.propagationfolder, 'init_%06d' % number)
        xyzfile = os.path.join(foldername, 'geom.xyz')
        if os.path.isdir(foldername):
            self.logger.info('Folder {0} already exists. Skipping folder'.format(foldername))
            return

        os.mkdir(foldername)

        copy2(self.sppinp, foldername)
        copy2(self.refgeom, foldername)
        cond = self.initconds.get_condition(number)
        
        self.initconds.write_xyz(xyzfile, number)


if __name__=="__main__":
    SetupSampling.from_commandline()
