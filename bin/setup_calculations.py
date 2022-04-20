import os
from shutil import copy2 as copy
#
from pysurf.setup import SetupBase
from pysurf.logger import get_logger
from pysurf.sampling import Sampling


class SetupComputation(SetupBase):

    folder = 'computations'
    subfolder = 'sp'

    _user_input = """
    # Database containing all the initial conditions
    sampling_db = :: existing_file
    """


    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        self.logger = get_logger('setup.log', 'setup')
        super().__init__(self.logger)
        sampling = Sampling.from_db(config['sampling_db'], logger=self.logger)
        #
        self.setup_folders(range(sampling.nconditions), config, sampling)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling):
        #name of new database
        qchemfile = os.path.join(foldername, 'init.xyz')
        #get info from old db and adjust
        condition = sampling.get_condition(number)
        with open(qchemfile, 'w') as fh:
            fh.write('$molecule\n')
            for idx, crd in zip(sampling.atomids, condition.crd):
                fh.write("%d %12.8f  %12.8f  %12.8f\n" % (idx, crd[0], crd[1], crd[2]))
            fh.write('$end\n')


if __name__=="__main__":
    SetupComputation.from_commandline()
