import os
from shutil import copy2 as copy
#
from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.setup import SetupBase
from pysurf.utils import exists_and_isfile
from pysurf.spp import SurfacePointProvider

from sp_calc import SinglePointCalculation

class SetupSpectrum(SetupBase):

    folder = 'spectrum'
    subfolder = 'condition'

    _questions = """
    # Number of conditions
    n_cond = :: int

    # Number of states
    nstates = :: int
   
    #Properties that should be calculated
    properties = ['energy', 'fosc'] :: list

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file

    # Filepath for the inputfile of the Single Point Calculation
    sp_calc = sp_calc.inp :: file
    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('setup_spectrum.log', 'setup_spectrum')
        logger.header('SETUP SPECTRUM', config)
        SetupBase.__init__(self, logger)
        #
        logger.info(f"Opening sampling database {config['sampling_db']}")
        sampling = Sampling.from_db(config['sampling_db'], logger=logger)

        if not exists_and_isfile(config['spp']):
            presets="""
                use_db = no
                """
            logger.info(f"Setting up SPP inputfile: {config['spp']}")
            SurfacePointProvider.generate_input(config['spp'], config=None, presets=presets)
        else:
            logger.info(f"Using SPP inputfile as it is")
            
        if not exists_and_isfile(config['sp_calc']):
            presets=f"""
                properties = {config['properties']}
                nstates = {config['nstates']}
                init_db = init.db
                """
            logger.info(f"Setting up inputfile for the single point calculations")
            SinglePointCalculation.generate_input(config['sp_calc'], config=None, presets=presets)
        else:
            logger.info(f"Using inputfile for the single point calculations as it is")

        logger.info("Starting to prepare the folders...")
        self.setup_folders(range(config['n_cond']), config, sampling)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling):
        copy(config['spp'], foldername)
        copy(config['sp_calc'], foldername)

        #name of new database
        initname = os.path.join(foldername, 'init.db')
        #get info from old db and adjust
        variables = sampling.info['variables']
        variables += config['properties']
        dimensions = sampling.info['dimensions']
        dimensions['nstates'] = config['nstates']
        dimensions['nactive'] = config['nstates']
        #setup new database 
        new_sampling = Sampling.create_db(initname, variables, dimensions, sampling.molecule, sampling.modes, model=sampling.model, sp=True)
        #copy condition to new db
        condition = sampling.get_condition(number)
        new_sampling.write_condition(condition, 0)

if __name__=="__main__":
    SetupSpectrum.from_commandline()
