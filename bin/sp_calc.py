from pysurf.logger import get_logger 
from pysurf.sampling import Sampling

from pysurf.utils import exists_and_isfile
from pysurf.spp import SurfacePointProvider
from pysurf.colt import Colt


class SinglePointCalculation(Colt):

    _questions = """
    # Number of excited states
    nstates = 3 :: int

    # Database containing the geometry
    init_db = init.db :: existing_file

    # Properties that should be calculated
    properties = [energy, gradient] :: list

    # SPP inputfile
    spp = spp.inp :: file
    """

    def __init__(self, config):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        logger = get_logger('sp_calc.log', 'sp_calc')
        #
        sampling = Sampling.from_db(config['init_db'])

        if not(exists_and_isfile(config['spp'])):
            SurfacePointProvider.generate_input(config['spp'])
   
        spp = SurfacePointProvider(config['spp'], 
                                  config['properties'],
                                  sampling.natoms,
                                  config['nstates'],
                                  sampling.atomids)

        crd = sampling.get_condition(1).crd
        
        res = spp.request(crd, config['properties'])
        print(res)
            
    @classmethod
    def from_config(cls, config):
        return cls(config)


if __name__=="__main__":
    SinglePointCalculation.from_commandline()
