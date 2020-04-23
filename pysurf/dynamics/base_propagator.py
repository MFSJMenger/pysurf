from abc import abstractmethod
import numpy as np

from pysurf.colt import PluginBase
from pysurf.random.shrandom import RandomNumberGeneratorNP
from pysurf.spp import SurfacePointProvider
from pysurf.utils import exists_and_isfile
from pysurf.logger import get_logger

from .dyn_db import DynDB

class PropagatorFactory(PluginBase):
    _plugins_storage = '_methods'
    _is_plugin_factory = True

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                            for name, method in cls._methods.items()})

class PropagatorBase(PropagatorFactory):
    _questions = 'inherited'
    subquestions: 'inherited'

    _register_plugin = False
    #Properties have to be defined for each Propagator
    properties = ['energy', 'gradient']

    @abstractmethod
    def run(self, nsteps, dt, *args, **kwargs):
        """propagation routine"""

    def __init__(self, spp_inp, natoms, nstates, atomids, mass, init_cond, logger=None):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        if logger is None:
            self.logger = get_logger('prop.log', 'propagator')
        else:
            self.logger = logger

        self.init = init_cond
        # setup SPP
        self.spp = SurfacePointProvider(spp_inp, self.properties, natoms, nstates, atomids)
        
        if exists_and_isfile('prop.db'):
            self.db = DynDB.from_dynamics('prop.db')
            if self.db.len > 0:
                self.restart = True
            else:
                self.restart = False
        else:
            if init_cond is None:
                self.logger.error('No initial condition or restart file for the Propagator')
            else:
                self.db = DynDB.create_db('prop.db', natoms, nstates, model=False, nmodes=0)
            self.restart = False
       
        if self.restart is True:
            self.mass = np.array(self.db['mass'])
        else:
            self.mass = np.array(mass)
            self.db.add_mass(self.mass)

    def init_random(self, seed=16661):
        self._random = RandomNumberGeneratorNP(seed)

    def random(self):
        return self._random.get()
