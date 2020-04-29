from abc import abstractmethod
import numpy as np

from pysurf.colt import PluginBase
from pysurf.spp import SurfacePointProvider
from pysurf.utils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.qctools.converter import Converter, time_converter

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

    def __init__(self, spp_inp, sampling, nstates, logger=None):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        self.nstates = nstates

        if logger is None:
            self.logger = get_logger('prop.log', 'propagator')
            info = {'spp inputfile': spp_inp}
            self.logger.header('PROPAGATION', info)
        else:
            self.logger = logger

        self.init = sampling.get_condition(0)
        # setup SPP
        self.spp = SurfacePointProvider(spp_inp, self.properties, sampling.natoms, nstates-1, sampling.atomids)
        
        if exists_and_isfile('prop.db'):
            self.db = DynDB.from_dynamics('prop.db')
            if len(self.db['crd']) > 0:
                self.restart = True
            else:
                self.restart = False
        else:
            if self.init is None:
                self.logger.error('No initial condition or restart file for the Propagator')
            else:
                self.db = DynDB.create_db('prop.db', sampling, nstates, self.properties)
            self.restart = False
       
        if self.restart is True:
            self.masses =  np.array(self.db['masses'])
        else:
            self.masses = sampling.masses

        self.t_converter = time_converter.get_converter(tin='au', tout='fs')

    def log_step(self, time, state, dE, ekin, epot, etot, infotext=None):
        self.logger.info(f"Time: {self.t_converter(time):6.2f} fs{'energy diff.:':>57} {dE:6.5f}")
        if infotext is not None: self.logger.info(infotext)
        self.logger.info(f"    {'curr. state':20}: {state:20}")
        self.logger.info(f"    {'ekin, epot, etot':20}: {ekin:10.6f}, {epot:10.6f}, {etot:10.6f}\n")
