from abc import abstractmethod
import numpy as np
import time

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
    def _run(self, nsteps, db, *args, **kwargs):

    def run(self, nsteps, dt, *args, **kwargs):
        """propagation routine"""
#       start timer
#       ....
        self.output_header(self)
        self.start_time = time.perf_counter()
        self._run(self, nsteps, dt, *args, **kwargs)

    def get_runtime(self):
        return (time.perf_counter() - self.start_time)

    def __init__(self, spp_inp, sampling, nstates, logger=None):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        self.nstates = nstates
        self.start_time = time.perf_counter()

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
        self.output = get_logger('prop.out')

    def output_header(self):
        self.output.info('#'+('='*119))
        self.output.info(f"#{'Step':^6}|{'Time':^6}|{'State':^7}|{'Energy':^28}|{'Gradient':^10}|{'Runtime':^9}|")
        self.output.info(f"#{:6}|{:^6}|{:^7}|{'kin':^6}|{'pot':^6}|{'tot':^6}|{'diff':^6}|{'RMS:^10}|{:^9}|")
        self.output.info(f"#{:6}|{[fs] :^6}|{:^7}|{'[au]':^21}|{'[au]':^10}|{'[sec]':^9}|")
        self.output.info('#' + ('='*119) + '\n')

    def output_step(self, step, time, state, dE, ekin, epot, etot):
        self.output.info(f"{step:8}{self.t_converter(time):6.2f}{ekin:7.5f}{epot:7.5f}{etot:7.5f}{ediff:7.5f}{:11}{self.get_runtime():9}")
