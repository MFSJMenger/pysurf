from abc import abstractmethod
import os
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
        """ propagation routine"""
        pass

    def run(self, nsteps, dt, *args, **kwargs):
        self.output_header()
        self.start_time = time.perf_counter()
        self._run(nsteps, dt, *args, **kwargs)

    def get_runtime(self):
        return (time.perf_counter() - self.start_time)

    def __init__(self, spp_inp, sampling, nstates, properties=None restart=True, logger=None):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        self.nstates = nstates
        self.start_time = time.perf_counter()
        self.sampling = sampling


        if logger is None:
            self.logger = get_logger('prop.log', 'propagator')
            info = {'spp inputfile': spp_inp}
            self.logger.header('PROPAGATION', info)
        else:
            self.logger = logger
        
        for prop in properties:
            if prop not in self.properties:
                self.properties += [prop]

        self.init = sampling.get_condition(0)
        # setup SPP
        if sampling.model is False:
            self.spp = SurfacePointProvider(spp_inp, self.properties, nstates, sampling.natoms, sampling.atomids)
        else:
            self.spp = SurfacePointProvider(spp_inp, self.properties, nstates, sampling.nmodes)
        

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
                self.create_new_db()
            self.restart = False

        if restart is False: self.restart = False

        if self.restart is True:
            self.masses =  np.array(self.db['masses'])
            mode_output = "a"
        else:
            self.masses = sampling.masses
            mode_output = "w"

        self.t_converter = time_converter.get_converter(tin='au', tout='fs')
        self.output = get_logger('prop.out', 'propagator_output', mode=mode_output)

    def create_new_db(self):
        name = 'prop.db'
        if exists_and_isfile(name): os.remove(name)
        self.db = DynDB.create_db(name, self.sampling, self.nstates, self.properties)

    def output_header(self):
        self.output.info('#'+('='*101))
        self.output.info(f"#{'Step':^9}|{'Time':^12}|{'State':^7}|{'Energy':^47}|{'Gradient':^11}|{'Runtime':^9}|")
        self.output.info(f"#{' ':9}|{' ':^12}|{' ':^7}|{'kin':^11}|{'pot':^11}|{'tot':^11}|{'diff':^11}|{'RMS':^11}|{' ':^9}|")
        self.output.info(f"#{' ':9}|{'[fs] ':^12}|{' ':^7}|{'[au]':^47}|{'[au]':^11}|{'[sec]':^9}|")
        self.output.info('#' + ('='*101) )

    def output_step(self, step, time, state, ekin, epot, etot, dE, grad=None):
        self.output.info(f"{step:>10}{self.t_converter(time):>13.2f}{state:>8} {ekin:>11.6f} {epot:>11.6f} {etot:>11.6f} {dE:>11.6f}{' ':12}{round(self.get_runtime(), 1):>10}")
