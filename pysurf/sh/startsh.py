import logging
import os
import numpy as np

from pysurf.utils.constants import fs2au
from pysurf.utils.strutils import split_str
from pysurf.sh.model import landau_zener_surfacehopping
from pysurf.wigner import WignerSampling

class StartSh():
    def __init__(self, config):
        self.logger = logging.getLogger('pysurf')

        # parse SURFACE HOPPING section
        if 'SURFACE HOPPING' in config.keys():
            try:
                self.restart = config.getboolean('SURFACE HOPPING', restart)
            except:
                self.logger.info('No information on restart, assume no restart!')
                self.restart = False

            # get number of steps
            try:
                self.nsteps = config.getint('SURFACE HOPPING', 'number of steps')
            except:
                self.logger.error('Could not get number of steps')
                exit()
            # get random seed for random number generator of hopping algorithm
            # only important for developpers and debugging
            try:
                self.seed = config.getint('SURFACE HOPPING', 'random seed')
            except:
                self.seed = None

            # get stepsize
            try:
                self.dt = config.getfloat('SURFACE HOPPING', 'stepsize')
                self.dt = self.dt * fs2au
            except:
                self.logger.info('Could not find stepsize. Assume 0.5 fs')
                self.dt = 0.5 * fs2au
            
            # get starting state
            try:
                self.initstate = config.getint('SURFACE HOPPING', 'initial state')
            except:
                if self.restart is False:
                    self.logger.error('No initial state is provided')
                    exit()

            # get number of trajectories
            try:
                self.ntrajs = config.getint('SURFACE HOPPING', 'number of trajectories')
            except:
                self.logger.info('Could not find number of trajectories')

            # get inputfile for SPP
            try:
                self.sppinput = config.get('SURFACE HOPPING', 'inputfile for SPP')
                if os.path.isfile(self.sppinput):
                    self.sppinput = os.path.abspath(self.sppinput)
                else:
                    self.logger.error('Could not find inputfile for SPP')
            except:
                self.logger.info('SPP inputfile not given')
                exit()

            # get hopping algorithm
            try:
                self.hopalg = config.get('SURFACE HOPPING', 'hopping algorithm')
            except:
                self.logger.info('No hopping algorithm provided, assume landau zener')
                self.hopalg = 'landau zener'

        else:
            self.logger.error('No section SURFACE HOPPING in inputfile!')
            exit()

        # parse INITIAL CONDITIONS section
        if ('INITIAL CONDITIONS' in config.keys()) and not(self.restart):
            self.init = False
            self.molden = False
            self.freqs = False
            try:
                self.sampling = config.get('INITIAL CONDITIONS', 'sampling')
            except:
                self.sampling = 'wigner'
                self.logger.info('No sampling information found in INITIAL CONDITIONS! '
                                 + 'Assume Wigner sampling!')
            # molden input for initial conditions
            try:
                self.molden = config.get('INITIAL CONDITIONS', 'from molden file')
                if os.path.isfile(self.molden):
                    self.molden = os.path.abspath(self.molden)
                else:
                    self.logger.error('Molden file path not correct for INITIAL CONDITIONS')
                self.init = True
            except:
                pass

            #get frequencies for initial conditions
            try:
                vfloat = np.vectorize(float)
                self.freqs = vfloat(split_str(config.get('INITIAL CONDITIONS', 'frequencies')))
                self.init = True
            except:
                pass

        elif not('INITIAL CONDITIONS' in config.keys()) and not(self.restart):
            self.logger.error('No INITIAL CONDITIONS section in inputfile!')
        
        if self.restart is False:
            initconds = self.get_init()
            self.start_sh(initconds)
    
    def get_init(self, filename='init_conds.db'):
        if self.init is False:
            self.logger.error('Not enough information provided for the initial sampling!')
            exit()
        if self.sampling == 'wigner':
            if self.molden is not False:
                sampling = WignerSampling.from_molden(self.molden)
                conditions = sampling.create_initial_conditions(filename, self.ntrajs)
            elif self.freqs is not False:
                sampling = WignerSampling.from_freqs(self.freqs)
                conditions = sampling.create_initial_conditions(filename, self.ntrajs, model=True)
            else:
                self.logger.error('No source for the sampling given, e.g. from molden file!')
        else:
            self.logger.error('Sampling method not known!')
            exit()
        return conditions

    def start_sh(self, initconds):
        for i in range(self.ntrajs):
            dirname = 'traj.' + str(i)
            if os.path.isdir(dirname):
                continue
            os.mkdir(dirname)
            os.chdir(dirname)
            print('Johannes')
            print(initconds.get_condition(i))

            if self.hopalg == 'landau zener':
                landau_zener_surfacehopping(initconds.get_condition(i),
                                            self.initstate,
                                            self.nsteps,
                                            self.seed,
                                            self.sppinput,
                                            self.dt)
            else:
                self.logger.error('Hopping algorithm not implemented!')
                exit()
            os.chdir('../')
            self.logger.info('Finished trajectory ' + str(i))

